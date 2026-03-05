#!/usr/bin/env python 

import os
import re
import tifffile
import string
import argparse
import random
import time
import numpy as np
import logging
import cv2 as cv
import copy
import gc
from tqdm import tqdm
from basicpy import BaSiC
from matplotlib import pyplot as plt
from scipy import ndimage as ndi
from tglow.io.compound_image_provider import CompoundImageProvider
from tglow.utils.tglow_utils import float_to_16bit_unint
from skimage.filters import threshold_otsu, gaussian, threshold_multiotsu, threshold_sauvola, rank
from skimage.morphology import disk
from skimage.transform import downscale_local_mean, resize, resize_local_mean
from numpy.polynomial import polynomial as P
from sklearn.linear_model import RidgeCV
from tglow.io.perkin_elmer_parser import PerkinElmerParser
import _tglow_patches  # noqa: F401 — applies upstream bug fixes to PerkinElmerParser


# Plot results from basicpy fit
def plot_basic_results(basic, filename):
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    im = axes[0].imshow(basic.flatfield)
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Flatfield")
    im = axes[1].imshow(basic.darkfield)
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Darkfield")
    axes[2].plot(basic.baseline)
    axes[2].set_xlabel("Frame")
    axes[2].set_ylabel("Baseline")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()

# Plot before and after of an imageset for a given index in the arrays
def plot_before_after(images, images_transformed, i, filename):
    
    images=images.astype(np.float32)
    images_transformed=images_transformed.astype(np.float32)
    images[images==0] = np.nan
    images_transformed[images_transformed==0] = np.nan

    cmap = plt.get_cmap('viridis').copy()

    # Set the color for NaN values (e.g., white)
    cmap.set_bad('darkgrey')

    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow(images[i], cmap=cmap)
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Original")
    im = axes[1].imshow(images_transformed[i], cmap=cmap)
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Corrected")
    fig.suptitle(f"frame {i}")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()

# Plot before and after for two pre-picked images, usefull for max projections
def plot_before_after_mp(images,images_transformed, filename, main="Max projected", zero_na=False):
    
    if zero_na:
        images=images.astype(np.float32)
        images_transformed=images_transformed.astype(np.float32)
        images[images==0] = np.nan
        images_transformed[images_transformed==0] = np.nan
    
    cmap = plt.get_cmap('viridis').copy()

    # Set the color for NaN values (e.g., white)
    cmap.set_bad('darkgrey')
    
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    im = axes[0].imshow(images, cmap=cmap)
    fig.colorbar(im, ax=axes[0])
    axes[0].set_title("Original")
    im = axes[1].imshow(images_transformed, cmap=cmap)
    fig.colorbar(im, ax=axes[1])
    axes[1].set_title("Corrected")
    fig.suptitle(f"{main}")
    fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close()


def plot_flatfield_evaluation(flatfield, image, image_corrected, filename):

    X=flatfield    
    Y=image
    Z=image_corrected
    
    f, ax = plt.subplots(2, 3, figsize=(12, 6))

    im = ax[0,0].imshow(X, cmap='inferno')
    ax[0,0].set_title('Flatfield')
    f.colorbar(im, ax=ax[0,0])

    im = ax[0,1].imshow(Y, cmap='inferno')
    ax[0,1].set_title('Raw image')
    f.colorbar(im, ax=ax[0,1])

    im = ax[0,2].imshow(Z, cmap='inferno')
    ax[0,2].set_title('Corrected image')
    f.colorbar(im, ax=ax[0,2])

    tmp = Z / np.mean(Z)

    im = ax[1,0].imshow(tmp, cmap='inferno')
    #im = ax[1,0].imshow(X/tmp, cmap='inferno')
    #ax[1,0].set_title('Flatfield / (corr / mean(corr))')
    ax[1,0].set_title('corr / mean(corr)')
    f.colorbar(im, ax=ax[1,0])

    # Uncorrected vs flatfield
    ax[1,1].scatter(X.flatten(), Y.flatten(), color="lightgrey")
    m, b = np.polyfit(X.flatten(),  Y.flatten(), 1)
    x_line = np.linspace(min(X.flatten()), max(X.flatten()), 100)
    ax[1,1].plot(x_line, np.repeat(np.mean(Y.flatten()), len(x_line)), color='black', linestyle='--')
    ax[1,1].plot(x_line, m * x_line + b, color='red', linestyle='--')
    ax[1,1].set_xlabel('Bin value of flatfield')
    ax[1,1].set_ylabel('Bin intensity')
    ax[1,1].set_title('Raw image')

    # Corrected vs X
    ax[1,2].scatter(X.flatten(), Z.flatten(), color="lightgrey")
    m, b = np.polyfit(X.flatten(),  Z.flatten(), 1)
    x_line = np.linspace(min(X.flatten()), max(X.flatten()), 100)
    ax[1,2].plot(x_line, np.repeat(np.mean(Z.flatten()), len(x_line)), color='black', linestyle='--')
    ax[1,2].plot(x_line, m * x_line + b, color='red', linestyle='--')
    ax[1,2].set_xlabel('Bin value of flatfield')
    ax[1,2].set_ylabel('Bin intensity')
    ax[1,2].set_title('Corrected image')

    plt.subplots_adjust(left=0.1, right=1.5, top=0.9, bottom=0.1)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class FlatFieldTrainer():
    
    def __init__(self, path, output_dir, output_prefix,  channel,  nimg, merge_n, tune, fit_darkfield, all_planes=False, max_project=False, plates=None, fields=None, planes=None, threshold=False, autosegment=False, plot=False, blacklist=None, pseudoreplicates=None, pseudoreplicates_test=None, blur=False, sigma=2, nimg_validate=0, bins=20):
        
        # Class to read images
        self.provider = CompoundImageProvider(path, nimg,  channel, blacklist, plates, fields, planes, pseudoreplicates, merge_n, max_project, all_planes)
        self.channel = channel
        
        self.nimg=nimg
        self.merge_n=merge_n
        self.nimg_validate=nimg_validate
        self.pseudoreplicates=pseudoreplicates
        self.pseudoreplicates_validate=pseudoreplicates_test
        self.all_planes=all_planes
        self.tune=tune
        self.fit_darkfield=fit_darkfield
        
        self.output_dir=output_dir
        self.output_prefix=output_prefix
    
        self.threshold=threshold
        self.autosegment=autosegment
        self.plot=plot
                
        self.blur=blur
        self.sigma=sigma
        self.bins=bins
        
        # Output
        self.out = f"{self.output_dir}/{self.output_prefix}_ch{str(self.channel)}"
        
            
    def train_basicpy(self):
        
        # Read random images possibly as compound
        training_imgs = self.provider.fetch_training_images()
        
        # Apply a gaussian blur to the images prior to running bp
        if self.blur:
            i=0
            while i < len(training_imgs): 
                training_imgs[i] = gaussian(training_imgs[i], sigma=self.sigma, preserve_range=True)
                i+=1
                
        # Gather the weights
        if self.threshold:
            weights=[]
            k=0
            for img in training_imgs:
                threshold = threshold_otsu(img)
                log.debug(f"Threshold for img {k} = {threshold}")
                weights.append(img>threshold)
                k = k+1
            weights = np.array(weights)
        else:
            weights = None
    
        # Convert into 3d numpy array
        merged = np.array(training_imgs)
        log.info(f"Read training files into array of shape {str(merged.shape)}")

        # Init basicpy object with
        basic = BaSiC(get_darkfield=self.fit_darkfield,
                      autosegment=self.autosegment,
                      working_size=128)

        log.debug("Tuning parameters:")
        log.debug(basic.resize_params)
        log.debug(basic.working_size)

        # Optimze parameters
        if self.tune:
            log.info("Tuning model")

            init_params = {"smoothness_flatfield": 2}
            search_space = {"smoothness_flatfield": list(np.logspace(-0.5, 1.5, 15))}

            if self.fit_darkfield:
                init_params.update({
                    "smoothness_darkfield":  1e-3, 
                    "sparse_cost_darkfield":  1e-3
                })
                
                search_space.update({
                    "smoothness_darkfield":  list(np.logspace(-0.5, 1.5, 15)),
                    "sparse_cost_darkfield":  list(np.logspace(-0.5, 1.5, 15))            
                })
            
            basic.autotune(merged,
                           fitting_weight=weights,
                           #n_iter=200,
                           #histogram_qmin=0.05,
                           #histogram_qmax=0.95,
                           early_stop_torelance=1e-7,
                           search_space=search_space,
                           init_params=init_params,
                           early_stop_n_iter_no_change=25)
        
        # Fit parameters
        log.info("Fitting model")
        basic.fit(merged,
                fitting_weight=weights)
        
        if not os.path.exists(self.out):
            os.makedirs(self.out)
    
        basic.save_model(self.out, overwrite=True)

        if self.plot:
            log.info("Plotting results")
            self.plot_results(basic, merged)
            
        del training_imgs
        del merged 
        gc.collect()
            
        if self.nimg_validate > 0:
            log.info("Running validation on new imageset")
            self.evaluate_flatfield(basic.flatfield, self.bins)
    
    
    def train_mean_filter(self):
        
        # Read random images possibly as compound
        training_imgs = self.provider.fetch_training_images()
        
        sum_image=None
        sum_mask=None

        for final_img in training_imgs:
            final_img = final_img.astype(np.float32) 
            final_img = final_img / np.iinfo(np.uint16).max
            
            if sum_image is None:
                sum_image = np.zeros_like(final_img, dtype=np.float32)
                sum_mask = np.ones_like(final_img, dtype=np.float32)

                sum_image = sum_image + final_img      
                
                if self.threshold: 
                    thresh = threshold_otsu(final_img)
                    sum_mask += (final_img > thresh)

        if self.threshold: 
            avg = sum_image / sum_mask
        else:
            avg = sum_image / len(training_imgs)
        
        # Convert to unit8
        avg_norm = avg / np.max(avg)
        avg_norm = avg_norm * np.iinfo(np.uint8).max
        avg_norm = avg_norm.astype(np.uint8)
        
        blurred = cv.GaussianBlur(avg_norm, (1001,1001),500)
        
        # Normalize to mean to derive scaling factors
        flatfield = blurred / np.mean(flatfield)
        
        self.create_output(flatfield,  np.array(training_imgs))
        
        if self.nimg_validate > 0:
            log.info("Running validation on new imageset")
            self.evaluate_flatfield(flatfield, self.bins)
        
        
    def train_polynomial(self, use_ridge, degree=3, one_model=False, log_trans=False):
    
        training_imgs = self.fetch_images(self.provider)

        downscale=4

        if one_model:
            
            tmp = downscale_local_mean(training_imgs[0], (downscale, downscale))
            Y, X = np.mgrid[:tmp.shape[0], :tmp.shape[1]]
            y_flat = Y.flatten()
            x_flat = X.flatten()
                
            if self.pseudoreplicates >0:
                nimg = self.pseudoreplicates
            else:
                nimg = self.nimg
    
            arr_size=tmp.shape[0] * tmp.shape[1]
                
            log.debug(f"Creating new flat array of size: {arr_size*nimg} of approx {(arr_size*nimg*2)/ float(1024*1024)} mb")
            x = np.empty(tmp.shape[0] * tmp.shape[1] * nimg, dtype=np.uint16)
            y = np.empty(tmp.shape[0] * tmp.shape[1] * nimg, dtype=np.uint16)
            z = np.empty(tmp.shape[0] * tmp.shape[1] * nimg, dtype=np.uint16)
            
            start = 0
            i =0
            for img in training_imgs:
                #log.info(f"Processing {i}")
                if downscale > 1:
                    log.debug(f"Downscaling images by {downscale}x")
                    cur_img = downscale_local_mean(img, (downscale, downscale))
                else:
                    cur_img = img
                #cur_flat = cur_img.flatten()
                
                #keep = cur_flat!=0
                x[start:start+arr_size] = x_flat
                y[start:start+arr_size] = y_flat
                z[start:start+arr_size] = cur_img.flatten()
                start += arr_size
                i +=1
                
            log.info(f"Read images into single flat array of length {len(z)}, last value x/y/i  {x[len(x)-1]}/ {y[len(y)-1]} / {z[len(z)-1]}")
            
            # Scale xy
            x = x/np.max(x)
            y = y/np.max(y)
            
            coeffs = self.fit_poly(x,y,z,
                          trim_quantile=0,
                          use_ridge=use_ridge,
                          remove_zeroes=True,
                          degree=degree,
                          log_trans=log_trans)
            
            # Scale xy
            X = X/np.max(X)
            Y = Y/np.max(Y)
            
            z_fit = P.polyval2d(X, Y, coeffs)
            
            if log_trans:
                z_fit = np.power(2, z_fit, out=z_fit)

            if downscale > 1:
                final_fit = resize(z_fit,
                    training_imgs[0].shape,
                    order=1,
                    mode="edge")
            else:
                final_fit = z_fit
        else:
            fit_sum = None
            i =0
            for img in training_imgs:
                log.info(f"Fitting poly to image {i}")
                
                if downscale > 1:
                    log.debug(f"Downscaling images by {downscale}x")
                    cur_img = downscale_local_mean(img, (downscale, downscale))
                else:
                    cur_img = img
                
                cur_ff = self.fit_poly_img(cur_img, use_ridge=use_ridge, remove_zeroes=self.threshold, degree=degree, trim_quantile=0, log_trans=log_trans)
                #if fit_sum is None:
                #    fit_sum = np.ones_like(cur_img)    
                if fit_sum is None:
                    fit_sum = cur_ff
                else:
                    fit_sum += cur_ff
                i += 1
                
            # Calculate the average over the images
            fit_avg = fit_sum / len(training_imgs)
            
            if downscale > 1:
                # Upsize with bilinear interpolation
                final_fit = resize(fit_avg,
                                training_imgs[0].shape,
                                order=1,
                                mode="edge")
            else:
                final_fit = fit_avg

        # Normalize flatfield to mean
        flatfield = final_fit / np.mean(final_fit)
        
        # For debugging
        #flatfield =  resize_local_mean(flatfield, (50, 50))
        #i=0
        #while i < len(training_imgs):
        #    training_imgs[i] = resize_local_mean(training_imgs[i], (50, 50))
        #    i+=1
        
        log.info("Done, saving output")
        
        # convert to 3d numpy array
        merged = np.array(training_imgs, copy=True)
        del training_imgs
        gc.collect()
        
        self.create_output(flatfield, merged)
        
        del merged
        gc.collect()

        
        if self.nimg_validate > 0:
            log.info("Running validation on new imageset")
            self.evaluate_flatfield(flatfield, self.bins)
        
        log.info("Done, successfully completed")


    def train_pe_index(self, pe_index):
        
        pe_channel = str(int(self.channel)+1)
        parser = PerkinElmerParser(pe_index)
        parser.parse_flatfields(channel=pe_channel)
        
        if (str(pe_channel) not in parser.flatfields.keys()):
            raise ValueError(f"Channel {pe_channel} (ch {self.channel}) not found in PE index among {parser.flatfields.keys()}")
        
        flatfield = parser.flatfields[pe_channel]['flatfield']
        
        # Normalize flatfield to mean
        #flatfield = flatfield / np.mean(flatfield)
        
        log.info("Done, saving output")
        self.plot = False
        self.create_output(flatfield, None)
        
        if self.nimg_validate > 0:
            log.info("Running validation on new imageset")
            self.evaluate_flatfield(flatfield, self.bins)
        
        log.info("Done, successfully completed")


    def fetch_images(self, provider):
        # Read random images possibly as compound
        training_imgs = provider.fetch_training_images()
        
        if self.blur:
            i=0
            pb = tqdm(total=len(training_imgs), desc='Gaussian blur', unit='image')
            while i < len(training_imgs): 
                training_imgs[i] = gaussian(training_imgs[i], sigma=self.sigma, preserve_range=True)
                i+=1
                pb.update(1)
            pb.close()
                
        if self.threshold:
            #size = round(training_imgs[0].shape[0] *0.25)
            #if size % 2 == 0:
            #    size = size-1
            #log.info(f"Local thresholding images with window {size}")
            i=0
            pb = tqdm(total=len(training_imgs), desc='Threshold', unit='image')
            while i < len(training_imgs): 
                
                # Remove outliers/saturated pixels prior to thresholding
                thresh_img=training_imgs[i]
                
                # Remove saturated pixels
                thresh_img[thresh_img == np.iinfo(thresh_img.dtype).max] = 0
                
                # Remove 1% of the most bright pixels
                #q=np.quantile(thresh_img, 0.99)
                #thresh_img[thresh_img > q] = 0
                
                # Localized sauvola threshold. With sparse images, gives issues
                #thresh = threshold_sauvola(training_imgs[i], size)
                
                # Multiotsu threshold, convert to float to allow proper histogramming
                # Otherwise it attempts to fit on a 65k intensity histogram
                thresh_img = thresh_img.astype(np.float32)
                unit16max = np.iinfo(training_imgs[i].dtype).max
                float32max = np.finfo(thresh_img.dtype).max
                
                thresh_img = (thresh_img/unit16max)*float32max
                thresh_raw = threshold_multiotsu(thresh_img, classes=2, nbins=5000)[0]
                thresh = (thresh_raw/float32max) * unit16max

                # Scale down the threshold by 2x to make it more lenient
                thresh = thresh/4
                
                #log.info(f"Otsu tresh/4: {thresh} ({thresh_raw})")
                
                # Regular otsu
                #thresh = threshold_otsu(training_imgs[i])
                
                training_imgs[i][training_imgs[i] < thresh] = 0
                i+=1
                pb.update(1)
            pb.close()
            
            #for img in training_imgs:
            #    log.debug(f"Kept {np.sum(img!=0)} / {img.shape[0]*img.shape[1]} pixels")
         
        return training_imgs

    
    def test(self, input):
        
        log.info("Loading previous flatfield")
        basic = BaSiC.load_model(input)
        
        if self.nimg_validate > 0:
            log.info("Running validation on new imageset")
            self.evaluate_flatfield(basic.flatfield, self.bins)

        log.info("Done, successfully completed")


    def fit_poly_img(self, image, scale_xy=True, trim_quantile=0, use_ridge=False, remove_zeroes=True, degree=2, log_trans=False):
        Y, X = np.mgrid[:image.shape[0], :image.shape[1]]
        Z = image
        
        x = X.flatten()
        y = Y.flatten()
        z = Z.flatten()
        
        if scale_xy:
            x = x / np.max(x)
            y = y / np.max(x)
        
        coeffs = self.fit_poly(x,y,z, trim_quantile=trim_quantile, use_ridge=use_ridge, remove_zeroes=remove_zeroes, degree=degree)
        
        if coeffs is None:
            log.warning("Model did not fit, returning even flatfield")
            return (np.zeros(image.shape))
        
        if scale_xy:
            Y = Y/np.max(Y)
            X = X/np.max(X)

        z_fit = P.polyval2d(X, Y, coeffs)
        
        if log_trans:
            z_fit = np.power(2, z_fit, out=z_fit)
        
        return(z_fit)

    
    def fit_poly(self, x, y, z, trim_quantile=0, use_ridge=False, remove_zeroes=True, degree=2, log_trans=False):
        
        nobs_in = len(z)

        if remove_zeroes:
            keep = z != 0
            x = x[keep]
            y = y[keep]
            z = z[keep]
            log.info(f"Removed {nobs_in - len(z)} zero values")

        if trim_quantile != 0:
            q = np.quantile(z, trim_quantile)
            keep = z < q
            log.info(f"Removed {len(z) - sum(keep)} >q{trim_quantile} (<{q}) values.")
            x = x[keep]
            y = y[keep]
            z = z[keep]
        
        #z = (z - np.mean(z)) / np.std(z)
        ptotal=len(z) / nobs_in
        log.info(f"Fitting moddel with {len(z)} values")
        if (ptotal < 0.05):
            log.warning("Fewer then 5% of pixel values left after filtering")
            #return None
    
        if log_trans:
            log.info(f"Log2 transforming intensity values")
            z = np.log2(z)
    
        if degree <= 0:
            # The cellprofiler model
            # 1 + x2 + y2  + x*y + x + y 
            #design_matrix = np.column_stack((np.ones(len(x)), x*x, y*y, x*y, x, y))
            
            # 4 degree poly PE
            # sourced: https://forum.image.sc/t/how-to-use-illumination-correction-file-from-perkin-elmer-operetta-microscope-outside-of-harmony-in-imagej-or-cellprofiler-instead/31650/8
            # 1 + [x + y] + [xy + x^2 + y^2] + [x^3 + x^2*y + y^2*x + y^3] + [x^4 + x^3*y + x^2*y^2 + y^3*x + y^4]
            design_matrix = np.column_stack((np.ones(len(x)),
                                            x, y, 
                                            x*y, x**2, y**2,
                                            x**3, (x**2)*y, (y**2)*x, y**3,
                                            x**4, (x**3)*y, (x**2)*(y**2), (y**3)*x, y**4))
        else:
            # Alternative full model, has all possible permutations
            # 1 + y + y^2 + x + x*y + x*y^2 + x^2 + x^2*y + x^2*y^2 ...
            design_matrix = P.polyvander2d(x, y, [degree, degree])
        
        if use_ridge:
            model = RidgeCV(alphas=np.logspace(-6, 6, 13), cv=10, fit_intercept=False)
            fit = model.fit(design_matrix, z)
            coef_tmp = fit.coef_
            log.info(f"Best alpha: {fit.alpha_} r2: {fit.best_score_}")   
        else:
            # OLS    
            coef_tmp, residuals, rank, s = np.linalg.lstsq(design_matrix, z, rcond=None)

        # CPR 2d poly
        # 2D matrix, where rows are the power of X and cols are the power of Y
        #coeffs = np.zeros((3, 3))
        #coeffs[0,0] = coef_tmp[0]
        #coeffs[2,0] = coef_tmp[1]
        #coeffs[0,2] = coef_tmp[2]
        #coeffs[1,1] = coef_tmp[3]
        #coeffs[1,0] = coef_tmp[4]
        #coeffs[0,1] = coef_tmp[5]

        if degree <= 0:
            # PE 4d poly
            # 1 + [x + y] + [xy + x^2 + y^2] + [x^3 + x^2*y + y^2*x + y^3] + [x^4 + x^3*y + x^2*y^2 + y^3*x + y^4]
            coeffs = np.zeros((5, 5))
            coeffs[0,0] = coef_tmp[0] # intercept
            coeffs[1,0] = coef_tmp[1] # x
            coeffs[0,1] = coef_tmp[2] # y
            coeffs[1,1] = coef_tmp[3] # xy
            coeffs[2,0] = coef_tmp[4] # x^2
            coeffs[0,2] = coef_tmp[5] # y^2
            coeffs[3,0] = coef_tmp[6] # x^3
            coeffs[2,1] = coef_tmp[7] # x^2*y 
            coeffs[1,2] = coef_tmp[8] # y^2*x 
            coeffs[0,3] = coef_tmp[9] # y^3
            coeffs[4,0] = coef_tmp[10] # x^4  
            coeffs[3,1] = coef_tmp[11] # x^3*y
            coeffs[2,2] = coef_tmp[12] # x^2*y^2
            coeffs[1,3] = coef_tmp[13] # y^3*x
            coeffs[0,4] = coef_tmp[14] # y^4
        else:
            # OR when using full model
            coeffs = coef_tmp.reshape((degree + 1, degree + 1))

        return(coeffs)

    
    def evaluate_flatfield(self, flatfield, bins=20):
        
        provider = copy.copy(self.provider)

        provider.pseudoreplicates = self.pseudoreplicates_validate
        provider.nimg = self.nimg_validate
        #provider.pseudoreplicates = 0
        #provider.merge_n = 1
        #provider.nimg = self.nimg_validate
    
        log.info("Reading images for validation")
        imgs = self.fetch_images(provider)
    
        I=np.array(imgs)
        
        X = flatfield
        #Y = provider.fetch_average_image(max_project=True)
        Y = np.array(imgs)
        Y = Y.astype(np.float32)
        Y[Y==0] = np.nan
        
        if self.threshold:
            plot_before_after(I, Y, 1, filename=f"{self.out}/model_evaluation_thresh_1_nimg_{self.nimg_validate}_nbin_{bins}.png")
            plot_before_after(I, Y, 2, filename=f"{self.out}/model_evaluation_thresh_2_nimg_{self.nimg_validate}_nbin_{bins}.png")
            plot_before_after(I, Y, 3, filename=f"{self.out}/model_evaluation_thresh_3_nimg_{self.nimg_validate}_nbin_{bins}.png")

        Y = np.nanmean(Y, axis=0)
        Z = Y
        Z = Z.astype(np.float32)
        Z = Z / X
        
        Y_plot = Y
        Y_plot[np.isnan(Y_plot)] = 0
        Z_plot = Z
        Z_plot[np.isnan(Z_plot)] = 0
        plot_before_after_mp(Y_plot, Z_plot, filename=f"{self.out}/model_evaluation_pre_post_nimg_{self.nimg_validate}_nbin_{bins}.png")
        
        # Reshape into bins and take the mean
        #X = X.reshape(bins, int(X.shape[0] / bins), bins, int(X.shape[1]/bins)).nanmean(3).nanmean(1)
        #Y = Y.reshape(bins, int(Y.shape[0] / bins), bins, int(Y.shape[1]/bins)).nanmean(3).nanmean(1)
        
        X = X.reshape(bins, int(X.shape[0] / bins), bins, int(X.shape[1]/bins))
        X = np.nanmean(X, axis=3)
        X = np.nanmean(X, axis=1)
        
        Y = Y.reshape(bins, int(Y.shape[0] / bins), bins, int(Y.shape[1]/bins))
        Y = np.nanmean(Y, axis=3)
        Y = np.nanmean(Y, axis=1)
        #X = resize_local_mean(X, (bins, bins))
        #Y = resize_local_mean(Y, (bins, bins))

        # Calculate the correction
        Z = Y
        Z = Z.astype(np.float32)
        Z = Z / X
                
        plot_flatfield_evaluation(X, Y, Z, filename=f"{self.out}/model_evaluation_nimg_{self.nimg_validate}_nbin_{bins}.png")
        

    def create_output(self, flatfield, merged):
        
        log.info(f"Flatfiled min: {np.min(flatfield)} max: {np.max(flatfield)} median: {np.median(flatfield)}")
        
        # Hack, but save the results as a basicpy object to keep it compatible with the pipeline
        basic = BaSiC()
        basic.flatfield = flatfield
        basic.darkfield = np.zeros_like(flatfield)
        
        if merged is None:
            basic.baseline = 1
        else:
            basic.baseline = np.mean(merged, axis=(1,2))
            
        if not os.path.exists(self.out):
            os.makedirs(self.out)
    
        basic.save_model(self.out, overwrite=True)
        plot_basic_results(basic, self.out + "/flat_and_darkfield.png")

        log.info("Output saved")
        
        if self.plot:
            self.plot_results(basic, merged)
            log.info("Plots saved")

        
    def apply_transform(self, basic, images):
        im_float = images.astype(np.float32)
        output = (im_float - basic.darkfield[np.newaxis]) / basic.flatfield[np.newaxis]
        return(output)
        
        
    def plot_results(self, basic, merged):
        
        plot_max = 5
        
        # Apply correction
        merged_corrected = basic.transform(merged)
        #merged_corrected = self.apply_transform(basic, merged)

        # Convert data back to original 16bit uint
        merged_corrected = float_to_16bit_unint(merged_corrected)
        
        # Plot results
        #plot_basic_results(basic, self.out + "/flat_and_darkfield.png")
        plot_before_after_mp(np.max(merged, axis=0), np.max(merged_corrected, axis=0), filename=self.out + "/all_imgs_max_proj_pre_post.png")
        
        # Plot before and after for first 5 images, or as many as are available
        if merged_corrected.shape[0] < plot_max:
            plot_max = merged_corrected.shape[0]
        
        i = 0
        while i < plot_max:
            plot_before_after(merged, merged_corrected, i, filename=self.out + f"/img{i}_pre_post.png")
            i += 1

    

if __name__ == "__main__":
    
    # CLI 
    parser = argparse.ArgumentParser(description="Train a basicpy model on raw HCI images orgnaized into <plate>/<row>/<col>/<field>.ome.tiff stacks with CZYX")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff')
    parser.add_argument('-o','--output', help='Output folder')
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--output_prefix', help='Output prefix name', default="flatfields")
    parser.add_argument('--no_tune', help="[BP only] Do not tune the basicpy model", action='store_true', default=False)
    parser.add_argument('--fit_darkfield', help="[BP only] Fit the darkfield component. For tglow not reccomended", action='store_true', default=False)
    parser.add_argument('--nimg', help="Number of random images to train on. Sampled with replacement", default=None, required=True)
    parser.add_argument('--nimg_test', help="Number of random images to evaluate the final flatfield on. Set to 0 to not evaluate.", default=0, required=False)
    parser.add_argument('--merge_n', help="Number of images to combine into one --nimg times", default=1)
    parser.add_argument('--pseudoreplicates', help="Number of pseudoreplicates to generate. Pseudoreplicate is a compound from --merge_n images samples from --nimg read images", default=0)
    parser.add_argument('--pseudoreplicates_test', help="As --pseudoreplicates, but for the testset", default=0)
    parser.add_argument('-p','--plate', help='Plate to process. Defaults to all detected plates.', nargs='+', default=None)
    parser.add_argument('-c','--channel', help="Channel number to correct, zero indexed", required=True)
    parser.add_argument('--fields', help='Fields to use. Defaults to use all fields.', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. Defaults to use all planes', nargs='+', default=None)
    #parser.add_argument('--gpu', help="Use the GPU", action='store_true', default=False)
    parser.add_argument('--all_planes', help="Instead of randomly picking one plane for a stack, use them all", action='store_true', default=False)
    parser.add_argument('--max_project', help="Calculate flatfields on max projections. Automatically activates --all_planes", action='store_true', default=False)
    parser.add_argument('--threshold', help="For mode MEAN | POLY, use otsu to only calculate mean in foreground regions. For BP Use Otsu threshold to set fitting_weight from basicpy to 1 for foreground and 0 for background. This conceptually the inverse of the autosegmentation", action='store_true', default=False)
    parser.add_argument('--autosegment', help="[BP only] Enable basicpy autosegment option", action='store_true', default=False)
    parser.add_argument('--blur', help="[BP|POLY] Apply a gaussian blur prior to fitting model.", action='store_true', default=False)
    parser.add_argument('--sigma', help="[BP|POLY] The sd of the gaussian kernel.", default=5)
    parser.add_argument('--plot', help="Plot basicpy results", action='store_true', default=False)
    parser.add_argument('--mode', help='The mode to run in BASICPY | MEAN | POLY | PE', default="BASICPY")
    parser.add_argument('--ridge', help="Use ridge regression instead of OLS for POLY", action='store_true', default=False)
    parser.add_argument('--bins', help="The number of 2d bins to use in validating. Reccomend 10-50, use lower number for sparse images", default=20, required=False)
    parser.add_argument('--pe_index', help="PerkinElmer index.xml file to extract flatfiels when --mode PE", default=None, required=False)
    parser.add_argument('--degree', help="[POLY only] The degree for the polynomial. If set to zero, uses PE poly (see docs)", default=0)
    parser.add_argument('--onemodel', help="[POLY only] Combine all images into one flat vector and fit on that", action='store_true', default=False)
    parser.add_argument('--flatfield', help="[TEST] Path to previous results", default=None)

    args = parser.parse_args()
    
    input=args.input
        
    # Set subfolders to process
    if args.plate == None:
        plates=[ name for name in os.listdir(input) if os.path.isdir(os.path.join(input, name)) ]
    else:
        plates=args.plate
    
    # If max projecting all_planes also needs to be true so the full stack is read
    # This is also forced in the constructor, but this is so its printed well here
    # TODO: make the printing a method
    if args.max_project:
        args.all_planes=True
    
    print("-----------------------------------------------------------")
    print("Input plates:\t" + str(plates))
    print("Input fields:\t" + str(args.fields))
    #print("Input planes:\t" + str(args.planes))    
    print("Input:\t\t" + input)
    print("Output:\t\t" + args.output)
    print("Output pre:\t" + args.output_prefix)
    print("channel:\t" + str(args.channel))
    print("n img:\t\t" + str(args.nimg))
    print("n img validate:\t" + str(args.nimg_test))
    print("tune basicpy:\t" + str(not args.no_tune))
    print("merge n:\t" + str(args.merge_n))
    print("darkfied:\t" + str(args.fit_darkfield))
    print("max project:\t" + str(args.max_project))
    print("threshold:\t" + str(args.threshold))
    print("autosegment:\t" + str(args.autosegment))
    print("plot:\t" + str(args.plot))
    print("blacklist:\t" + str(args.blacklist))
    print("pseudoreplicates:\t" + str(args.pseudoreplicates))
    print("pseudoreplicates test:\t" + str(args.pseudoreplicates_test))
    print("gaussian blur:\t" + str(args.blur))
    print("sigma:\t" + str(args.sigma))
    print("bins:\t" + str(args.bins))
    #print("degree:\t" + str(args.degree))

    print("-----------------------------------------------------------")

    trainer = FlatFieldTrainer(path=input,
                            output_dir=args.output,
                            output_prefix=args.output_prefix,
                            channel=args.channel,
                            nimg=int(args.nimg),
                            tune=not args.no_tune,
                            merge_n=int(args.merge_n),
                            fit_darkfield=args.fit_darkfield,
                            all_planes=args.all_planes,
                            max_project=args.max_project,
                            plates=plates,
                            fields=args.fields,
                            threshold=args.threshold,
                            autosegment=args.autosegment,
                            plot=args.plot,
                            blacklist=args.blacklist,
                            pseudoreplicates=int(args.pseudoreplicates),
                            pseudoreplicates_test=int(args.pseudoreplicates_test),
                            blur=args.blur,
                            sigma=float(args.sigma),
                            nimg_validate=int(args.nimg_test),
                            bins=int(args.bins))
    
    # Test mode
    if (args.flatfield is not None):
        trainer.test(args.flatfield)
        exit(0)
    
    # Train modes
    if (args.mode == "BASICPY"):                 
        trainer.train_basicpy()
    elif (args.mode == "MEAN"):
        log.warning("MEAN method is very dodgy and I do not reccomend it, Use it at your own risk!")
        trainer.train_mean_filter()
    elif (args.mode == "POLY"):
        #trainer.train_polynomial(int(args.degree))
        trainer.train_polynomial(use_ridge=args.ridge, degree=int(args.degree), one_model=args.onemodel, log_trans=False)
    elif (args.mode == "PE"):
        
        if args.pe_index is None:
            log.error("When --mode PE, must provide --pe_index to extract flatfields from")
            exit(1)
        
        trainer.train_pe_index(pe_index=args.pe_index)
    else:
        raise ValueError(f"{args.mode}, not a valid mode")






    
    
    
    