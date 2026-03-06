#!/usr/bin/env python 

from cellpose import models, io
import numpy as np
import time
import logging
import argparse
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_utils import float_to_16bit_unint
from skimage.morphology import disk, ball, closing, binary_erosion
from skimage.filters import threshold_otsu, median, rank, threshold_sauvola
from skimage.transform import downscale_local_mean, resize
import os
import math
import tifffile
import skimage.morphology as mp

# Cellpose logger
#logger = io.logger_setup()
logging.getLogger("cellpose").setLevel(logging.INFO)

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def close_per_object(label_image, disk):
    #background = label_image == 0
    #filled = np.zeros_like(label_image)
    for label in np.unique(label_image):
        if label == 0:
            continue
        mask = label_image == label
        mask_filled = closing(mask, disk)
        label_image[mask_filled] = label
    #filled[background] = 0
    return label_image

class CellposeRunner():
    
    def __init__(self, path, plate, output, model, model_nucl, nucl_channel, other_channel, diameter, diameter_nucl,
                do_3d, anisotropy=None, min_cell_area=None, min_nucl_area=None, fields=None, plot=False, cell_flow_thresh=0.4,
                cell_prob_thresh=0, use_nucl_for_declump=True, nucl_power=None, cell_power=None, post_process=True, downsample=None,
                nucl_flow_thresh=0.4, nucl_prob_thresh=0):
        
        self.path=path
        self.plate=plate
        
        self.anisotropy=anisotropy
        self.do_3d=do_3d
        self.fields=fields
        self.__init_reader__()
        
        self.output=output
        self.model=model
        self.model_nucl=model_nucl,
        self.nucl_channel=nucl_channel
        self.other_channel=other_channel
        self.diameter=diameter
        self.diameter_nucl=diameter_nucl
        
        self.cell_flow_thresh=cell_flow_thresh
        self.cell_prob_thresh=cell_prob_thresh
        self.nucl_flow_thresh=nucl_flow_thresh
        self.nucl_prob_thresh=nucl_prob_thresh
        
        self.downsample=downsample

        if self.downsample is not None:
            self.diameter = round(self.diameter / self.downsample)
            
            if self.diameter_nucl is not None:
                self.diameter_nucl = round(self.diameter_nucl / self.downsample)
            
            if self.anisotropy is not None:
                self.anisotropy=self.anisotropy / self.downsample

        # Minimum object area in pixels
        if min_cell_area is None:
            if self.do_3d:
                self.min_size = 4/3 * math.pi * (self.diameter/6)**3
            else:
                self.min_size =  math.pi * ((self.diameter/6))**2
        else:
            self.min_size=int(min_cell_area)
        
        if min_nucl_area is None:
            
            self.min_size_nucl = None
            if self.nucl_channel is not None:
                if self.diameter_nucl is not None:
                    if self.do_3d:
                    # self.min_size_nucl = 4 * math.pi * (self.diameter_nucl/6)**2
                        self.min_size_nucl = (4/3) * math.pi * (self.diameter_nucl/6)**3
                    else:
                        #self.min_size_nucl = 2 * math.pi * ((self.diameter_nucl/6))
                        self.min_size_nucl =  math.pi * (self.diameter_nucl/6)**2
                else:
                    raise Exception("Must set --diameter_nucl when supplying nucleus channel")
        else:
            self.min_size_nucl=int(min_nucl_area)
            
        self.use_nucl_for_declump=use_nucl_for_declump

        log.info(f"Set cell diam:{self.diameter} min area/volume: {self.min_size}")
        log.info(f"Set nucl diam:{self.diameter_nucl} min area/volume: {self.min_size_nucl}")

        self.plot=plot
        
        self.cell_power=cell_power
        self.nucl_power=nucl_power
        self.post_process=post_process
        
    def __init_reader__(self):
        
        self.reader = AICSImageReader(self.path,
                                      plates_filter=self.plate,
                                      fields_filter=self.fields)
        
        if self.anisotropy is None and self.do_3d is True:
            log.info("Estimating anisotropy based on first image")
            
            img1 = self.reader.images[self.plate[0]][0]
            meta = self.reader.get_img(img1)            
            
            if meta.physical_pixel_sizes is not None:
                px = meta.physical_pixel_sizes
                log.info(f"Pixel sizes from ome {px}")
                if px[0] is not None and px[1] is not None and px[1] != 0:
                    self.anisotropy = px[0] / px[1]
                    log.info(f"Estimated anisotropy {self.anisotropy}")
                else:
                    self.anisotropy = 1.0
                    log.warning("Physical pixel sizes incomplete or missing Z/Y values; defaulting anisotropy to 1.0")
            else:
                log.error("Must provide --anisotropy when running in 3d or ome metadata must contain PhysicalPixelSizes")
                raise TypeError("Must provide --anisotropy when running in 3d or ome metadata must contain PhysicalPixelSizes")
        
    def run_model(self, query):
        
        log.debug(f"Processing query {query.to_string()}, row: {query.row},row: {query.get_row_letter()}, col: {query.col}, well: {query.get_well_id()}")
        
        # Read channel with cell outlines
        query.channel = self.other_channel
        log.debug(f"Processing cell channel: {query.channel}")
        
        data_cell = self.reader.read_image(query)
        
        # If max projection is true
        if not self.do_3d:
            data_cell = np.max(data_cell, axis=0)
            data_cell = data_cell[np.newaxis,:,:]
            log.debug(f"cell shape: {data_cell.shape}")
            
        # Read channel with nucleus signal
        if self.nucl_channel:
            query.channel = self.nucl_channel
            log.debug(f"Processing nucl channel: {query.channel}")

            data_nucl = self.reader.read_image(query)
            
            if not self.do_3d:
                data_nucl = np.max(data_nucl, axis=0)
                data_nucl = data_nucl[np.newaxis,:,:]
                log.debug(f"nucl shape: {data_nucl.shape}")
        
        if self.nucl_channel and self.use_nucl_for_declump:
            # Combine into stack
            img = np.stack([data_cell, data_nucl], axis=-1)
            channel_axis=len(img.shape)
            
            # define CHANNELS to run segementation on
            channels = [1,2]
            # channels = [cytoplasm, nucleus]
        else:
            # In the case there is no nucleus channel, make sure there is a channel axis,
            #  so it is always the same shape
            img = data_cell[:,:,:, np.newaxis]
            channels = [0,0]
            channel_axis=None
            
        log.debug(f"Read images in stack of shape {img.shape}")
        log.debug(f"Model: {self.model}")
        log.debug(f"Model nucl: {self.model_nucl[0]}")
        
        img = img.astype(np.float32)

        # Optionally pre-process the data
        cellprob_thresh = self.cell_prob_thresh
        if self.cell_power is not None:
            img = img / np.iinfo(np.uint16).max#np.max(nucl)
            img = img**self.cell_power
            #img = img / np.max(img)
            #img = img / np.percentile(img, 99)
            #img = float_to_16bit_unint(img * np.iinfo(np.uint16).max)
            cellprob_thresh = self.cellprob_thresh-self.cell_power
        
        # Optionally downsample
        image_dims = img.shape
        if self.downsample is not None:
            img = downscale_local_mean(img, (1, self.downsample, self.downsample, 1))
            log.info(f"Downsampled image to {img.shape}, running cellpose with d {self.diameter}")
        
        # Normalize between 1st and 99.99th percentile
        for channel in range(0,img.shape[3]):
            img[:,:,:,channel] = img[:,:,:,channel] - np.percentile(img[:,:,:,channel], 1)
            img[:,:,:,channel] = img[:,:,:,channel] / np.percentile(img[:,:,:,channel], 99.99)
        img[img < 0] = 0

        # Run cellpose
        start_time = time.time()
        masks, flows, styles, diams = self.model.eval(img if self.do_3d else img[0,:,:,:],
                                                diameter=self.diameter,
                                                channels=channels,
                                                channel_axis=channel_axis,
                                                do_3D=self.do_3d,
                                                min_size=self.min_size,
                                                anisotropy=self.anisotropy,
                                                flow_threshold=self.cell_flow_thresh,
                                                cellprob_threshold=cellprob_thresh,
                                                normalize=False)
        log.info("Cell mask running time %s seconds" % round(time.time() - start_time))
        log.info(f"Identfied {np.max(masks)} cells")
        masks = masks if self.do_3d else masks[np.newaxis,:,:]

        # Scale back up
        if self.downsample is not None:
            masks = resize(masks, image_dims[0:3], order=0)
            log.info(f"Upscaled masks back to {masks.shape}")
            
            for plane in range(0, masks.shape[0]):
                masks[plane, :,:] = median(masks[plane, :,:], disk(2))
                
            log.info(f"Smoothed edges with median filter {masks.shape}")

        if not os.path.exists(f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}/"):
            os.makedirs(f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}")

        # save results
        tifffile.imwrite(f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}/{query.field}_cell_mask_d{str(self.diameter)}_ch{self.other_channel}_cp_masks.tiff", masks, compression="zlib")
        
        start_time = time.time()
        
        if self.nucl_channel is not None:
            nucl = data_nucl
            nucl = nucl.astype(np.float32)

            # Binarize cell masks
            masks_bin = (masks > 0).astype(np.uint16)
            
            # Optionally pre process
            cellprob_thresh = self.nucl_prob_thresh
            if self.nucl_power is not None:
                nucl = nucl / np.iinfo(np.uint16).max#np.max(nucl)
                nucl = nucl**self.nucl_power
                # Normalize against the max percentile inside cell objects
                # to avoid normalizing against very overexposed areas.
                tmp = nucl * masks_bin
                nucl = nucl / np.max(tmp)
               # nucl = float_to_16bit_unint(nucl * np.iinfo(np.uint16).max)
                cellprob_thresh = self.nucl_prob_thresh-self.nucl_power

            # Optiionally downsample
            image_dims = nucl.shape
            if self.downsample is not None:
                nucl = downscale_local_mean(nucl, (1, self.downsample, self.downsample))
                log.info(f"Downsampled nucleus to {nucl.shape}, running cellpose with d {self.diameter_nucl}")
            
            # Normalize to 1st and 99.9th percentile
            nucl = nucl - np.percentile(nucl, 1)
            nucl = nucl / np.percentile(nucl, 99.99)
            nucl[nucl < 0] = 0

            nucl_masks, flows, styles, diams = self.model_nucl[0].eval(nucl if self.do_3d else nucl[0,:,:],
                                                    diameter=self.diameter_nucl,
                                                    channels=[0,0],
                                                    do_3D=self.do_3d,
                                                    anisotropy=self.anisotropy,
                                                    min_size=self.min_size_nucl,
                                                    flow_threshold=self.nucl_flow_thresh,
                                                    cellprob_threshold=cellprob_thresh,
                                                    normalize=False)
            log.info("Nucleus running time %s seconds" % round(time.time() - start_time))
            log.info(f"Identfied {np.max(masks)} nuclei")
            
            nucl_masks = nucl_masks if self.do_3d else nucl_masks[np.newaxis,:,:]

            # Close up small gaps when raising to a power
            if self.post_process and self.do_3d:
                # Set a global threshold for nuclei, to remove bits where cellpose fits to the background
                #nucl_thresh = (nucl > threshold_otsu(nucl)).astype(np.uint16)       
                # Replaced this with a local threshold instead
                nucl_thresh =(nucl > threshold_sauvola(nucl, 51)).astype(np.uint16)
                
                disk_size = 10 if self.downsample is None else math.floor(10/self.downsample)

                # Run a 2d closing operation on the threshold mask
                for plane in range(0, nucl_thresh.shape[0]):
                    log.debug(f"Closing {plane}")
                    nucl_thresh[plane,:,:] = closing(nucl_thresh[plane,:,:], disk(disk_size))
                
                # Remove areas of the nucleus mask which don't pass thresh
                nucl_masks = nucl_masks * nucl_thresh
         
                # Remove the area not within the cell
                #nucl_masks = nucl_masks * masks_bin
                
            if self.downsample is not None:
                nucl_masks = resize(nucl_masks, image_dims, order=0)
                log.info(f"Upscaled masks back to {nucl_masks.shape}")
                
                for plane in range(0, nucl_masks.shape[0]):
                    nucl_masks[plane, :,:] = median(nucl_masks[plane, :,:], disk(2))
                
                log.info(f"Smoothed edges with median filter {nucl_masks.shape}")

            tifffile.imwrite(f"{self.output}/{query.plate}/{query.get_row_letter()}/{query.col}/{query.field}_nucl_mask_d{self.diameter_nucl}_ch{self.nucl_channel}_cp_masks.tiff", nucl_masks, compression="zlib")



if __name__ == "__main__":
    
    # CLI 
    parser = argparse.ArgumentParser(description="Run cellpose on images organized into <plate>/<row>/<col>/<field>.ome.tiff stacks with CZYX")
    parser.add_argument('-i','--input', help='Base dir to input organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('-p','--plate', help='Plates to process (at least one)', nargs='+', required=True)
    parser.add_argument('-w','--well', help='Wells to process (at least one)', nargs='+', required=True)
    parser.add_argument('-o','--output', help='Output folder', default="./")
    parser.add_argument('--nucl_channel', help='Channel for nucleus signal', required=False, default=None)
    parser.add_argument('--cell_channel', help='Channel for cell segmentation signal', required=True)
    parser.add_argument('--model', help='Cellpose model', default="cyto2")
    parser.add_argument('--gpu', help="Use the GPU", action='store_true', default=False)
    parser.add_argument('--anisotropy', help="Ratio between z / xy resolution", default=None)
    parser.add_argument('--diameter', help="Estimated cellsize", default=None)
    parser.add_argument('--diameter_nucl', help="Estimated nucleus size", default=None)
    parser.add_argument('--min_cell_area', help="Minimal area or volume of cells. Defaults to area of circle or volume of sphere 1/6th of --diameter", default=None)
    parser.add_argument('--min_nucl_area', help="Minimal area of volume nuclei. Defaults to area of circle or volume of sphere 1/6th of --diameter_nucl", default=None)
    parser.add_argument('--no_3d', help="Don't run in 3d mode", action='store_true', default=False)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    parser.add_argument('--plot', help="Plot overlay masks", action='store_true', default=False)
    parser.add_argument('--cell_flow_threshold', help="Cellpose flow threshold for cells", default=0.4)
    parser.add_argument('--nucl_flow_threshold', help="Cellpose flow threshold for nuclei", default=0.4)
    parser.add_argument('--cell_prob_threshold', help="Cellpose cell probability threshold for cells", default=0)
    parser.add_argument('--nucl_prob_threshold', help="Cellpose cell probability threshold for nuclei", default=0)
    parser.add_argument('--nucl_power', help="Raises the nucleus image to this power as a form of soft thresholding. Sets cellprob threshold to -6.", default=None)
    parser.add_argument('--cell_power', help="Raises the cell image to this power as a form of soft thresholding. Sets cellprob threshold to -6.", default=None)
    parser.add_argument('--dont_use_nucl_for_declump', help="Fit nucleus masks, but do not supply as a 2nd channel to model.", action='store_true', default=False)
    parser.add_argument('--dont_post_process', help="Output the raw celllpose masks, without constraints on nuclei and filtering [3d only]", action='store_true', default=False)
    parser.add_argument('--downsample', help="Downsample the images when running cellpose by this factor, upscale the masks using NN.", default=None)

    args = parser.parse_args()
    
    # Init cellpose model
    model = models.Cellpose(gpu=args.gpu, model_type=args.model)
    #model_nucl = models.Cellpose(gpu=args.gpu, model_type="cyto2")
    model_nucl = model
 
    if not args.no_3d and args.diameter is None:
        e = "Must provide --diameter if running in 3d mode"
        log.error(e)
        raise TypeError(e)
    
    if not args.no_3d and (args.diameter_nucl is None and args.nucl_channel is not None):
        e = "Must provide --diameter_nucl if running in 3d mode and --nucl_channel is specified"
        log.error(e)
        raise TypeError(e)
    
    if args.no_3d:
        log.info("Max projecting image stacks before running cellpose")
        
    # Cellpose runner class
    runner = CellposeRunner(args.input,
                            args.plate,
                            args.output,
                            model,
                            model_nucl,
                            args.nucl_channel,
                            args.cell_channel,
                            int(args.diameter) if args.diameter is not None else args.diameter,
                            int(args.diameter_nucl) if args.diameter_nucl is not None else args.diameter_nucl,
                            not args.no_3d,
                            args.anisotropy,
                            args.min_cell_area,
                            args.min_nucl_area,
                            args.fields,
                            args.plot,
                            float(args.cell_flow_threshold),
                            float(args.cell_prob_threshold),
                            not args.dont_use_nucl_for_declump,
                            float(args.nucl_power) if args.nucl_power is not None else args.nucl_power,
                            float(args.cell_power) if args.cell_power is not None else args.cell_power,
                            not args.dont_post_process,
                            int(args.downsample) if args.downsample is not None else args.downsample,
                            float(args.nucl_flow_threshold),
                            float(args.nucl_prob_threshold),
)
    
    # Loop, ideally one plate and well at the time is supplied, but can run all
    for plate in args.plate:
        for well in args.well:
            row, col = ImageQuery.well_id_to_index(well)    
            #log.debug(f"Detected row_col {row_col} for well {well} in {plate}")
            #log.debug(f"{runner.reader.index.keys()}")

            #fields = runner.reader.index[plate][str(row_col[0])][str(row_col[1])].keys()
            fields = runner.reader.fields[plate][str(row)][str(col)]
            log.info(f"Detected fields {fields} for well {well} in {plate}")
            
            for field in fields:
                q = ImageQuery(plate, row, col, field)
                log.info(f"Running for {q.to_string()}")
                runner.run_model(q)