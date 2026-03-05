#!/usr/bin/env python 

import os
import re
import tifffile
import string
import argparse
import json
import numpy as np
from basicpy import BaSiC
from pathlib import Path
from skimage import transform, io, exposure, color, util, filters, measure
from matplotlib import pyplot as plt
#from pystackreg import StackReg
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from scipy.ndimage import shift
from skimage.segmentation import clear_border
from skimage.registration import phase_cross_correlation
from matplotlib import pyplot as plt
import numpy as np
from skimage import data, util, transform, feature, measure, filters, metrics

import pickle
import copy
import logging
import pandas as pd
from tglow.io.tglow_io import AICSImageReader
from tglow.io.image_query import ImageQuery
from tglow.utils.tglow_plot import composite_images, plot_registration_imgs
from tglow.utils.tglow_utils import float_to_16bit_unint, apply_registration_cv

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#------------------------------------------------------------------
# Main runner
class Registration:
    
    # Constructor, and parse arguments
    def __init__(self, args):
          
        self.input=args.input
        self.output=args.output
        self.ref_channel_eval = None
        self.qry_channel_eval = None
        self.mode = args.mode
        
        # Plates to process
        self.ref_plate=args.plate
        if len(self.ref_plate) != 1:
            raise RuntimeError("Must only supply one plate for registration as the reference plate")
        self.ref_plate=args.plate[0]
        
        # Channels for alignment and plates to merge
        self.plates_merge=args.plate_merge   

        if self.plates_merge is not None:
            # Reader for plates
            self.plate_reader=AICSImageReader(self.input, plates_filter=[self.ref_plate] + self.plates_merge, fields_filter=args.fields)
        
            if (args.ref_channel is None) | (args.qry_channel is None):
                raise RuntimeError("Need to supply reference and query channels for registration")
            
            self.ref_channel=int(args.ref_channel[0])
            
            if not args.ref_channel_eval is None:
                self.ref_channel_eval = int(args.ref_channel_eval[0])
            
            #self.qry_channel=int(args.qry_channel[0])
            self.qry_channel={}
            self.qry_channel_eval={}
             
            i=0
            for cur_plate in self.plates_merge:
                self.qry_channel[cur_plate] = int(args.qry_channel[i])
                
                if not args.qry_channel_eval is None:
                    self.qry_channel_eval[cur_plate] = int(args.qry_channel_eval[i])
                i+=1      
        else:
            raise RuntimeError("Need to supply --plates_merge for running registration")
            
        
        # Registration options
        self.plot=args.plot
        self.save_reg=True
        #self.transform=StackReg.TRANSLATION
        self.transform='TRANSLATION'
        self.eval_merge=args.eval_merge
        
        if self.eval_merge:
            if (args.ref_channel_eval is None) | (args.qry_channel_eval is None):
                raise RuntimeError("Need to supply reference and query evaluation channels for registration evaluation")      
            
        # Check number of fields matches
        self.fields = self.plate_reader.fields[self.ref_plate]
        
        for plate_merge in self.plates_merge:
            
            if len(self.plate_reader.fields[plate_merge]) != len(self.fields):
                raise RuntimeError("Number of fields in merge plate must match number of fields in reference plate")     
        
        self.offset_x = int(args.offset_x)
        self.offset_y = int(args.offset_y)
        
        
    # Main loops
    def run(self, well):
    
        row, col = ImageQuery.well_id_to_index(well)
        
         # Construct the output dir
        reg_dir=f"{self.output}/{self.ref_plate}/{ImageQuery.ID_TO_ROW[str(row)]}/{col}/"

        alignment_matrices={}
        
        if not os.path.exists(reg_dir):
            os.makedirs(reg_dir)
            log.info(f"Folder created: {reg_dir}")

        # if self.eval_merge:
        #     log.info(f"Calculating pearson correlation between ref and qry for channels {self.ref_channel_eval} and {self.qry_channel_eval}")
        #     outfile = f"{reg_dir}/registration_eval_pearson_correlations.tsv"
        #     file = open(outfile, "w")
        #     file.write("well\trow\tcol\tfield\tref_plate\tqry_plate\tref_ch_eval\tqry_ch_eval\tref_ch\tqry_ch\tr_before\tp_before\tr_after\tt_after\n")
        #     file.close()#
            
        log.info(f"Running alignment between plates on max projections")
        
        if self.eval_merge:
            cols = ['well', 'row', 'col', 'field', 'ref_plate', 'qry_plate', 'ref_ch', 'qry_ch', 'offset_x', "offset_y", 'r_b', 'r_a', 'nobj_ref', 'nobj_qry', 'nobj_qry_a', 'pobj_ba', 'ppix_ba', 'ppix_ab']
            self.eval_data = pd.DataFrame(columns=cols)
            
        for field in self.plate_reader.fields[self.ref_plate][str(row)][str(col)]:
            cur_iq = ImageQuery(self.ref_plate, row, col, field)
            
            alignment_matrices = self.align_field(cur_iq, reg_dir)
            
            if alignment_matrices is None:
                continue
            
            if self.save_reg:
                alignment_matrices['transform'] = self.transform
                pickle.dump(alignment_matrices, open(f"{reg_dir}/{field}.pickle", "wb")) 
            
            if self.eval_merge:
                self.eval_field(cur_iq, alignment_matrices)
        
        if self.eval_merge:
            outfile = f"{reg_dir}/registration_eval.tsv"
            self.eval_data.to_csv(outfile, sep="\t", index=False)
        
    # Perform registration between planes
    def align_field(self, iq, reg_dir):
            
        ref_iq = ImageQuery(self.ref_plate, iq.row, iq.col, iq.field, self.ref_channel)
        
        #log.debug(ref_iq.to_string())
        ref_stack = self.plate_reader.read_image(ref_iq)
    
        # Always max project, as we don't want different alignments between planes
        ref_stack = np.max(ref_stack, axis=0)
        
        alignment_matrices={}
        i=0
        # Loop over plates to merge onto ref
        for plate_merge in self.plates_merge:
            log.info(f"Processing field {iq.field} plate {plate_merge}")
            
            if iq.field not in self.plate_reader.fields[plate_merge][str(iq.row)][str(iq.col)]:
                log.warning(f"Missing field {iq.field} in merge plate {plate_merge}, skipping")
                return None

            if i==0:
                alignment_matrices[plate_merge] = {} 
                i+=1
                
            qry_iq = ImageQuery(plate_merge, iq.row, iq.col, iq.field, self.qry_channel[plate_merge])
            
            qry_stack = self.plate_reader.read_image(qry_iq)
            
            # Always max project, as we don't want different alignments between planes
            qry_stack = np.max(qry_stack, axis=0)
            
            # Pre offset on a known factor in case things are very off.
            if self.offset_x !=0 | self.offset_y !=0:
                log.info(f"Applying offsets x {self.offset_x}, y {self.offset_y}")
                qry_stack_in = shift(qry_stack, (self.offset_y, self.offset_x))
            else:
                qry_stack_in = qry_stack

            # Calculate threshold and mask only pos region in ref to align
            ref_mask = ref_stack > threshold_otsu(ref_stack)
            qry_mask = qry_stack > threshold_otsu(qry_stack)

            #ref_stack = ref_stack * (ref_stack > threshold_otsu(ref_stack))
            #qry_stack_in = qry_stack_in * (qry_stack_in > threshold_otsu(qry_stack_in))
            if self.mode == "STACKREG":
                # Calculate the offsets
                raise NotImplementedError("STACKREG mode is not currently supported, as pystackreg causes a depency conflict in pandas")
                #sr = StackReg(self.transform)
                #align_mat = sr.register(ref_stack, qry_stack_in)
            elif self.mode == "CROSS":
                align_mat = np.eye(3, dtype=np.float64)
                #offset, _, _ = phase_cross_correlation(ref_stack, qry_stack_in, reference_mask=ref_mask, return_error='always')
                offset, _, _ = phase_cross_correlation(ref_stack, qry_stack_in, reference_mask=ref_mask, moving_mask=qry_mask)
                # Offsets are inverted compared to stackreg
                align_mat[0,2] = -offset[1]
                align_mat[1,2] = -offset[0]                                    
            else:
                raise ValueError(f"--mode must be STACKREG | CROSS not {self.mode}")

            # Apply the constant offsets in x and y      
            align_mat[0,2] = align_mat[0, 2] - self.offset_x
            align_mat[1,2] = align_mat[1, 2] - self.offset_y
             
            # Save the offsets
            alignment_matrices[plate_merge] = align_mat

            # Plot 
            if self.plot:
                # Apply the offesets
                #sr = StackReg(self.transform)
                #qry_stack_reg = sr.transform(qry_stack, tmat=align_mat)
                #tform = transform.AffineTransform(matrix=align_mat)
                #qry_stack_reg = transform.warp(qry_stack, tform, order=0, preserve_range=True)
                
                # It is a bit silly to now apply it again with the inverse, but following
                # the stackreg convention, it should be inverted
                # the above should all be identical, but this is what is used downstream, and the opencv version
                # is way quicker then skimage
                qry_stack_reg = apply_registration_cv(qry_stack, np.linalg.inv(align_mat))
                
                before_reg = composite_images([ref_stack, qry_stack])
                after_reg = composite_images([ref_stack, qry_stack_reg])
                
                outfile = f"{reg_dir}/{iq.field}_{plate_merge}_refch{self.ref_channel}_qrych{self.qry_channel[plate_merge]}.png"
                plot_registration_imgs(before_reg, after_reg, filename=outfile)
        
        return alignment_matrices

                
    def eval_field(self, iq, alignment_mats):
        
        i = len(self.eval_data)
        for plate_merge in self.plates_merge:
            
            # Calculate the correlations between the evaluation channels before and after registering
            a = np.max(self.plate_reader.read_image(ImageQuery(self.ref_plate, iq.row, iq.col, iq.field, self.ref_channel_eval)), axis=0)
            c = np.max(self.plate_reader.read_image(ImageQuery(plate_merge, iq.row, iq.col, iq.field, self.qry_channel_eval[plate_merge])), axis=0)  
            
            #sr = StackReg(self.transform)
            #b = sr.transform(c, tmat=alignment_mats[plate_merge])
            b = apply_registration_cv(c, np.linalg.inv(alignment_mats[plate_merge]))

            a = float_to_16bit_unint(a)
            b = float_to_16bit_unint(b)
            c = float_to_16bit_unint(c)

            alab, abin = self.get_objects(a)
            blab, bbin = self.get_objects(b)
            clab, cbin = self.get_objects(c)

            log.info(f"Ref image has {np.max(alab)} objects, query has {np.max(blab)}. {round((np.max(blab)/np.max(alab))*100, 2)}%")
            
            # Intersect between a and b
            ab_ol = (abin * bbin)
            ab_ol_perc = (np.sum(ab_ol) / np.sum(abin)) *100

            ba_ol_perc = (np.sum(ab_ol) / np.sum(bbin)) *100
            
            ac_ol = (abin * cbin)
            ac_ol_perc = (np.sum(ac_ol) / np.sum(abin)) *100
            log.info(f"before/after {round(ac_ol_perc, 2)}/{round(ab_ol_perc, 2)}% of pixels in query mask overlap with ref")

            apcc, apval = measure.pearson_corr_coeff(a*abin, c*abin)
            log.info(f"Before Pearson R: {apcc:0.3g}, p-val: {apval:0.3g}")

            bpcc, bpval = measure.pearson_corr_coeff(a*abin, b*abin)
            log.info(f"After Pearson R {bpcc:0.3g}, p-val: {bpval:0.3g}")
            
            self.eval_data.at[i, "well"] = iq.get_well_id()
            self.eval_data.at[i, "row"] = iq.get_row_letter()
            self.eval_data.at[i, "col"] = iq.col
            self.eval_data.at[i, "field"] = iq.field
            self.eval_data.at[i, "ref_plate"] = self.ref_plate
            self.eval_data.at[i, "qry_plate"] = plate_merge
            self.eval_data.at[i, "field"] = iq.field
            self.eval_data.at[i, "ref_ch"] = self.ref_channel
            self.eval_data.at[i, "qry_ch"] = self.qry_channel[plate_merge]
            self.eval_data.at[i, "offset_x"] = alignment_mats[plate_merge][0, 2]
            self.eval_data.at[i, "offset_y"] = alignment_mats[plate_merge][1, 2]  
            self.eval_data.at[i, "r_b"] = apcc
            self.eval_data.at[i, "r_a"] = bpcc
            self.eval_data.at[i, "nobj_ref"] = np.max(alab)
            self.eval_data.at[i, "nobj_qry"] = np.max(clab)
            self.eval_data.at[i, "nobj_qry_a"] = np.max(blab)
            self.eval_data.at[i, "pobj_ba"] = (np.max(blab)/np.max(alab))*100
            self.eval_data.at[i, "ppix_ba"] = ab_ol_perc
            self.eval_data.at[i, "ppix_ab"] = ba_ol_perc
            
            i+=1
            # line=f"{iq.get_well_id()}\t{iq.row}\t{iq.col}\t{iq.field}\t{self.ref_plate}\t{plate_merge}\t{self.ref_channel_eval}\t{self.qry_channel_eval[plate_merge]}\t{self.ref_channel}\t{self.qry_channel[plate_merge]}\t{apcc}\t{apval}\t{bpcc}\t{bpval}\n"
            # file = open(outfile, 'a')
            # file.write(line)
            # file.close()

    def get_objects(self, img, min_thresh=0.02):
        thresh = threshold_otsu(img)
        
        if np.max(img) < (np.iinfo(img.dtype).max * min_thresh):
            log.warning(f"Max in image smaller then {min_thresh}% of intensity range ({(np.iinfo(img.dtype).max * min_thresh)}, otsu: {thresh}), skipping")
            bw = np.zeros_like(img)
        else:
            bw = closing(img > thresh, square(3))
            bw = clear_border(bw)

        # label image regions
        label_image = label(bw)
        
        return label_image, bw


    # Load a previously generated alignment matrix in {plate: <matrix>} dict
    def load_field(iq, reg_dir):
        
        pickle_path = f"{reg_dir}/{iq.plate}/{ImageQuery.ID_TO_ROW[iq.row]}/{iq.col}/{iq.field}.pickle"
        alignment_matrices = pickle.load(pickle_path)

        return alignment_matrices
    
    def printParams(self):
        
        log.info(f"Input plates:\t{str(self.ref_plate)}")
        log.info(f"Input fields:\t{str(self.fields)}")
        #log.info(f"Input planes:\t{str(self.planes)}")   
        #log.info(f"Input channel:\t{str(self.channels)}")     
        log.info(f"Input:\t\t{self.input}")
        log.info(f"Output:\t\t{self.output}")    

        log.info(f"----------------------------")    
        log.info(f"Merge plates:\t{str(self.plates_merge)}")     
        #log.info(f"Merge channels:\t{str(self.channels_merge)}")  
        log.info(f"Ref channel:\t{str(self.ref_channel)}")  
        log.info(f"Qry channel:\t{str(self.qry_channel)}")
        log.info(f"Ref channel eval:\t{str(self.ref_channel_eval)}")  
        log.info(f"Qry channel eval:\t{str(self.qry_channel_eval)}")  
        log.info(f"Transform:\t{str(self.transform)}")        
        log.info(f"Save registr:\t{str(self.save_reg)}")     
        log.info(f"Plot:\t\t{str(self.plot)}")     


# Main loop
if __name__ == "__main__":
            
    # CLI 
    parser = argparse.ArgumentParser(description="Merge raw Phenix tiff files organized per well into per channel tiffs, one page per z-stack")
    parser.add_argument('-w', '--well', help='Well ID to merge')
    parser.add_argument('-i','--input', help='Base dir to raw input')
    parser.add_argument('-o','--output', help='Output prefix')
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+')
    parser.add_argument('-r','--ref_channel', help='Reference channel for registration (nucleus)', nargs=1, default=None)
    parser.add_argument('-q','--qry_channel', help='Query channel for registration (nucleus)', nargs="+", default=None)
    parser.add_argument('-m','--plate_merge', help='Plates to align using pystackreg (e.g. multiple cycles of the same plate). Output will be added as additional channels.', nargs='+', default=None)
    parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. <plane #1> | [<plane #1> <plane #2> <plane #n>]', nargs='+', default=None)
    parser.add_argument('--plot', help="[OPTIONAL] Plot the before and after registration", action='store_true', default=False)
    parser.add_argument('--eval_merge', help="[OPTIONAL] Calculate correlation between two channels before and after", action='store_true', default=False)
    parser.add_argument('--ref_channel_eval', help='[OPTIONAL] Reference channel for registration evaluation (nucleus)', nargs=1, default=None)
    parser.add_argument('--qry_channel_eval', help='[OPTIONAL] Query channel for registration evaluation (nucleus)', nargs='+', default=None)
    parser.add_argument('--offset_x', help='[OPTIONAL] Constant offset to apply in x (if registration is very off)', default=0)
    parser.add_argument('--offset_y', help='[OPTIONAL] Constant offset to apply in y (if registration is very off)', default=0)
    parser.add_argument('--mode', help='Algorithm for registration, STACKREG | CROSS. Stackreg currently disabled ', default="CROSS")

    args = parser.parse_args()
  
    # args.input = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/MUNCHKIN/pipeline/results/images"
    # args.well = "B08"
    # args.output = "/lustre/scratch125/humgen/projects/cell_activation_tc/projects/MUNCHKIN/pipeline/scripts/misc/testing_2"
    # args.plate = ["250309_230505-V"]
    # args.plate_merge=["250314_015149-V"]
    # args.ref_channel=[2]
    # args.qry_channel=[2]
    # args.ref_channel_eval=[2]
    # args.qry_channel_eval=[2]
    # args.plot=True
    # args.offset_x=0
    # args.offset_y=250
    # args.fields=['1']
    # args.eval_merge=False
    
    log.info("-----------------------------------------------------------")
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        log.info(f"Folder created: {args.output}")
    
    runner = Registration(args)
    log.info(f"Well:\t{args.well}")
    runner.printParams()
    log.info("-----------------------------------------------------------")
    
    runner.run(args.well)

    #main(well, input, output, ref_channel, qry_channel, write_zstack, write_max_projection, basicpy_models, uint32)







