#!/usr/bin/env python 

import numpy as np
import time
import logging
import argparse
import pandas as pd
import glob
import os
import warnings

from skimage.measure import regionprops_table
from skimage.filters import threshold_otsu
from tglow.io.processed_image_provider import ProcessedImageProvider
from tglow.io.tglow_io import BlacklistReader, AICSImageReader, ControllistReader, ControlRecord
from tqdm import tqdm

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def intensity_mean_nan(mask, intensity):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(intensity[mask])

class MaskedControlIntensityCalculator():
    
    def __init__(self, controls, input, obj_mask_dir, obj_mask_pattern, output, max_project, plates=None, plate_merge=None, flatfields=None, registration_dir=None, mask_channels=None, mask_dir=None,mask_pattern=None, blacklist=None, uint32=False, threshold=False):
        
        self.controls = ControllistReader(controls, plates_filter=plates, blacklist=blacklist).read_controlist()
        self.provider = ProcessedImageProvider(path=input,
                                               blacklist=blacklist,
                                               plate=plates,
                                               plate_merge=plate_merge,
                                               registration_dir=registration_dir,
                                               flatfields=flatfields,
                                               scaling_factors=None,
                                               scaling_slope=None,
                                               scaling_bias=None,
                                               mask_channels=mask_channels,
                                               mask_dir=mask_dir,
                                               mask_pattern=mask_pattern,
                                               uint32=uint32)
                    
        self.cell_mask_reader = AICSImageReader(path=obj_mask_dir, plates_filter=plates, blacklist=blacklist, pattern=obj_mask_pattern)
        self.max_project = max_project
        self.output=output
        self.threshold=threshold

    def make_intensity_summary(self):
        
        df = pd.DataFrame(columns=["plate", "channel", "nobjects", 'nfields', 'nfields_total', 'nwells', 'name', 'mean_threshold', "mean_mean_object_intensity", "mean_integrated_object_intensity", "mean_max_object_intensity"])
        
        i = 0    
        for plate in self.object_intensity['plate'].unique():
            for channel in self.object_intensity['channel'].unique():
                cur_df = self.object_intensity.loc[(self.object_intensity['plate'] == plate) & (self.object_intensity['channel'] == channel)].copy()
                
                # Derrive integrated intensity
                cur_df.loc[:,'intensity_integrated'] = cur_df['area'] * cur_df['intensity_mean']
                cur_df.loc[:,'well_label_field'] = cur_df['label'].astype(str) + "_" + cur_df['well'] + "_" + cur_df['field'].astype(str)
                cur_df.loc[:,'well_field'] = cur_df['well'] + "_" + cur_df['field'].astype(str)

                df.at[i, "plate"] = plate
                df.at[i, "channel"] = channel
                df.at[i, "nobjects"] = len(cur_df['well_label_field'].unique())
                df.at[i, "nfields"] = len(cur_df['field'].unique())
                df.at[i, "nfields_total"] = len(cur_df['well_field'].unique())
                df.at[i, "nwells"] = len(cur_df['well'].unique())
                df.at[i, "name"] = ','.join(map(str, cur_df['name'].unique().tolist()))
                df.at[i, "mean_threshold"] =  cur_df['threshold'].mean()
                df.at[i, "mean_mean_object_intensity"] = cur_df['intensity_mean_nan'].mean()
                df.at[i, "mean_mean_object_intensity_nothresh"] = cur_df['intensity_mean'].mean()
                df.at[i, "mean_integrated_object_intensity"] = cur_df['intensity_integrated'].mean()
                df.at[i, "mean_max_object_intensity"] =  cur_df['intensity_max'].mean()
                i+=1
        
        df.to_csv(f"{self.output}/plate_level_control_intensity_summary.tsv", sep='\t', index=False)

    def fetch_dummy_object_intensity(self):  
        intensity_summary=None
        for control in self.controls:
            for field in self.provider.plate_reader.fields[control.plate][str(control.row)][str(control.col)]:
                for channel in control.channels:
                    cols=["label", "centroid", "area", "intensity_max", "intensity_mean", "intensity_min"]
                    df = pd.DataFrame([[1] * len(cols)], columns=cols)
                    df[cols]=1
                    df['plate'] = control.plate
                    df['well'] = control.well
                    df['field'] = field
                    df['channel'] = channel
                    df['name'] = control.name
                    if intensity_summary is None:
                        intensity_summary = df
                    else:
                        intensity_summary = pd.concat([intensity_summary, df], axis=0)
        self.object_intensity=intensity_summary
        
    def fetch_object_intensity(self):
        intensity_summary=None
        
        pb = tqdm(total=len(self.controls * 9), desc='Reading control images', unit='image')
        for control in self.controls:
            
            for field in self.provider.plate_reader.fields[control.plate][str(control.row)][str(control.col)]:
                img = self.provider.fetch_image(control.get_query(field))
                mask = self.cell_mask_reader.read_image(control.get_query(field))
                pb.update(1)

                # Max project
                if self.max_project:
                    img = np.max(img, axis=1)
                    mask = np.max(mask, axis=1)[0,:,]
                else:
                    mask = mask[0,]
                
                for channel in control.channels:
                    
                    if self.threshold:
                        cur_img = img[channel,].astype(np.float32)
                        thresh = threshold_otsu(cur_img)
                        cur_img[cur_img < thresh] = np.nan
                        
                    else:
                        thresh = 0
                        cur_img = img[channel,]
                    
                    props = regionprops_table(mask,
                                              cur_img,
                                              properties=["label", "centroid", "area", "intensity_max", "intensity_mean", "intensity_min"],
                                              extra_properties=[intensity_mean_nan])
                    
                    df = pd.DataFrame(props)                    
                    log.info(f"Threshold for channel {channel} is {thresh} with {df.shape[0]} objects and mean intensity {df['intensity_mean_nan'].mean()}")
                    #log.info(f"Identified {df.shape[0]} objects")
                    df['plate'] = control.plate
                    df['well'] = control.well
                    df['field'] = field
                    df['channel'] = channel
                    df['name'] = control.name
                    df['threshold'] = thresh

                    if intensity_summary is None:
                        intensity_summary = df
                    else:
                        intensity_summary = pd.concat([intensity_summary, df], axis=0)

        pb.close()
        self.object_intensity=intensity_summary
                
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate per plate & per channel scaling factors based on intensity in foreground region")
    
    # Common options
    parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    parser.add_argument('-o','--output', help='Output prefix', required=True)
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+', required=True)
    parser.add_argument('-m','--plate_merge', help='Plates to combine as multiple cycles of the same plate. Output will be seqeuntially added as additional channels in the order provided', nargs='+', default=None)
    parser.add_argument('--registration_dir', help="Path to registration root storing <plate>/<row>/<col>/<field>.pickle", default=None)
    parser.add_argument('--mask_dir', help="Path to mask root storing masks <plate>/<row>/<col>/<field><mask_pattern>.tiff", default=None)
    parser.add_argument('--mask_pattern', help="The pattern to discover masks, defaults to match run_cellpose.py nucleus. <pattern>.tiff", default="*_nucl_mask_*_cp_masks.tiff")
    parser.add_argument('--mask_channels', help='Channels to mask, output will get 2 extra channels per mask_channel, in the order provided. This is applied after registering and flatfield correction, so use the final channel ids', nargs='+', default=None)
    #parser.add_argument('--fields', help='Fields to use. <field #1> | [<field #1> <field #2> <field #n>]', nargs='+', default=None)
    #parser.add_argument('--planes', help='Z planes to use. <plane #1> | [<plane #1> <plane #2> <plane #n>]', nargs='+', default=None)
    parser.add_argument('--flatfields', help='Basicpy model dir for a channel. If merging channels ids are assigned seqeuntially for extra channels <plate>_ch<channel>=/path/to/model ', nargs='+', default=None)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    parser.add_argument('--uint32', help="Write as 32 bit unsigned integer instead of clipping to 16 bit uint after applying basicpy model", action='store_true', default=False)

    # Specific options
    parser.add_argument('-c','--controls', help='TSV file with <plate> <well> <channel> <name> describing control wells to use per plate', required=True)
    parser.add_argument('--obj_mask_dir', help='Base dir to cell masks organized <plate>/<row>/<col>/<field>.ome.tiff', required=True)
    parser.add_argument('--obj_mask_pattern', help="The pattern to object discover masks, defaults to match run_cellpose.py cell. <pattern>.tiff", default="*_cell_mask_*_cp_masks.tiff")
    parser.add_argument('--max_project', help="Calculate flatfields on max projections. Automatically activates --all_planes", action='store_true', default=False)
    parser.add_argument('--dummy_mode', help="Set all intensities to 1 (equal scaling)", action='store_true', default=False)
    parser.add_argument('--threshold', help="Should background regions be ignored when calculating the object means", action='store_true', default=False)

    args = parser.parse_args()


    for plate in args.plate:
        
        if not os.path.exists(f"{args.output}/{plate}"):
            os.makedirs(f"{args.output}/{plate}")
            log.info(f"Folder created: {args.output}/{plate}")
        
        scaler = MaskedControlIntensityCalculator(controls=args.controls,
                        input=args.input,
                        plates=[plate],
                        plate_merge=args.plate_merge,
                        output=f"{args.output}/{plate}",
                        blacklist=args.blacklist,
                        mask_dir=args.mask_dir,
                        mask_pattern=args.mask_pattern,
                        mask_channels=args.mask_channels,
                        obj_mask_dir=args.obj_mask_dir,
                        obj_mask_pattern=args.obj_mask_pattern,
                        registration_dir=args.registration_dir,
                        flatfields=args.flatfields,
                        uint32=args.uint32,
                        max_project=args.max_project,
                        threshold=args.threshold)

        if args.dummy_mode:
            scaler.fetch_dummy_object_intensity()
        else:
            scaler.fetch_object_intensity()
            
        scaler.make_intensity_summary()
