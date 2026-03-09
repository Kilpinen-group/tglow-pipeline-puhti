#!/usr/bin/env python 

import numpy as np
import time
import logging
import argparse
import pandas as pd
import glob
import os
from tglow.io.tglow_io import BlacklistReader

# Logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def plot_grouped_density(df, group_col, value_col, filename, vline=None):
    """
    Plots density (KDE) curves for each group in a DataFrame using matplotlib.

    Parameters:
        df (pd.DataFrame): The input DataFrame.
        group_col (str): The column name to group by.
        value_col (str): The column name to plot the density of.
    """
    groups = df[group_col].unique()
    plt.figure(figsize=(10, 6))

    for group in groups:
        data = df[df[group_col] == group][value_col].dropna()
        if len(data) < 2:
            continue  # KDE needs at least 2 data points
        kde = gaussian_kde(data)
        x_min, x_max = data.min(), data.max()
        x_vals = np.linspace(x_min, x_max, 200)
        plt.plot(x_vals, kde(x_vals), label=str(group))

    if vline is not None:
        plt.axvline(x=vline, color='black', linestyle='--', linewidth=2)

    plt.xlabel(value_col)
    plt.ylabel('Density')
    plt.title(f'Density Plot of {value_col} Grouped by {group_col}')
    plt.legend(title=group_col)
    plt.tight_layout()
    plt.savefig(filename)



class ScalingCalculator():
    
    
    def __init__(self, path, pattern, output, path_control=None, pattern_control=None, plate=None, blacklist=None, plate_groups=None, mask_channels=None):

        self.output=output
        self.path=path

        if mask_channels is not None:
            self.mask_channels = {}
            
            for val in mask_channels:
                keypair = val.split("=")
                log.info(f"Adding channel mask {keypair}")
                if keypair[0] not in self.mask_channels:
                    self.mask_channels[keypair[0]]=set()
                self.mask_channels[keypair[0]].add(int(keypair[1]))    
        else:
            self.mask_channels = None    

        # Build list of files
        files = glob.glob(f"{path}/**/{pattern}", recursive=True)
        
        plates = glob.glob(f"{path}/**", recursive=False)
        plates = [os.path.basename(plate) for plate in plates if os.path.isdir(plate)]
        log.info(f"Indexed {len(files)} files over {len(plates)} plates")

        # Read the blacklist
        if blacklist is not None:
            bl_reader = BlacklistReader(blacklist)
            bl = bl_reader.read_blacklist_as_prc()
            log.info(f"Read blacklist with {len(bl)} patterns")
            log.info(f"Blacklist consists of patterns: {bl[:3] if len(bl)>=3 else bl}")

            #files = [file for file in files if not reg]  
            files = [file for file in files if not any(pattern in file for pattern in bl)]
            log.info(f"Filtered using blacklist to {len(files)} files")

        if plate is not None:
            log.info(f"Platelist consists of {len(plate)} patterns {plate[:3] if len(plate)>=3 else plate}")
            files = [file for file in files if any(pattern in file for pattern in args.plate)]
            plates = [plate for plate in plates if any(pattern in plate for pattern in args.plate)]
            log.info(f"Filtered using platelist to {len(files)} files and {len(plates)} plates")
        
        # Set the final filelist
        self.files = files
        self.plates = plates
        
        if path_control is not None:
            self.control_files = glob.glob(f"{path_control}/**/{pattern_control}", recursive=True)
            log.info(f"Indexed {len(self.control_files)} control intensity files")

        # Define the plate groups (one cycle of imaging = a plate group)
        # This ensures channels are scaled witin a cycle and not accross cycles.
        self.plate_groups = []
        if plate_groups is not None:
            df = pd.read_csv(plate_groups, sep="\t")
            
            print(df.head())
                    
            # Extract the first column (ref plate)
            self.plate_groups.append(df["reference_plate"].tolist())

            # Extract the third column (query plates)
            qry_plates = df["query_plates"].tolist()

            # Split the 3rd col on comma, and add each subsequent col
            # as a plategroup
            i = 0
            for qry in qry_plates:
                cur_qry = qry.split(",") 
                
                for j in range(0, len(cur_qry)):
                    if i == 0:
                        self.plate_groups.append([])
                        
                    self.plate_groups[j+1].append(cur_qry[j])
            
                i += 1
        else:
            self.plate_groups.append(plates)
                
        self.main_df=None
        self.control_df=None
        
    def read_intensity_files(self):
        log.info("Reading intensity files")

        # Read per well intensity CSV files
        for file in self.files:
            #log.debug(f"Reading: {file}")
            cur_df = pd.read_csv(file, sep="\t")

            # Override plate column with the directory name to handle cases where
            # the plate ID embedded in the OME-TIFF metadata differs from the
            # directory name used by the pipeline (e.g. microscope acquisition IDs).
            plate_dir = os.path.relpath(file, self.path).split(os.sep)[0]
            if cur_df['plate'].iloc[0] != plate_dir:
                log.warning(f"Plate name mismatch in {file}: '{cur_df['plate'].iloc[0]}' -> '{plate_dir}'")
                cur_df['plate'] = plate_dir

            if self.main_df is not None:
                self.main_df = pd.concat((self.main_df, cur_df))
            else:
                self.main_df = cur_df
    
    def read_control_intensities(self):
        
        if self.control_files is None:
            raise ValueError("Must provide control_files to read control files")
        
        for file in self.control_files:
            cur_df = pd.read_csv(file, sep="\t")
            
            if self.control_df is not None:
                self.control_df = pd.concat((self.control_df, cur_df))
            else:
                self.control_df = cur_df
        
    def calculate_plate_offsets(self, q1, q2, use_col="mean_mean_object_intensity"):
        
        for channel in self.control_df['channel'].unique():
            cur_df = self.control_df[self.control_df['channel'] == channel]
            self.control_df.loc[self.control_df['channel'] == channel, 'plate_offset'] = cur_df[use_col] / cur_df[use_col].min()
        
        tmp = self.control_df[['plate', 'channel', 'plate_offset']]
        
        merged = pd.merge(self.channel_index, tmp, left_on=['ref_plate', 'channel'], right_on=['plate', 'channel'], validate="one_to_one", suffixes=[None, "_y"], how="outer")
        
        #merged['max_scale_total'] = merged['max_scale_total'].fillna(1)
        
        # Set the default to 1 (no scaling)
        # This plate offset descibes how the intensity varies in the controls
        merged['plate_offset'] = merged['plate_offset'].fillna(1)
        
        # Normalize the base plate scale. This descibes the saturation point in each plate after normalization
        merged['max_scale_plate_norm'] =  merged['max_scale_plate'] / merged['plate_offset']
        
        # Derrive the base scale, the point where after equalling plates, the intensity is highest
        # this forms the constant offset applied to all plates
        for channel in merged['channel'].unique():
            cur_df = merged[merged['channel'] == channel]
            merged.loc[merged['channel'] == channel, 'base_scale'] = cur_df['max_scale_plate_norm'].max()
        
        # Set the final scale factor so that the plate with the lowest intensity signal is scaled the most upwards
        # and the plates with the highest intensity as scaled the most downwards, while maintaining the intensity
        # range most optimally based on the base_scale
        merged['scale_factor'] = merged['base_scale'] * merged['plate_offset']
        
        # Set the default to 1 (no scaling)
        merged['scale_factor'] = merged['scale_factor'].fillna(1)
        
        for plate in merged['plate'].unique():
            cur_df = merged[merged['plate'] == plate].copy()
            
            for channel in cur_df['orig_channel'].unique():
                cur_intensity = self.main_df[(self.main_df['channel'] == int(channel)) & (self.main_df['plate'] == plate)]
                
                merged.loc[(merged['plate'] == plate) & (merged['orig_channel'] == channel), 'observed_at_q'] = np.percentile(cur_intensity[q1], q2)

        merged['predicted_at_q_post_scale'] = merged['observed_at_q'] / merged['scale_factor']
        
        self.channel_index = merged
        
        #add to the channel_index
            
    # Build an index of what the new image channels will look like given these settings    
    def build_channel_index(self):
    
        if self.main_df is None:
            self.read_intensity_files()
        # Peak at an image for channel dims and names
        #img = self.plate_reader.get_img(self.plate_reader.images[self.plates[0]][0])      
        
        final_df=None
        
        plate_idx = 0
        for plate in self.plate_groups[0]:
            df = pd.DataFrame(columns=["ref_plate", "plate", "cycle", "channel", "name", "orig_channel", "orig_name"])

            cycle = 1
            channel_id = 0
            tmp_df = self.main_df[self.main_df['plate'].isin([plate])].copy()
            if (len(tmp_df) == 0):
                plate_idx += 1
                continue
            
            # Loop over cycle one channels
            for channel in sorted([int(x) for x in tmp_df['channel'].unique()]):
                df.at[channel, "ref_plate"] = plate
                df.at[channel, "plate"] = plate
                df.at[channel, "channel"] = channel
                df.at[channel, "name"] = ""
                df.at[channel, "cycle"] = cycle
                df.at[channel, "orig_channel"] = channel
                df.at[channel, "orig_name"] = ""

                channel_id += 1  
                            
            # Add the channel names and IDs in the merged plate
            if len(self.plate_groups) != 1:
                
                for plate_group in self.plate_groups[1:]:
                    
                    cycle += 1
                    cur_plate = plate_group[plate_idx]
                    log.debug(cur_plate)
                    tmp_df = self.main_df[self.main_df['plate'].isin([cur_plate])].copy()

                    for channel in sorted([int(x) for x in tmp_df['channel'].unique()]):
                        df.at[channel_id, "ref_plate"] = plate
                        df.at[channel_id, "plate"] = cur_plate
                        df.at[channel_id, "channel"] = channel_id
                        df.at[channel_id, "name"] = ""
                        df.at[channel_id, "cycle"] = cycle
                        df.at[channel_id, "orig_channel"] = channel
                        df.at[channel_id, "orig_name"] = ""
                        channel_id += 1
                        
            # If there are channel to mask
            if self.mask_channels is not None:
                for mask_channel in self.mask_channels[plate]:
                        mask_channel = int(mask_channel)
                        df.at[channel_id, "ref_plate"] = df.iloc[mask_channel]["ref_plate"]
                        df.at[channel_id, "plate"] = df.iloc[mask_channel]["plate"]
                        df.at[channel_id, "channel"] = channel_id
                        df.at[channel_id, "name"] = f"ch{channel_id} - {df.iloc[mask_channel]['name']} mask inclusive"
                        df.at[channel_id, "cycle"] = df.iloc[mask_channel]["cycle"]
                        df.at[channel_id, "orig_channel"] = df.iloc[mask_channel]["channel"]
                        df.at[channel_id, "orig_name"] = df.iloc[mask_channel]["name"]
                        channel_id += 1
                    
                        df.at[channel_id, "ref_plate"] = df.iloc[mask_channel]["ref_plate"]
                        df.at[channel_id, "plate"] = df.iloc[mask_channel]["plate"]
                        df.at[channel_id, "channel"] = channel_id
                        df.at[channel_id, "name"] = f"ch{channel_id} - {df.iloc[mask_channel]['name']} mask exclusive"
                        df.at[channel_id, "cycle"] = df.iloc[mask_channel]["cycle"]
                        df.at[channel_id, "orig_channel"] = df.iloc[mask_channel]["channel"]
                        df.at[channel_id, "orig_name"] = df.iloc[mask_channel]["name"]
                        channel_id += 1

            plate_idx += 1

            if final_df is None:
                final_df=df
            else:
                final_df=pd.concat([final_df, df], axis=0)
        
        
        #df.final_df = final_df["plate"] + "_ch" + str(final_df["orig_channel"])
        self.channel_index=final_df
        
    def save_channel_index(self):
        self.channel_index.to_csv(f"{self.output}/channel_index_with_scaling.tsv", sep='\t', index=False)
              
    def save_scaling_factors(self, scale):
        scaling_factors = [x for x in self.channel_index["ref_plate"] + "_ch" +  self.channel_index["channel"].astype(str) + "=" + self.channel_index[scale].astype(str)]
        # Write scaling factors
        info_file=open(f"{self.output}/scaling_factors.txt", 'w')
        info_file.write(" ".join(scaling_factors))
        info_file.flush()
        info_file.close()
    
    def apply_scaling_factor(self, plot_col, scale):
        
        lookup = {}
        for index, row in self.channel_index.iterrows():
            lookup[f"{row['plate']}_ch{row['orig_channel']}"] = row[scale]
            
        self.main_df['lookup_key'] = self.main_df['plate'].astype(str) + "_ch" + self.main_df['channel'].astype(str)
        self.main_df['scale_value'] = self.main_df['lookup_key'].map(lookup)
        self.main_df[plot_col + "_scaled"] = self.main_df[plot_col] / self.main_df['scale_value']
    
    def save_density_plots(self, plot_col, vline=None):
        
        if not os.path.exists(f"{self.output}/histograms/{plot_col}"):
            os.makedirs(f"{self.output}/histograms/{plot_col}")
            log.info(f"Folder created: {self.output}/histograms/{plot_col}")
                
        for cycle in self.channel_index['cycle'].unique():
            plates = self.channel_index[self.channel_index['cycle'] == cycle]['plate'].unique()
        
            df = self.main_df[self.main_df['plate'].isin(plates)]
            
            for channel in df['channel'].unique():                
                cur_df = df[df['channel'] == channel]
                plot_grouped_density(cur_df, "plate", plot_col, f"{self.output}/histograms/{plot_col}/cycle{cycle}_ch{channel}_{plot_col}_intensity_histograms.png", vline=vline)
            
    def calculate_quantile_scaling_factors(self, q1, q2, scale_max):

        log.info("Considering plate groups:")
        log.info(self.plate_groups)
        
        indices = []
        stats_max=[]
        stats_mean=[]
        stats_median=[]
        scale_factors = {}
        scale_factors_plate = {}

        # Quantiles to report on
        quantiles = [0, 0.01, 0.1, 25, 50, 75, 95, 99, 99.9, 99.99, 99.999, 99.9999, 100]
        
        i = 0
        for plate_group in self.plate_groups:
            
            log.debug(f"current plategroup: {plate_group}")
            
            # Subset to the plates in the group to normalize
            plate_df = self.main_df[self.main_df["plate"].isin(plate_group)].copy()
            log.debug(f"plate_df {plate_df.head}")

            # Find available channels
            for channel in  plate_df["channel"].unique(): 
                
                indices.append(f"group{i}_ch{channel}")
            
                # Subset to the current channel
                cur_df = plate_df[plate_df["channel"] == channel].copy()
            
                # Calculate the scale factor
                scale_factor = np.percentile(cur_df[q1], q2) / scale_max
                
                # Max
                cur_stats = np.percentile(cur_df["q100"], quantiles).tolist()
                cur_stats.append(np.mean(cur_df["q100"]))
                stats_max.append(cur_stats)
                
                # Mean
                cur_stats = np.percentile(cur_df["mean"], quantiles).tolist()
                cur_stats.append(np.mean(cur_df["mean"]))
                stats_mean.append(cur_stats)
                
                # Median
                cur_stats = np.percentile(cur_df["q50"], quantiles).tolist()
                cur_stats.append(np.mean(cur_df["q50"]))
                stats_median.append(cur_stats)
                
                # Scale factors
                for plate in plate_df["plate"].unique():
                    scale_factors[f"{plate}_ch{channel}"]=scale_factor
                    
                    cur_plate_df = plate_df[plate_df['plate'] == plate].copy()
                    cur_plate_df = cur_plate_df[cur_plate_df["channel"] == channel]
                    plate_scale_factor=np.percentile(cur_plate_df[q1], q2) / scale_max
                    scale_factors_plate[f"{plate}_ch{channel}"] = plate_scale_factor 
                    
            i += 1
        
        orig_ch_index=[x for x in self.channel_index["plate"] + "_ch" + self.channel_index["orig_channel"].astype(str)]
        self.channel_index["max_scale_total"] = [scale_factors[key] for key in orig_ch_index]
        self.channel_index["max_scale_plate"] = [scale_factors_plate[key] for key in orig_ch_index]

        raw_scale_factors = [x for x in self.channel_index["ref_plate"] + "_ch" +  self.channel_index["channel"].astype(str) + "=" + self.channel_index["max_scale_total"].astype(str)]
        
        # Write scaling factors
        info_file=open(f"{self.output}/raw_scaling_factors.txt", 'w')
        info_file.write(" ".join(raw_scale_factors))
        info_file.flush()
        info_file.close()

        # Write the summary file        
        header = ['type', 'channel'] + [f"q{str(x)}" for x in quantiles] + ["mean"]
        
        mean_file=open(f"{self.output}/intensity_summary.tsv", 'w')
        mean_file.write("\t".join(header))
        mean_file.write("\n")
        
        for ch in range(len(indices)):
            mean_file.write(f"mean_intensity\t{indices[ch]}\t")
            mean_file.write("\t".join([str(num) for num in stats_mean[ch]]))
            mean_file.write("\n")
            
            mean_file.write(f"median_intensity\t{indices[ch]}\t")
            mean_file.write("\t".join([str(num) for num in stats_median[ch]]))
            mean_file.write("\n")
        
            mean_file.write(f"max_intensity\t{indices[ch]}\t")
            mean_file.write("\t".join([str(num) for num in stats_max[ch]]))
            mean_file.write("\n")
        
        mean_file.flush()
        mean_file.close()
        
        log.info("Done")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Calculate per channel scaling factors based on previously calculated intensity_stats.tsv")
    parser.add_argument('-i','--input', help='Base dir to raw input', required=True)
    parser.add_argument('-o','--output', help='Output prefix', required=True)
    parser.add_argument('-p','--plate', help='Subfolder in raw dir to process', nargs='+', required=True)
    parser.add_argument('--registration_dir', help="Path to registration root storing <plate>/<row>/<col>/<field>.pickle", default=None)
    parser.add_argument('--mask_channels', help='Channels to mask, output will get 2 extra channels per mask_channel, in the order provided. This is applied after registering and flatfield correction, so use the final channel ids', nargs='+', default=None)
    parser.add_argument('--blacklist', help='TSV file with "<plate>  <well>" on each row descrbing what to ignore', default=None)
    
    # Specific arguments
    parser.add_argument('--control_dir', help='Output dir from masked_control_intensity_calculator.py used for plate+channel specific offsets to scale factor organized in <control_dir>/<plate>/<pattern_control>', default=None)
    parser.add_argument('--plate_groups', help='File describing how plates are grouped, used for multicycle runs to ensure scaling is done per cycles. Use registration_manifest.tsv from tglow pipeline', default=None)
    parser.add_argument('--q1', help='Quantile 1, the quantile in an image', default="q99.9999")
    parser.add_argument('--q2', help='Quantile 2, the quantile over all images in input 0-100', default="99")
    parser.add_argument('--pattern', help='File pattern of sumstats tsv', default="intensity_stats.tsv")
    parser.add_argument('--pattern_control', help='File pattern of the control intensties', default="plate_level_control_intensity_summary.tsv")
    parser.add_argument('--scale_max', help='The max possible value of the input images, this will be the max value of the output', default=65535)
    args = parser.parse_args()
    
    args.q2 = float(args.q2)
    args.scale_max = int(args.scale_max)
        
    calculator = ScalingCalculator(path=args.input,
                                    path_control=args.control_dir,
                                    pattern=args.pattern,
                                    pattern_control=args.pattern_control,
                                    output=args.output,
                                    plate=args.plate,
                                    mask_channels=args.mask_channels,
                                    blacklist=args.blacklist,
                                    plate_groups=args.plate_groups)
    
    
    # parser = argparse.ArgumentParser(description="Calculate per channel scaling factors based on previously calculated intensity_stats.tsv")
    # args = parser.parse_args()
    # args.q1="q99.9999"
    # args.q2=99
    # args.scale_max=65535

    # calculator = ScalingCalculator(path="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/decon",
    #                                pattern="intensity_stats.tsv",
    #                                output="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/testing",
    #                                plate=None,
    #                                mask_channels=[0],
    #                                path_control="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/results/scaling/offsets",
    #                                pattern_control="plate_level_control_intensity_summary.tsv",
    #                                blacklist="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/scripts/blacklist.tsv",
    #                                plate_groups="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/DRUG_PERTURB/pipeline/scripts/manifest_registration.tsv")
        
    
    calculator.read_intensity_files()
    calculator.build_channel_index()
    
    if args.control_dir is not None:
        calculator.read_control_intensities()
        
    calculator.calculate_quantile_scaling_factors(args.q1, args.q2, args.scale_max)

    if args.control_dir is not None:
        calculator.calculate_plate_offsets(args.q1, args.q2)
    
    calculator.save_channel_index()
    calculator.save_density_plots(args.q1)
    
    if args.control_dir is not None:
        if calculator.channel_index['plate'].isnull().values.any():
            msg = "NA's found in scaling factors, this should not happen if control intensties are matched, check all arguments are specififed the same between control intensities and this"
            log.error(msg)
            raise RuntimeError(msg)
    
        calculator.save_scaling_factors("scale_factor")
        calculator.apply_scaling_factor(args.q1, "scale_factor")
        
    else:
        if calculator.channel_index['plate'].isnull().values.any():
            msg = "NA's found in scaling factors, usually due to misspecification of arguments, check 'channel_index_with_scaling.tsv' for issues."
            log.error(msg)
            raise RuntimeError(msg)     
        calculator.save_scaling_factors("max_scale_total")
    
        calculator.apply_scaling_factor(args.q1, "max_scale_total")
    
    calculator.save_density_plots(args.q1+"_scaled", vline=np.iinfo(np.uint16).max)

    #calculator.apply_scaling_factor("q99", "scale_factor")
    #calculator.save_density_plots("q99_scaled", vline=np.iinfo(np.uint16).max)



    #calculator.channel_index[["cycle","channel","max_scale_total","plate_offset", "scale_factor"]]