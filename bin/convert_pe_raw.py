#!/usr/bin/env python 

import os
import logging
import _tglow_patches  # noqa: F401 — applies upstream bug fixes to PerkinElmerParser
from tglow.io.image_query import ImageQuery
import tglow.io.tglow_io as tglow_io
import argparse
from aicsimageio.types import PhysicalPixelSizes


# Setup logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def main(input_file, output_path, wells):

    log.info("Initializing reader")
    pe_reader = tglow_io.PerkinElmerRawReader(input_file, os.path.dirname(input_file))
    pe_reader.pe_index

    pe_well_index={}
    i = 0
    for well in pe_reader.pe_index.wells:
        pe_well_index[well["id"]] = i
        i+=1
        
    rows = []
    cols = []

    if wells is None:
        for well in pe_reader.pe_index.wells:
            rows.append(well["row"])
            cols.append(well["col"])
    else:
        for well in wells:
            row, col = ImageQuery.well_id_to_index(well)
            rows.append(row)
            cols.append(col)
            
    # Retrieve the channel names to add to the OME metadata    
    # Channel info
    if pe_reader.pe_index.channels is not None:
        channel_names = [f"ch{channel['id']} - {channel['name']}" for channel in pe_reader.pe_index.channels]
        i = 0
        for channel in pe_reader.pe_index.channels:
            i+=1
    else:
        channel_names = None
        
    # ZYX pixel sizes
    resolution=pe_reader.pixel_sizes

    if resolution is not None:
        physical_pixel_sizes=PhysicalPixelSizes(resolution[0], resolution[1], resolution[2])
    else:
        physical_pixel_sizes=None

        
    # AICS writer
    writer = tglow_io.AICSImageWriter(output_path,
                                        channel_names=channel_names,
                                        physical_pixel_sizes=physical_pixel_sizes)
                                            
    
    # Loop over the selected wells
    for row, col in zip(rows, cols):
        
        idx = pe_well_index[f"{str(row).zfill(2)}{str(col).zfill(2)}"]
        
        fields = set(img["field"] for img in pe_reader.pe_index.wells[idx]["images"])
        
        for field in fields:
            q = ImageQuery(pe_reader.pe_index.plate["name"],
                        row,
                        col,
                        field)
            log.info(f"Processing {q.to_string()}")
            
            stack = pe_reader.read_stack(q)
             
            if (stack is not None):           
                writer.write_stack(stack, q)
            else:
                log.warning(f"Stack is None for {q.to_string()}. Skipping writing this stack. This usually happens if an export is incomplete.") 
                
        # Write intensity statistics        
        writer.write_image_stats( ImageQuery(pe_reader.pe_index.plate["name"], row, col, 1))
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Copy raw tiles from Perkin Elmer export to /plate/row/col/field.ome.tiff where field.ome.tiff')
    parser.add_argument('--input_file', type=str, required=False, help='Path to the PE Index file')
    parser.add_argument('--output_path', type=str, required=True, help='Path to the output directory where to copy the files and creaate the folder structure')
    parser.add_argument('-w', '--well', help='Well IDs to merge', default=None, nargs='+')
    
    try:
       args = parser.parse_args()
    except:
       parser.print_help()
       exit(1)

    main(args.input_file, args.output_path, args.well)