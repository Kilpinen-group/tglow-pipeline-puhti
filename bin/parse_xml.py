#!/usr/bin/env python 

# Adapted from Martin Prete's scripts
# Works with both Opera Phenix and Operetta indices
####################################################################
# Script invocation:
# python actions/parse_xml.py \
#    --input_file "/path/to/Images/Index.xml" \
#    --output_path "/path/to/output"
####################################################################
# Store Index.xml information in JSON format for easier parsing
# and re-use of the pate/well/images metadata.
####################################################################

import os
import re
import json
import logging
import argparse
from xml.etree import ElementTree as ET
from tglow.io.perkin_elmer_parser import PerkinElmerParser
from datetime import datetime

# Upstream bug fix (tglow-core 0.1.2): parse_wells() crashes with KeyError when the
# Wells section of the XML references image IDs that are absent from the Images section
# (e.g. unacquired z-planes in best-focus acquisition mode).
# Fix: skip image references whose ID is not in self.images.
def _parse_wells_patched(self):
    import logging as _logging
    _log = _logging.getLogger(__name__)
    _log.info("[+] Reading Wells metadata")
    self.wells = []
    for well in self.xml.findall("./PE:Wells/PE:Well", self.NS):
        w = {
            "id":  well.find("./PE:id",  self.NS).text,
            "row": int(well.find("./PE:Row", self.NS).text),
            "col": int(well.find("./PE:Col", self.NS).text),
            "images": [
                self.images[wi.attrib["id"]]
                for wi in well.findall("./PE:Image", self.NS)
                if wi.attrib["id"] in self.images
            ],
        }
        self.wells.append(w)
    _log.info(f" └ Wells: {len(self.wells)}")

PerkinElmerParser.parse_wells = _parse_wells_patched


# Upstream bug fix (tglow-core 0.1.2): estimate_pixel_sizes() crashes when the
# acquisition has only one z-plane (img1 stays None, then dereferenced unconditionally).
# Fix: return None early when no second plane is present (2D / best-focus data).
def _estimate_pixel_sizes_patched(self):
    import logging as _logging
    _log = _logging.getLogger(__name__)
    img0 = None
    img1 = None
    for img in self.wells[0]["images"]:
        if img["plane"] == '1':
            img0 = img
        if img["plane"] == '2':
            img1 = img
            break
    if img0 is None or img1 is None:
        _log.warning("Could not estimate z resolution: fewer than 2 planes found")
        zres = None
    else:
        zres = abs(img0["position"]["z"]["value"] - img1["position"]["z"]["value"])
        zunit = img1["position"]["z"]["unit"]
        if zunit == "m":
            zres = zres * 1e6
        else:
            _log.warning(f"Could not estimate z resolution, invalid unit {zunit}")
            zres = None
    if self.channels is not None and len(self.channels) > 0:
        yres = self.channels[0]["image_resolution"]["y"]["value"]
        yunit = self.channels[0]["image_resolution"]["y"]["unit"]
        yres = yres * 1e6 if yunit == "m" else None
        xres = self.channels[0]["image_resolution"]["x"]["value"]
        xunit = self.channels[0]["image_resolution"]["x"]["unit"]
        xres = xres * 1e6 if xunit == "m" else None
    else:
        yres = None
        xres = None
    if zres is not None and yres is not None and xres is not None:
        return [zres, yres, xres]
    return None

PerkinElmerParser.estimate_pixel_sizes = _estimate_pixel_sizes_patched
import dateutil
#from tglow.io.image_query import ImageQuery

# Setup logging
logging.basicConfig(format='%(asctime)s %(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse Perkin Elmer index XML file to JSON file')
    parser.add_argument('--input_file', type=str, required=False, help='Path to the PE Index file')
    parser.add_argument('--output_path', type=str, required=True, help='Path to the output directory where to storer the JSON file')
    parser.add_argument('--to_manifest', required=False, action='store_true', default=False, help='Store the output as a <well> <plate> <index> tsv file')
    
    # for debugging
   #args.input_file="/lustre/scratch125/humgen/projects/cell_activation_tc/projects/KITTY/pipeline/results/images/jm52_KITTY_20230929_4days/Index.xml"
    #args.output_path="/software/teamtrynka/installs/tglow-pipeline/bin/output/testing"
    #args.to_manifest=True
    try:
       args = parser.parse_args()
    except:
       parser.print_help()
       exit(1)

    pe_data = PerkinElmerParser(args.input_file)
    
    if args.to_manifest:
        pe_data.write_manifest(args.output_path)
        pe_data.save(args.output_path)
        
        # Channel info
        info_file=open(f"{args.output_path}/acquisition_info.txt", 'w')
        
        info_file.write("#-----------------------------------------------------------------------------------\n")
        info_file.write("# Channel list:\n")
        if pe_data.channels is not None:
            
            i = 0
            info_file.write(f"index\tpe_id\tname\tnew_name\tresolution_zyx\texcitation\temission\texposure\n")
            for channel in pe_data.channels:
            
                res=channel['size']
                res = (len(pe_data.planes), res[1], res[0])
            
                info_file.write(f"{i}\t{channel['id']}\t{channel['name']}\tch{channel['id']} - {channel['name']}\t{res}\t{channel['main_excitation_wavelength']}\t{channel['main_emission_wavelength']}\t{channel['exposure_time']}\n")
                i+=1
        else:
            channel_names = None
            info_file.write("None detected\n")
            
        info_file.write("\n")
        
        # Objective info
        info_file.write("#-----------------------------------------------------------------------------------\n")
        info_file.write("# Objective info:\n")
        if pe_data.channels is not None:
            channel = pe_data.channels[0]
            info_file.write(f"{'acquisition_type:':<35} {channel['acquisition_type']}\n")
            info_file.write(f"{'binning_x:':<35} {channel['binning_x']}\n")
            info_file.write(f"{'binning_y:':<35} {channel['binning_y']}\n")
            info_file.write(f"{'objective_na:':<35} {channel['objective_na']}\n")
            info_file.write(f"{'objective_magnification:':<35} {channel['objective_magnification']}\n")
        else:
            channel_names = None
            info_file.write("None detected\n")
        
        info_file.write("\n")

        # ZYX pixel sizes
        resolution = pe_data.estimate_pixel_sizes()
        
        info_file.write("#-----------------------------------------------------------------------------------\n")
        info_file.write("# Pixel sizes: \n")
        if resolution is not None:
            info_file.write(f"{'ZYX:':<35}{resolution}\n")
            info_file.write(f"{'Estimated anisotropy:':<35}{resolution[0] / resolution[1]}\n")
        else:
            physical_pixel_sizes=None
            info_file.write("None detected\n")
                    
        info_file.write("\n")

        # Plate info
        
        info_file.write("#-----------------------------------------------------------------------------------\n")
        info_file.write("# Plate info: \n")
        
        if pe_data.wells is not None:
            datetimes = []
            channels = set()
            fields = set()
            planes = set()
            
            for well in pe_data.wells:
                for image in well["images"]:
                    date_string = image["time_abs"]
                    #datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S %z")
                    #datetimes.append(datetime.fromisoformat(date_string))
                    datetimes.append(dateutil.parser.parse(date_string))
                    channels.add(image['channel'])
                    fields.add(image['field'])
                    planes.add(image['plane'])

            start = min(datetimes)
            end = max(datetimes)
            delta = end - start
            per_well = delta / len(pe_data.wells)
            per_image = delta / len(datetimes)
            per_field = per_image * len(planes) * len(channels)
            
            info_file.write(f"{'Plate name:':<35}{pe_data.plate['name']}\n")
            info_file.write(f"{'Plate type:':<35}{pe_data.plate['type']}\n")
            info_file.write(f"{'Number of images:':<35}{len(datetimes)}\n")
            info_file.write(f"{'Number of wells:':<35}{len(pe_data.wells)}\n")
            info_file.write(f"{'Number of images per well:':<35}{len(fields) * len(channels) * len(planes)}\n")
            info_file.write(f"{'Number of fields per well:':<35}{len(fields)}\n")
            info_file.write(f"{'Number of channels per image:':<35}{len(channels)}\n")
            info_file.write(f"{'Number of planes per image:':<35}{len(planes)}\n")

            info_file.write("\n")
            info_file.write("#-----------------------------------------------------------------------------------\n")
            info_file.write("# Plate timings: \n")

            info_file.write(f"{'Plate start:':<25}{ start.strftime('%Y-%m-%d %H:%M:%S')}\n")
            info_file.write(f"{'Plate end:':<25}{ end.strftime('%Y-%m-%d %H:%M:%S')}\n")
            info_file.write(f"{'Total time (H:M:S)':<25}{str(delta)}\n")
            info_file.write(f"{'Time per well (H:M:S)':<25}{str(per_well)}\n")
            info_file.write(f"{'Time per field (H:M:S)':<25}{str(per_field)}\n")


        info_file.flush()
        info_file.close() 
        
    else:
        pe_data.save(args.output_path)

