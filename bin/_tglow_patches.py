"""
Upstream bug fixes for tglow-core 0.1.2 applied via monkey-patching.
Import this module early in any script that instantiates PerkinElmerParser.
"""
import logging as _logging
from tglow.io.perkin_elmer_parser import PerkinElmerParser as _P

_log = _logging.getLogger(__name__)


def _parse_wells_patched(self):
    """Skip well image references whose IDs are absent from the acquired images dict.

    Bug: Wells section in Harmony XML lists all planned image slots (e.g. 5 z-planes)
    but only the acquired images (e.g. plane 1 in best-focus mode) appear in the Images
    section. The original code does an unconditional dict lookup, raising KeyError for
    unacquired planes.
    """
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


def _estimate_pixel_sizes_patched(self):
    """Return None gracefully when fewer than 2 z-planes are present.

    Bug: original code unconditionally dereferences img1 (the second-plane image) to
    compute z-step size, but img1 stays None for single-plane (2D / best-focus) data.
    """
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
        zres = zres * 1e6 if zunit == "m" else None
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


_P.parse_wells = _parse_wells_patched
_P.estimate_pixel_sizes = _estimate_pixel_sizes_patched
