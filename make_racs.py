#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from casatasks import split
from casacore.tables import table, taql
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
import shutil
from tqdm.auto import tqdm

def main(
    emu_dir: Path,
    integration_time: u.Quantity = 15 * u.min,
    field_prefix: str = "EMU",
):
    mss = list(emu_dir.glob("*.ms"))
    for ms in tqdm(mss):
        print(f"Parsing {ms.name}")
        with table((ms / "POINTING").as_posix(), ack=False) as tab:
            # Get ELEVATION for antenna 0
            sels = taql("select ELEVATION, TIME from $tab where ANTENNA_ID == 0 and TRACKING == True")
            elevation = sels.getcol("ELEVATION") * u.deg
            times = Time(sels.getcol("TIME") * u.s, format="mjd")
        
        # Find times around meridian
        peak_time = times[np.argmax(elevation)]
        start = (peak_time - integration_time/2)
        end = (peak_time + integration_time/2)
        outf = Path(ms.as_posix().replace(field_prefix, "RACS"))

        # Split out the data to new MS
        print(f"Splitting {ms.name} to {outf.name}")
        if outf.exists():
            print(f"Removing old {outf.name}")
            shutil.rmtree(outf.as_posix())
        split(
            vis=ms.as_posix(),
            outputvis=outf.as_posix(),
            timerange=f"{start.iso.replace('-', '/').replace(' ', '/')}~{end.iso.replace('-', '/').replace(' ', '/')}",
            datacolumn="all",
        )
        # Rename the FIELD table
        with table((outf / "FIELD").as_posix(), ack=False, readonly=False) as tab:
            names = tab.getcol("NAME")
            new_names = []
            for name in names:
                new_names.append(name.replace(field_prefix, "RACS"))
            tab.putcol("NAME", new_names)

def cli():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("emu_dir", type=Path)
    args = parser.parse_args()
    main(
        emu_dir=Path(args.emu_dir),
    )

if __name__ == "__main__":
    cli()