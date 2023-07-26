#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert ASKAP long track data to 'RACS' data
"""
import shutil
from pathlib import Path
from typing import Tuple

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from casacore.tables import table, taql
from casatasks import split
from tqdm.auto import tqdm


def find_meridian_time(
    ms: Path, integration_time: u.Quantity = 15 * u.min
) -> Tuple[Time, Time]:
    """Find the time around the meridian for a given MS and integration time

    Args:
        ms (Path): MS to find meridian time for
        integration_time (u.Quantity, optional): Integration time. Defaults to 15*u.min.

    Returns:
        Tuple[Time, Time]: Start and end times around meridian
    """
    with table((ms / "POINTING").as_posix(), ack=False) as tab:
        # Get ELEVATION for antenna 0
        sels = taql(
            "select ELEVATION, TIME from $tab where ANTENNA_ID == 0 and TRACKING == True"
        )
        elevation = sels.getcol("ELEVATION") * u.deg
        times = Time(sels.getcol("TIME") * u.s, format="mjd")

    # Find times around meridian
    peak_time = times[np.argmax(elevation)]
    start = peak_time - integration_time / 2
    end = peak_time + integration_time / 2
    return start, end


def long_to_racs(
    ms: Path,
    start: Time,
    end: Time,
    field_prefix: str = "EMU",
):
    """Convert a long MS to a RACS MS

    Args:
        ms (Path): Long track MS
        start (Time): Start time
        end (Time): End time
        field_prefix (str, optional): Field name prefix. Defaults to "EMU".
    """
    outf = Path(ms.as_posix().replace(field_prefix, "RACS"))
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


def make_racs(
    data_dir: Path,
    integration_time: u.Quantity = 15 * u.min,
    field_prefix: str = "EMU",
):
    mss = list(data_dir.glob("*.ms"))
    for ms in tqdm(mss):
        print(f"Parsing {ms.name}")
        # Find the times around the meridian
        start, end = find_meridian_time(ms, integration_time=integration_time)
        # Split out the data to new MS
        long_to_racs(
            ms=ms,
            start=start,
            end=end,
            field_prefix=field_prefix,
        )


def cli():
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("data_dir", type=Path)
    parser.add_argument(
        "--integration_time", type=float, default=15, help="Integration time in minutes"
    )
    parser.add_argument(
        "--field_prefix", type=str, default="EMU", help="Field name prefix"
    )
    args = parser.parse_args()
    make_racs(
        data_dir=Path(args.data_dir),
        integration_time=args.integration_time * u.min,
        field_prefix=args.field_prefix,
    )


if __name__ == "__main__":
    cli()
