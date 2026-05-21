#!/usr/bin/env python3
"""
extract_vtk_profiles.py

Automates the manual ParaView workflow:
  Plot Over Line → z-axis at domain centre, resolution 64 → Save Data
  (write timesteps, precision 10, add metadata, add time step, add time)

For every .pvd found in cfvtk/, all timesteps are read, a z-axis profile is
extracted at the x/y midpoint of the domain, and the result is written to a
CSV in the same directory as this script.

Output filename convention (matching the existing CSV files):
  cf_{cpu|gpu}_mfconst_{nu}_dat_{last_idx}_{last_time}.csv
  cf_{cpu|gpu}_mftempdep_dat_{last_idx}_{last_time}.csv

Usage:
    python apps/extract_vtk_profiles.py                          # all simulations
    python apps/extract_vtk_profiles.py tempdep                  # one sim by name fragment
    python apps/extract_vtk_profiles.py "const_0.1"              # one sim by fragment
    python apps/extract_vtk_profiles.py --vtk-dir path/to/cfvtk --out-dir path/
"""

import re
import base64
import struct
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd


# ── VTK binary decoding ───────────────────────────────────────────────────────

_VTK_DTYPE: dict[str, str] = {
    "Float32": "f4",
    "Float64": "f8",
    "Int32":   "i4",
    "UInt32":  "u4",
    "Int64":   "i8",
    "UInt64":  "u8",
}


def _decode_binary_array(text: str, vtk_type: str) -> np.ndarray:
    """Decode a VTK XML binary (base64 + uint32 LE length header) DataArray."""
    raw = base64.b64decode("".join(text.split()))
    n_bytes = struct.unpack_from("<I", raw, 0)[0]
    dtype = np.dtype("<" + _VTK_DTYPE[vtk_type])
    return np.frombuffer(raw, dtype=dtype, offset=4, count=n_bytes // dtype.itemsize)


# ── VTI piece reader ──────────────────────────────────────────────────────────

def _read_vti(vti_path: Path) -> dict:
    """
    Parse a .vti file.

    Returns:
        extent  — (x0, x1, y0, y1, z0, z1) global node extent
        arrays  — {name: ndarray shaped (nz, ny, nx) or (nz, ny, nx, ncomp)}
    """
    tree = ET.parse(vti_path)
    img   = tree.getroot().find("ImageData")
    piece = img.find("Piece")

    ext = list(map(int, img.attrib["WholeExtent"].split()))
    x0, x1, y0, y1, z0, z1 = ext
    nx, ny, nz = x1 - x0, y1 - y0, z1 - z0

    arrays: dict[str, np.ndarray] = {}
    for da in piece.find("CellData").findall("DataArray"):
        name   = da.attrib["Name"]
        ncomp  = int(da.attrib.get("NumberOfComponents", 1))
        vtype  = da.attrib.get("Type") or da.attrib.get("type", "Float64")
        flat   = _decode_binary_array(da.text, vtype)
        # VTK ImageData stores cell data x-varies-fastest → reshape as (nz, ny, nx)
        if ncomp > 1:
            arrays[name] = flat.reshape(nz, ny, nx, ncomp)
        else:
            arrays[name] = flat.reshape(nz, ny, nx)

    return {"extent": (x0, x1, y0, y1, z0, z1), "arrays": arrays}


# ── PVTI / PVD helpers ────────────────────────────────────────────────────────

def _parse_pvd(pvd_path: Path) -> list[tuple[int, int, Path]]:
    """Return [(step_idx, timestep, pvti_path), ...] sorted by timestep."""
    tree = ET.parse(pvd_path)
    result = []
    for idx, ds in enumerate(tree.findall(".//DataSet")):
        result.append((
            idx,
            int(ds.attrib["timestep"]),
            pvd_path.parent / ds.attrib["file"],
        ))
    return result


def _pvti_whole_extent(pvti_path: Path) -> tuple[int, int, int, int, int, int]:
    tree = ET.parse(pvti_path)
    e = list(map(int, tree.find(".//PImageData").attrib["WholeExtent"].split()))
    return tuple(e)


def _pvti_pieces(pvti_path: Path) -> list[tuple[tuple, Path]]:
    """Return [(node_extent, vti_path), ...] for all pieces in a PVTI."""
    tree = ET.parse(pvti_path)
    base = pvti_path.parent
    return [
        (tuple(map(int, p.attrib["Extent"].split())), base / p.attrib["Source"])
        for p in tree.findall(".//Piece")
    ]


# ── Z-profile extraction ──────────────────────────────────────────────────────

def _extract_z_profile(pvti_path: Path) -> dict[str, np.ndarray]:
    """
    Extract a z-axis line profile at the x/y centre of the domain.

    Replicates ParaView's "Plot Over Line" with:
      - line from (nx/2, ny/2, 0) to (nx/2, ny/2, nz)
      - resolution = nz  →  nz+1 sample points at integer z = 0 .. nz

    The probe point at integer z=k falls on the lower face of cell iz=k, which
    VTK resolves to cell iz=k (except z=nz which clamps to the last cell).

    Returns dict of 1-D arrays of length nz+1:
        density, temperature, viscosity, velocity (shape (nz+1, 3)),
        arc_length, points_x, points_y, points_z
    """
    x0g, x1g, y0g, y1g, z0g, z1g = _pvti_whole_extent(pvti_path)
    nx_g = x1g - x0g  # e.g. 128
    ny_g = y1g - y0g  # e.g.  64
    nz_g = z1g - z0g  # e.g.  64

    # Cell just below the x/y midpoint (matching Points:0 = nx_g/2, Points:1 = ny_g/2)
    probe_ix = nx_g // 2 - 1   # e.g. 63
    probe_iy = ny_g // 2 - 1   # e.g. 31

    density     = np.empty(nz_g)
    temperature = np.empty(nz_g)
    viscosity   = np.empty(nz_g)
    velocity    = np.empty((nz_g, 3))

    for (x0, x1, y0, y1, z0, z1), vti_path in _pvti_pieces(pvti_path):
        if not (x0 <= probe_ix < x1 and y0 <= probe_iy < y1):
            continue
        piece  = _read_vti(vti_path)
        lx     = probe_ix - x0
        ly     = probe_iy - y0
        gs, ge = z0 - z0g, z1 - z0g   # global z slice
        arr    = piece["arrays"]
        density    [gs:ge] = arr["density"]    [:, ly, lx]
        temperature[gs:ge] = arr["temperature"][:, ly, lx]
        viscosity  [gs:ge] = arr["viscosity"]  [:, ly, lx]
        velocity   [gs:ge] = arr["velocity"]   [:, ly, lx, :]

    # Pad to nz_g+1 points: last cell is repeated (z=nz_g clamps to cell iz=nz_g-1)
    def _pad1(a: np.ndarray) -> np.ndarray:
        return np.concatenate([a, a[-1:]])

    n = nz_g + 1
    return {
        "density":     _pad1(density),
        "temperature": _pad1(temperature),
        "viscosity":   _pad1(viscosity),
        "velocity":    np.concatenate([velocity, velocity[-1:]], axis=0),
        "arc_length":  np.arange(n, dtype=float),
        "points_x":    np.full(n, float(probe_ix + 1)),  # e.g. 64.0
        "points_y":    np.full(n, float(probe_iy + 1)),  # e.g. 32.0
        "points_z":    np.arange(n, dtype=float),
    }


# ── Output naming ─────────────────────────────────────────────────────────────

def _output_name(sim_name: str, last_idx: int, last_time: int) -> str:
    """Map PVD stem → CSV filename matching the existing naming convention."""
    m_const = re.match(r"couette_flow_(cpu|gpu)_const_([\d.]+)_\d+x\d+x\d+$", sim_name)
    m_tdep  = re.match(r"couette_flow_(cpu|gpu)_tempdep_\d+x\d+x\d+$",         sim_name)
    if m_const:
        plat   = m_const.group(1)
        nu_str = str(float(m_const.group(2)))   # "0.1000" → "0.1"
        return f"cf_{plat}_mfconst_{nu_str}_dat_{last_idx}_{last_time}.csv"
    if m_tdep:
        return f"cf_{m_tdep.group(1)}_mftempdep_dat_{last_idx}_{last_time}.csv"
    return f"cf_{sim_name}_dat_{last_idx}_{last_time}.csv"


# ── Per-simulation processing ─────────────────────────────────────────────────

def _process_pvd(pvd_path: Path, out_dir: Path) -> None:
    sim_name = pvd_path.stem
    print(f"\n{sim_name}", flush=True)

    entries = _parse_pvd(pvd_path)
    if not entries:
        print("  no timesteps — skipping")
        return

    last_idx, last_time = entries[-1][0], entries[-1][1]
    out_path = out_dir / _output_name(sim_name, last_idx, last_time)

    rows: list[dict] = []
    n = len(entries)
    for i, (step_idx, timestep, pvti_path) in enumerate(entries, 1):
        if not pvti_path.exists():
            print(f"  [{i}/{n}] missing {pvti_path.name} — skipping")
            continue

        prof = _extract_z_profile(pvti_path)
        n_pts = len(prof["density"])
        for j in range(n_pts):
            rows.append({
                "TimeStep":          step_idx,
                "Time":              timestep,
                "density":           prof["density"][j],
                "temperature":       prof["temperature"][j],
                "velocity:0":        prof["velocity"][j, 0],
                "velocity:1":        prof["velocity"][j, 1],
                "velocity:2":        prof["velocity"][j, 2],
                "viscosity":         prof["viscosity"][j],
                "vtkValidPointMask": 1,
                "arc_length":        prof["arc_length"][j],
                "Points:0":          prof["points_x"][j],
                "Points:1":          prof["points_y"][j],
                "Points:2":          prof["points_z"][j],
            })

        if i % 25 == 0 or i == n:
            print(f"  [{i}/{n}]", flush=True)

    if not rows:
        print("  no data written")
        return

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False, float_format="%.10g")
    print(f"  → {out_path.name}  ({n} timesteps × {n_pts} pts = {len(df)} rows)")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sim", nargs="?", metavar="SIM",
                    help="name or glob pattern of a single simulation to process "
                         "(e.g. 'couette_flow_cpu_const_0.1*' or 'tempdep'); "
                         "omit to process all simulations")
    ap.add_argument("--vtk-dir", type=Path,
                    default=Path(__file__).parent / "cfvtk",
                    help="directory containing .pvd files  (default: cfvtk/)")
    ap.add_argument("--out-dir", type=Path,
                    default=Path(__file__).parent,
                    help="output directory for CSV files  (default: next to this script)")
    args = ap.parse_args()

    if args.sim:
        pattern = args.sim if args.sim.endswith(".pvd") else f"*{args.sim}*.pvd"
        pvd_files = sorted(args.vtk_dir.glob(pattern))
        if not pvd_files:
            print(f"No .pvd files matching '{pattern}' in {args.vtk_dir}")
            return
    else:
        pvd_files = sorted(args.vtk_dir.glob("*.pvd"))
        if not pvd_files:
            print(f"No .pvd files found in {args.vtk_dir}")
            return

    print(f"Found {len(pvd_files)} simulation(s) in {args.vtk_dir}")
    for pvd in pvd_files:
        _process_pvd(pvd, args.out_dir)
    print("\nDone.")


if __name__ == "__main__":
    main()
