"""
interpolator.py
Loads PDB structures into a headless PyMOL session, aligns them to a common
spatial frame, reduces them to a shared atom set, and builds a morph
trajectory by simple linear interpolation of atomic coordinates.

The native cmd.morph command is NOT used: it relies on RigiMOL, which is only
available in the paid "Incentive" PyMOL and raises IncentiveOnlyException in
the open-source build declared in requirements.txt.  Linear interpolation of
matched coordinates works in every PyMOL edition and is the "simple
interpolation" suggested in the assignment for bridging dissimilar frames.
"""

import pymol
import numpy as np
from pymol import cmd


# ── Constants ─────────────────────────────────────────────────────────────────

_DEFAULT_MORPH_STEPS = 30
_MORPH_OBJECT        = "morph_trajectory"


# ── Auxiliary functions ───────────────────────────────────────────────────────

def _load_structures(pdb_entries: list) -> list:
    """
    Load each (filepath, chain) entry as a named PyMOL object and clean it so
    only one comparable copy of the protein backbone/side chains remains:
    the requested chain is kept and waters, hydrogens, hetero-atoms (ligands,
    ions, cofactors) and minor alternate conformations are discarded.

    Object names follow the pattern struct_000, struct_001, …
    Returns only the names of objects that were successfully loaded and
    contain at least one atom after filtering.
    """
    object_names = []

    for i, entry in enumerate(pdb_entries):
        filepath, chain = entry
        obj_name        = f"struct_{i:03d}"
        loaded_ok       = True

        try:
            cmd.load(filepath, obj_name)
            cmd.remove(f"({obj_name} and not chain {chain})")
            cmd.remove(f"{obj_name} and (solvent or hydro or not polymer or not alt ''+A)")
            cmd.alter(obj_name, "alt=''")
            n_atoms = cmd.count_atoms(obj_name)
            if n_atoms == 0:
                print(
                    f"[_load_structures] Chain {chain} not found in '{filepath}'; "
                    "skipping."
                )
                loaded_ok = False
            else:
                print(
                    f"[_load_structures] '{filepath}' chain {chain} "
                    f"→ '{obj_name}' ({n_atoms} atoms)"
                )
        except Exception as exc:
            print(f"[_load_structures] Could not load '{filepath}': {exc}")
            loaded_ok = False

        if loaded_ok:
            object_names.append(obj_name)

    return object_names


def _align_structures(object_names: list) -> dict:
    """
    Align every structure to the first one (reference) via cmd.align.

    Returns a dict {object_name: rmsd_float | None} for all objects.
    RMSD for the reference is set to 0.0; None indicates a failed alignment.
    """
    reference  = object_names[0]
    rmsd_table = {reference: 0.0}

    for obj in object_names[1:]:
        try:
            # cmd.align returns (rmsd, n_atoms, n_cycles, n_rejected, n_raw, raw_rmsd)
            result          = cmd.align(obj, reference)
            rmsd_table[obj] = result[0]
            print(
                f"[_align_structures] {obj} → {reference}: "
                f"RMSD = {result[0]:.3f} Å  ({result[1]} atoms paired)"
            )
        except Exception as exc:
            print(f"[_align_structures] Alignment failed for '{obj}': {exc}")
            rmsd_table[obj] = None

    return rmsd_table


def _order_by_rmsd(object_names: list, rmsd_table: dict) -> list:
    """
    Return object names sorted by ascending RMSD from the reference structure
    so that consecutive frames represent the smallest possible structural jump.

    Objects whose alignment failed (RMSD = None) are placed at the end.
    The reference (RMSD = 0.0) always sorts first.
    """
    infinity  = float("inf")
    sort_keys = [
        rmsd_table[name] if (name in rmsd_table and rmsd_table[name] is not None)
        else infinity
        for name in object_names
    ]
    return [name for _, name in sorted(zip(sort_keys, object_names))]


def _atom_keys(obj_name: str) -> set:
    """Return the set of (resi, name) atom identifiers present in *obj_name*."""
    keys = set()
    cmd.iterate(obj_name, "keys.add((resi, name))", space={"keys": keys})
    return keys


def _reduce_to_common_atoms(object_names: list, min_overlap: float = 0.5) -> list:
    """
    Make all structures share an identical, one-to-one atom set so that
    coordinate interpolation is well defined: atom i of one frame and atom i
    of the next are guaranteed to be the same atom of the same residue.

    Done in two stages, both anchored on the reference (object_names[0], the
    structure the user queried):

      1. Drop any structure that shares fewer than *min_overlap* of the
         reference atoms.  This discards fragments, structures with a different
         residue numbering, or a wrongly-picked chain — any one of which would
         otherwise empty the common set and abort the whole animation.
      2. In every surviving structure keep only the (resi, name) atoms present
         in ALL survivors, then sort for a canonical atom order.

    The b-factor column is used as a scratch keep/discard flag (coloring is by
    residue index, so clobbering it is harmless).  Dropped objects are removed
    from the PyMOL session.

    Returns the surviving object names in their original order (reference
    first), or a list shorter than 2 when interpolation is not possible.
    """
    reference = object_names[0]
    ref_keys  = _atom_keys(reference)
    ref_size  = max(len(ref_keys), 1)
    survivors = [reference]

    for obj in object_names[1:]:
        overlap = len(ref_keys & _atom_keys(obj)) / ref_size
        if overlap >= min_overlap:
            survivors.append(obj)
        else:
            print(
                f"[_reduce_to_common_atoms] Dropping '{obj}': only "
                f"{overlap * 100:.0f}% atom overlap with the reference."
            )
            cmd.delete(obj)

    common = set(ref_keys)
    for obj in survivors:
        common = common & _atom_keys(obj)

    n_common = len(common)
    enough   = (len(survivors) >= 2 and n_common > 0)

    if enough:
        for obj in survivors:
            cmd.alter(obj, "b = 1.0 if (resi, name) in common else 0.0",
                      space={"common": common})
            cmd.remove(f"{obj} and b < 0.5")
            cmd.sort(obj)
        print(
            f"[_reduce_to_common_atoms] {len(survivors)} compatible structures, "
            f"{n_common} common atoms each."
        )
        result = survivors
    else:
        print(
            f"[_reduce_to_common_atoms] Not enough compatible structures "
            f"({len(survivors)} survivor(s), {n_common} common atoms)."
        )
        result = []

    return result


def _build_linear_morph(object_names: list, steps: int) -> str:
    """
    Build a multi-state trajectory object by linear interpolation of atomic
    coordinates between consecutive structures.

    Every input structure becomes an exact keyframe; *steps* intermediate
    frames are inserted between each consecutive pair, so two structures that
    differ a lot are bridged by a smooth, simple interpolation instead of an
    abrupt jump.  The interpolation of frame coordinates A and B at fraction
    t is (1 - t)·A + t·B.

    Requires that all objects already share an identical atom set and order
    (see _reduce_to_common_atoms).  Returns the trajectory object name on
    success, or '' if coordinates could not be read.
    """
    morph_name = ""
    template   = object_names[0]
    coords     = [cmd.get_coords(obj, 1) for obj in object_names]
    all_ok     = all(c is not None for c in coords)

    if all_ok:
        cmd.create(_MORPH_OBJECT, template, 1, 1)
        cmd.load_coords(coords[0], _MORPH_OBJECT, 1)
        state = 1

        for k in range(len(coords) - 1):
            start = coords[k]
            end   = coords[k + 1]
            for s in range(1, steps + 1):
                fraction     = float(s) / float(steps)
                interpolated = start * (1.0 - fraction) + end * fraction
                state        = state + 1
                cmd.create(_MORPH_OBJECT, template, 1, state)
                cmd.load_coords(interpolated, _MORPH_OBJECT, state)

        cmd.dss(_MORPH_OBJECT)

        # Remove the per-structure helper objects so only the trajectory
        # remains in the session; otherwise the animator would render the
        # static input structures overlaid on top of the morph.
        for obj in object_names:
            cmd.delete(obj)

        morph_name = _MORPH_OBJECT
        print(
            f"[_build_linear_morph] '{morph_name}' built — "
            f"{cmd.count_states(morph_name)} frames "
            f"({len(object_names)} keyframes, {steps} steps between each)."
        )
    else:
        print("[_build_linear_morph] Could not read coordinates; morph aborted.")

    return morph_name


# ── Public API ────────────────────────────────────────────────────────────────

def interpolate_structures(
    pdb_entries: list,
    morph_steps: int = _DEFAULT_MORPH_STEPS,
) -> str:
    """
    Launch PyMOL in headless mode, load PDB structures (one chain each),
    align them to a common spatial frame, order them by ascending RMSD,
    reduce them to a shared atom set, and build a linearly interpolated morph
    trajectory.

    The PyMOL session remains open so animator.animate_morph can operate on
    the resulting object without re-initializing the application.

    Parameters
    ----------
    pdb_entries : List of (filepath, chain) tuples.  At least 2 required.
    morph_steps : Number of interpolated frames between each consecutive pair
                  of structures (default 30).

    Returns
    -------
    Name of the PyMOL morph object (non-empty string) on success,
    or '' if PyMOL could not be launched, too few structures were loaded, or
    the structures shared no common atoms.
    """
    morph_name   = ""
    pymol_ok     = False
    object_names = []

    try:
        pymol.finish_launching(["pymol", "-cq"])
        pymol_ok = True
    except Exception as exc:
        print(f"[interpolate_structures] PyMOL launch error: {exc}")

    if pymol_ok:
        cmd.reinitialize()
        object_names = _load_structures(pdb_entries)

        if len(object_names) >= 2:
            rmsd_table    = _align_structures(object_names)
            ordered_names = _order_by_rmsd(object_names, rmsd_table)
            survivors     = _reduce_to_common_atoms(ordered_names)

            if len(survivors) >= 2:
                morph_name = _build_linear_morph(survivors, morph_steps)
            else:
                print(
                    "[interpolate_structures] Fewer than 2 compatible structures "
                    "after filtering; cannot interpolate."
                )
        else:
            print(
                f"[interpolate_structures] At least 2 valid structures required; "
                f"only {len(object_names)} loaded successfully."
            )

    return morph_name
