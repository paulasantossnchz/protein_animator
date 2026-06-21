"""
animator.py
Applies cartoon visual settings to a PyMOL morph trajectory, renders each
state as a PNG frame via cmd.png, and assembles the frames into an animated
GIF using Pillow.  Assumes a live PyMOL session started by interpolator.py.
"""

import os
import pymol
from pymol import cmd
from PIL import Image


# ── Constants ─────────────────────────────────────────────────────────────────

_DEFAULT_OUTPUT_DIR = os.path.join("data", "output_animations")
_FRAME_WIDTH        = 800
_FRAME_HEIGHT       = 600
_DEFAULT_FPS        = 25


# ── Auxiliary functions ───────────────────────────────────────────────────────

def _apply_visual_settings(morph_name: str) -> None:
    """
    Configure the PyMOL viewport for a clean cartoon animation.

    Representation : cartoon with fancy helices and smooth loops.
    Coloring        : residue-index spectrum (blue → white → red).
    Background      : black, with antialiasing enabled.
    Camera          : oriented along the principal axes, zoomed to fit.
    """
    cmd.hide("everything", morph_name)
    cmd.show("cartoon",    morph_name)

    cmd.set("cartoon_fancy_helices",       1, morph_name)
    cmd.set("cartoon_smooth_loops",        1, morph_name)
    cmd.set("cartoon_cylindrical_helices", 0, morph_name)

    cmd.spectrum("resi", "blue_white_red", morph_name)

    cmd.bg_color("black")
    cmd.set("antialias",           2)
    cmd.set("ray_opaque_background", 1)
    cmd.set("all_states",          0)

    cmd.orient(morph_name)
    cmd.zoom(morph_name, buffer=5.0)

    return None


def _export_frames(morph_name: str, frames_dir: str) -> list:
    """
    Render each state of *morph_name* as a PNG file inside *frames_dir*.

    Uses cmd.set("state", n) to navigate through states, then cmd.png to
    capture the viewport without ray-tracing for speed.

    Returns the ordered list of paths that were successfully written.
    """
    os.makedirs(frames_dir, exist_ok=True)
    n_states    = cmd.count_states(morph_name)
    frame_paths = []

    for state_idx in range(1, n_states + 1):
        frame_path = os.path.join(frames_dir, f"frame_{state_idx:04d}.png")
        cmd.set("state", state_idx)
        try:
            cmd.png(frame_path, width=_FRAME_WIDTH, height=_FRAME_HEIGHT, ray=0, quiet=1)
            frame_paths.append(frame_path)
        except Exception as exc:
            print(f"[_export_frames] Frame {state_idx:04d} export failed: {exc}")

    print(f"[_export_frames] {len(frame_paths)}/{n_states} frames written to '{frames_dir}'.")
    return frame_paths


def _assemble_gif(frame_paths: list, output_path: str, fps: int) -> bool:
    """
    Combine a sequence of PNG files into a looping animated GIF via Pillow.

    Each image is loaded and converted to RGB (dropping the alpha channel
    produced by PyMOL) before being appended; this detaches it from the file
    handle and, crucially, prevents Pillow from assigning a transparent
    palette index, which otherwise makes the first frame flash with a
    transparent/white background on every loop.

    Returns True on success, False if no frames are available or saving fails.
    """
    assembled = False
    frames    = []

    for fpath in frame_paths:
        if os.path.exists(fpath):
            img = Image.open(fpath)
            img.load()                       # force full read into memory
            frames.append(img.convert("RGB"))  # drop alpha → solid opaque GIF
            img.close()

    if len(frames) > 0:
        frame_duration_ms = max(1, int(1000 / fps))
        try:
            frames[0].save(
                output_path,
                save_all=True,
                append_images=frames[1:],
                optimize=False,
                duration=frame_duration_ms,
                loop=0,
            )
            assembled = True
            print(
                f"[_assemble_gif] {len(frames)}-frame GIF ({fps} fps) "
                f"saved → {output_path}"
            )
        except Exception as exc:
            print(f"[_assemble_gif] GIF save failed: {exc}")
    else:
        print("[_assemble_gif] No frames found; GIF not created.")

    return assembled


def _cleanup_frames(frame_paths: list, frames_dir: str) -> None:
    """Remove temporary PNG files and attempt to delete their directory."""
    for fpath in frame_paths:
        if os.path.exists(fpath):
            os.remove(fpath)

    try:
        os.rmdir(frames_dir)
    except Exception as exc:
        print(f"[_cleanup_frames] Could not remove temp dir '{frames_dir}': {exc}")

    return None


# ── Public API ────────────────────────────────────────────────────────────────

def animate_morph(
    morph_name: str,
    output_dir: str = _DEFAULT_OUTPUT_DIR,
    fps: int        = _DEFAULT_FPS,
) -> str:
    """
    Apply cartoon visuals to a PyMOL morph object and export an animated GIF.

    Expects an active PyMOL session (started by interpolator.interpolate_structures)
    with *morph_name* already loaded.  Does NOT reinitialize PyMOL.

    Workflow
    --------
    1. Verify that *morph_name* exists and has states.
    2. Configure cartoon representation and spectrum coloring.
    3. Render each state to a temporary PNG via cmd.png (no ray-tracing).
    4. Assemble the PNGs into a looping GIF with Pillow.
    5. Delete temporary frames.

    Parameters
    ----------
    morph_name : Name of the PyMOL morph/trajectory object.
    output_dir : Destination directory for the GIF (created if absent).
    fps        : Animation playback speed in frames per second (default 25).

    Returns
    -------
    Absolute (or relative) path of the saved GIF on success, or '' on failure.
    """
    output_path  = ""
    frames_dir   = os.path.join(output_dir, "_tmp_frames")
    gif_path     = os.path.join(output_dir, f"{morph_name}_animation.gif")
    morph_exists = False

    os.makedirs(output_dir, exist_ok=True)

    try:
        n_states     = cmd.count_states(morph_name)
        morph_exists = (n_states > 0)
        print(f"[animate_morph] Found '{morph_name}' with {n_states} states.")
    except Exception as exc:
        print(f"[animate_morph] Object '{morph_name}' not found in PyMOL session: {exc}")

    if morph_exists:
        try:
            _apply_visual_settings(morph_name)
        except Exception as exc:
            print(f"[animate_morph] Visual settings partially failed (continuing): {exc}")

        frame_paths = _export_frames(morph_name, frames_dir)
        assembled   = _assemble_gif(frame_paths, gif_path, fps)
        _cleanup_frames(frame_paths, frames_dir)

        if assembled:
            output_path = gif_path

    return output_path
