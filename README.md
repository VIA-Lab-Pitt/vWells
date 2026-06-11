<p align="center">
  <img src="assets/treetrack-logo-light.svg" alt="Treetrack" width="500">
</p>

<p align="center">
  A pipeline for segmenting and visualizing brain vasculature from MRI scans.
</p>

Treetrack is a set of software for semi-automated binary segmentation of 3D medical images using "variance wells" (vWells). We combine statistical shortest path and region-growing algorithms with 3D visualization overlaid on the source image to enable efficient user-guided segmentation, validation, and analysis. There is a simplified version in 2d available at the following page: https://github.com/SatyajBhargava/2024-2D-vWell-Algorithm- 

There are two ways to use Treetrack:

- vWell Viewer (installer located in releases) — a standalone app for inspecting a scan and its saved vessel segmentation in 3D. No Python required.
- Full program (from source) — the complete Treetrack tool for creating and editing segmentations, run from the Python source.

---

# vWell Viewer (Standalone App)

A standalone viewer for inspecting a scan and its completed segmentation in 3D. You select a scan file, and the viewer shows the scan slices together with any vessel segmentation that was previously created and saved in the main Treetrack tool. The viewer only *displays* segmentations — it does not create or edit them.

You do not need Python or any packages installed to run it.

## Installing the viewer

### Windows

1. Download `Treetrack_Setup.exe`.
2. Double-click it and follow the prompts to choose an install location and create shortcuts.
3. Launch from the **Start Menu** or **desktop shortcut**.

### macOS

1. Download and double-click the `.dmg` to mount it.
2. Drag **vWell Viewer** into your **Applications** folder, then eject the `.dmg`.
3. Launch from **Applications** or **Launchpad**.

> **First launch:** because this is unsigned research software, Windows and macOS show a security warning the first time you run it. This is expected and only needs clearing once — see [Troubleshooting &amp; Security](#troubleshooting--security).

## What you need before you start

The viewer works with three things, all living **in the same folder**:

| File / folder | What it is | Required? |
|---|---|---|
| `YourScan.nii` (or `.nii.gz`) | The brain scan, in NIfTI format. This is the file you select. | **Yes** |
| `YourScan_LOPRF2.npz` | The saved segmentation (the vessels). Created in the main Treetrack tool with the **S** key. | Needed to see vessels |
| `YourScan_DATA_RF2/` | A folder of cached vWell data. Created automatically on first run. | Optional (speeds up loading) |

The companion files all start with the **same base name as the scan**. If you rename the scan, rename them to match, or the viewer won't find them.

### Example folder layout

```
Test Brain/
├── ExtractedCOW.nii            ←  select this one
├── ExtractedCOW_LOPRF2.npz     ←  the saved segmentation
└── ExtractedCOW_DATA_RF2/      ←  cached data (created/used automatically)
        ExtractedCOW_SERF2.npz
        ExtractedCOW_SERF2_DD0.pkl
        ...
```

> The material included in the "Segmentation Data" folder follows this file structure. The vessels you see come from a segmentation completed by our lab using main Treetrack program.

## Loading a scan

1. **Launch the viewer** (see [Installing the viewer](#installing-the-viewer) above).
2. A **file-selection dialog titled "Select Input Image"** opens. A console/terminal window may also appear showing progress messages — leave it open while you work.
3. **Navigate to your scan folder** and select the `.nii` (or `.nii.gz`) file, then click **Open**.
4. The viewer loads/prepares the vWell data. This may take a moment.
5. The viewer then looks for the saved segmentation (`..._LOPRF2.npz`) in the scan's folder:
   - **Found** → the 3D vessels load automatically and the console prints `Loaded existing segmentation.`
   - **Not found** → the console prints `No saved segmentation found.` and you'll see only the scan slices (no vessels).
6. The **3D viewer window** opens.

## Controls

A summary of these controls is shown on screen inside the viewer. Press **?** to hide or show it.

### Keyboard

| Key | Action |
|---|---|
| **Up / Down** | Move the slice up / down |
| **Left / Right** | Switch the displayed image (input scan → vWell mean → variance) |
| **x / y / z** | View the slice along the X, Y, or Z axis |
| **D** | Show / hide the slice plane |
| **V** | Show / hide the 3D vessels |
| **o** | Toggle vessel opacity (50% ↔ 100%) |
| **?** | Show / hide the on-screen command list |

### Mouse

| Action | Result |
|---|---|
| Left-drag (off the plane) | Rotate the view |
| Middle-drag (off the plane) | Pan / translate |
| Right-drag (off the plane) | Zoom |
| Middle-drag (on the plane) | Slide the slice plane through the volume |
| Right-drag (on the plane) | Adjust brightness / contrast |
| Scroll wheel | Zoom |

## Notes

- Keep the console/terminal window open while using the viewer; closing it closes the program.
- This build uses a fixed resampling factor of **2** (hence the `RF2` / `_LOPRF2` naming).
- Input scans should be NIfTI (`.nii` or `.nii.gz`).

---

# Troubleshooting &amp; Security

## Getting past the first-launch security warning

The viewer is unsigned research software, so the operating system warns you the first time you run it. This does **not** mean the app is unsafe — it only means the publisher isn't registered with Microsoft or Apple. You clear it once.

### Windows (SmartScreen)

1. When you run `Treetrack_Setup.exe`, a blue **"Windows protected your PC"** dialog appears.
2. Click the **More info** link (small text below the message).
3. Click the **Run anyway** button that appears.
4. If prompted by **User Account Control**, click **Yes**.

**Antivirus false positives:** self-contained executables are sometimes flagged incorrectly. If the installer or app disappears or won't run, open **Windows Security → Virus & threat protection → Protection history**, find the quarantined item, and choose **Restore** (or add an **exclusion** for the install folder). On managed machines, IT may need to approve it.

### macOS (Gatekeeper)

1. Double-click the app once. You'll get a *"developer cannot be verified"* message — click **Done** or **Cancel** to dismiss it. (This registers the attempt.)
2. Open **System Settings** → **Privacy & Security**, then scroll down to the **Security** section.
3. Find the line noting that the app *was blocked* and click **Open Anyway**. (This button only appears for about an hour after you tried to open the app; if it's gone, double-click the app again.)
4. Confirm with **Open Anyway** in the popup, then enter your password or Touch ID.

On older macOS you can instead **Control-click** (or right-click) the app in Finder and choose **Open**. If macOS says the app *"is damaged,"* that's the download quarantine flag — remove it in Terminal:

```bash
xattr -dr com.apple.quarantine "/Applications/vWell Viewer.app"
```

## Viewer troubleshooting

**No vessels appear.**
The `..._LOPRF2.npz` segmentation file isn't in the scan's folder, or its name doesn't match the scan's base name. Confirm the file is present and named correctly (e.g. `ExtractedCOW.nii` → `ExtractedCOW_LOPRF2.npz`). If the segmentation was never saved in the main tool, there's nothing to display.

**The dialog flashes and closes, or nothing opens.**
Check the console/terminal window for an error message. The most common causes are selecting a file that isn't a readable NIfTI image, or a missing/corrupt companion file.

**The slice looks washed out or too dark.**
Right-drag while hovering over the slice plane to adjust brightness/contrast, or use Left/Right to switch between the input, mean, and variance images.

**Vessels are hidden.**
Press **V** to toggle them back on, and **o** to switch between opaque and semi-transparent rendering.

---

# Running from Source (Full Program)

The full Treetrack program — including creating and editing segmentations — runs from the Python source. `tt_Main.py` is the full segmentation editor; `tt_Viewer.py` is the read-only viewer (the same tool packaged as the standalone app above).

## Installing the libraries

Treetrack runs on Python 3. Install the dependencies with:

```bash
pip install -r requirements.txt
```

> **Note:** This is only needed to run the `.py` source files. The packaged viewer bundles everything.

### Python packages

| Package | Purpose |
| --- | --- |
| `numpy`, `scipy`, `pandas` | Core scientific computing stack |
| `itk` | Medical image I/O and processing |
| `scikit-image` | Morphological operations (imported as `skimage.morphology`) |
| `vtk` | 3D visualization |
| `matplotlib` | Plotting |
| `numba` | JIT compilation for the fast graph / region-grow routines |
| `xlsxwriter` | Excel output engine for pandas (`pd.ExcelWriter(..., engine='xlsxwriter')`) — not imported directly, but required at runtime when saving results |

### Standard library

These ship with Python and need no installation: `os`, `pickle`, `time`, `collections`, `heapq`, `tkinter`.

`tkinter` (used for the file-open dialog) is included with Python on Windows and macOS. On Linux, if it's missing, install it via your system package manager:

```bash
sudo apt install python3-tk
```

## Running the program

From the project root, with your dependencies installed:

```bash
python tt_Main.py       # full segmentation editor
python tt_Viewer.py     # read-only viewer
```

Both open the same **"Select Input Image"** file dialog described under [Loading a scan](#loading-a-scan). The editor adds segmentation creation and editing on top of the viewer's display features.
