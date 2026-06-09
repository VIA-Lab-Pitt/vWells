<p align="center">
  <img src="assets/treetrack-logo-light.svg" alt="Treetrack" width="500">
</p>

## Requirements

Treetrack runs on Python 3. To run the source code, install the dependencies with:

```bash
pip install -r requirements.txt
```

> **Note:** This is only needed to run the `.py` source files. The packaged executable bundles everything.

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

A standalone viewer for inspecting a brain scan and its segmented vessels in 3D. You select a scan file, and the viewer shows the scan slices together with any vessel segmentation that was previously created and saved in the main Treetrack tool.

This manual covers the executable build — you do **not** need Python or any packages installed to run it.

---

## 1. What you need before you start

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

> The vessels you see come from a segmentation made earlier in the main Treetrack program. The viewer only *displays* segmentations — it does not create or edit them.

---

## 2. Loading a scan

1. **Launch the executable** (double-click it).
2. A **file-selection dialog titled "Select Input Image"** opens. Depending on the build, a console/terminal window may also appear showing progress messages — leave it open while you work.
3. **Navigate to your scan folder** and select the `.nii` (or `.nii.gz`) file, then click **Open**.
4. The viewer prepares the vWell data:
   - **First time on a scan:** it computes the vWell structure from scratch. This can take a moment; progress prints in the console. When it finishes it saves a `..._DATA_RF2` folder so future loads are fast. Don't close the window during this step.
   - **After that:** it reuses the cached `..._DATA_RF2` folder and loads quickly.
5. The viewer then looks for the saved segmentation (`..._LOPRF2.npz`) in the scan's folder:
   - **Found** → the 3D vessels load automatically and the console prints `Loaded existing segmentation.`
   - **Not found** → the console prints `No saved segmentation found.` and you'll see only the scan slices (no vessels).
6. The **3D viewer window** opens.

---

## 3. Controls

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

---

## 4. Troubleshooting

**No vessels appear.**
The `..._LOPRF2.npz` segmentation file isn't in the scan's folder, or its name doesn't match the scan's base name. Confirm the file is present and named correctly (e.g. `ExtractedCOW.nii` → `ExtractedCOW_LOPRF2.npz`). If the segmentation was never saved in the main tool, there's nothing to display.

**The first load takes a moment.**
That's expected the first time — it's computing the vWell structure. It caches the result in `..._DATA_RF2`, so subsequent loads of the same scan are faster. To force a fresh recompute, delete that folder.

**The dialog flashes and closes, or nothing opens.**
Check the console/terminal window for an error message. The most common causes are selecting a file that isn't a readable NIfTI image, or a missing/corrupt companion file.

**The slice looks washed out or too dark.**
Right-drag while hovering over the slice plane to adjust brightness/contrast, or use Left/Right to switch between the input, mean, and variance images.

**Vessels are hidden.**
Press **V** to toggle them back on, and **o** to change opacity to 100%.

---

## 5. Notes

- Keep the console/terminal window open while using the viewer; closing it closes the program.
- This build uses a fixed resampling factor of **2** (hence the `RF2` / `_LOPRF2` naming).
- Input scans should be NIfTI (`.nii` or `.nii.gz`).
