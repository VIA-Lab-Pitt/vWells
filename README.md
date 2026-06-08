## Requirements

Treetrack runs on Python 3. To run the source code, install the dependencies with:

```bash
pip install -r requirements.txt
```

> **Note:** This is only needed to run the `.py` source files. The packaged executable bundles everything and needs no installation.

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
