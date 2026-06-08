# Packages required to run the Treetrack source code.
# Install with:  pip install -r requirements.txt
#
# (Not needed for the packaged executable — only for running the .py files.)

# Core scientific stack
numpy
scipy
pandas

# Medical image I/O and processing
itk
scikit-image          # imported as skimage.morphology

# 3D visualization and plotting
vtk
matplotlib

# JIT compilation for the fast graph/region-grow routines
numba

# Excel output engine used by pandas (pd.ExcelWriter(..., engine='xlsxwriter'))
# Not imported directly, but required at runtime when saving results.
xlsxwriter

# --- Standard library (no install needed) ---
# os, pickle, time, collections, heapq, tkinter
#
# tkinter is used for the file-open dialog. It ships with Python on Windows
# and macOS. On Linux, if it's missing, install it via the system package
# manager, e.g.:  sudo apt install python3-tk
