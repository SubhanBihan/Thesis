import numpy as np
import tkinter as tk
from tkinter import filedialog

# 1. Create a hidden root window so only the file dialog appears
root = tk.Tk()
root.withdraw()

# 2. Open the file browser
file_path = filedialog.askopenfilename(
    title="Select your .dat file",
    filetypes=[("Data files", "*.dat"), ("All files", "*.*")]
)

# 3. Process the file if the user didn't hit 'Cancel'
if file_path:
    # Adjust dtype to 'float64' if your data is double-precision
    data = np.fromfile(file_path, dtype='float32')
    
    print(f"Successfully loaded: {file_path}")
    print(f"Data preview:\n{data[:10]}...\n{data[9991:10000]}") # Shows the first 10 elements
else:
    print("No file selected.")