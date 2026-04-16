import os
import numpy as np
from tkinter import Tk, filedialog
from pathlib import Path

def convert_binary_data():
    # Hide the main Tkinter window
    root = Tk()
    root.withdraw()

    # 1. Ask for the source directory
    print("Please select the source parent folder...")
    src_root = filedialog.askdirectory(title="Select Source Folder")
    
    if not src_root:
        print("No folder selected. Exiting.")
        return

    src_path = Path(src_root)
    # 2. Create the 'alt' parent folder name
    dst_path = src_path.parent / f"{src_path.name}_f64"
    
    print(f"Source: {src_path}")
    print(f"Destination: {dst_path}")
    print("Processing...")

    count = 0
    # 3. Recursively find all files
    for file in src_path.rglob("*.dat"):
        # Determine the relative path to recreate the tree
        rel_path = file.relative_to(src_path)
        target_file = dst_path / rel_path
        
        # Create subdirectories in the destination if they don't exist
        target_file.parent.mkdir(parents=True, exist_ok=True)

        try:
            # 4. Read f32 data
            # We use NumPy's fromfile which is extremely fast for binary data
            data_f32 = np.fromfile(file, dtype=np.float32)

            # 5. Convert to f64
            data_f64 = data_f32.astype(np.float64)

            # 6. Write as f64 binary
            data_f64.tofile(target_file)
            count += 1
            print(f"Converted: {rel_path}")
            
        except Exception as e:
            print(f"Failed to convert {file}: {e}")

    print(f"\nFinished! Processed {count} files.")
    print(f"Your converted data is located at: {dst_path}")

if __name__ == "__main__":
    convert_binary_data()