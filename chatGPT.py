import h5py
import numpy as np

def visit_datasets(name, obj):
    if isinstance(obj, h5py.Dataset):
        dtype = obj.dtype
        shape = obj.shape

        # Number of elements in dataset
        n_elements = np.prod(shape) if shape else 1

        # Size of each element in bytes
        itemsize = dtype.itemsize

        # Total size in bytes
        size_bytes = n_elements * itemsize

        print(f"Dataset: {name}")
        print(f"  dtype: {dtype}")
        print(f"  shape: {shape}")
        print(f"  total size: {size_bytes/1024**2:.2f} MB ({size_bytes:,} bytes)")
        print()

def summarize_hdf5(filename):
    with h5py.File(filename, "r") as f:
        f.visititems(visit_datasets)

if __name__ == "__main__":
    hdf_file = "mem_data.hdf"  # replace with your filename
    summarize_hdf5(hdf_file)
