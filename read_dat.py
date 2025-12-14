import h5py

with h5py.File("libra_data/mem_data.hdf", 'r') as f:
    print(list(f.keys()))

with h5py.File("libra_data/mem_data.hdf", 'r') as f:
    #q = list(f["q/data"][args.itraj * args.iskip + args.istart, 0, i:])
    Etot = f["Etot_ave/data"][:]
    t = f["time/data"][:]
    q = f["q/data"][:]
    dt = f["timestep/data"][:]

print(type(Etot))
print(Etot.shape)

print(type(t))
print(t.shape)

print(type(q))
print(q.shape)

print(type(dt))
print(dt.shape)

print(dt)
print(t)


