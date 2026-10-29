# Generate the photopeaks file to be used in the energy calibration


#path = '/home/kale-chen/Documents/PET/MDA_04112026/run1singles.ldat'
path = '/home/kale-chen/Documents/PET/MiniPET/singles/old/Run118_singles.ldat'
#path = 'data/testcoinc.bin'
#path = 'data/testtarget29.bin'
#path = 'data/joint_filter_testtarget29.bin'
row_size = 16

import numpy as np
import os
from tqdm import tqdm
import energycal

num_rows = os.path.getsize(path) // row_size
dtypes = np.dtype([('time', np.int64), 
                    ('energy', np.float32), 
                    ('id', np.int32)])
data = np.memmap(path, dtype=dtypes, mode='r', shape=(num_rows,))
print(data[0:3])
print(num_rows)

# Highest hit energy collection
chunk_size = 10_000_000  # or any size suitable for memory usage
totalevents = 0
t0 = data[0]['time']
energies = [[] for _ in range(256)]
for start in tqdm(range(0, len(data), chunk_size)):
    chunk = np.array(data[start:start + chunk_size])
    chunk_ids = chunk['id'] % 1024
    chunk_ids_adj = np.where(chunk_ids < 512, chunk_ids - 256, chunk_ids - 384)
    chunk['id'] = chunk_ids_adj
    times = chunk['time']
    id_vals = chunk['id']
    dt = np.diff(times, prepend=times[0])
    group_boundaries = np.where(dt >= 10000)[0]
    group_starts = np.append(0, group_boundaries)
    group_ends = np.append(group_boundaries, len(chunk))
    for start_idx, end_idx in zip(group_starts[:-1], group_starts[1:]):
        group_ids = id_vals[start_idx:end_idx]
        if group_ids.size == 0:
            continue
        group_energies = chunk['energy'][start_idx:end_idx]
        max_energy_idx = np.argmax(group_energies)
        max_energy_id = group_ids[max_energy_idx]
        max_energy = group_energies[max_energy_idx]
        energies[max_energy_id].append(max_energy)
        totalevents += 1
print(totalevents)
photopeaks = np.zeros(256)
for id in range(256):
    photopeaks[id] = energycal.fitphotopeak(energies[id])[0][1]
photopeaks.tofile(os.path.join(os.path.dirname(__file__), 'photopeaks.bin')) # Read this later for singles analysis
