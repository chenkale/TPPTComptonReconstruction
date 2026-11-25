# Preprocess gammaPath.text to remove non interacting events
# based on what events are present in bout.bin

# Optionally, filter for events with detector hits <= 2
limit2 = True


import numpy as np
import pandas as pd
from tqdm import tqdm
import struct
import matplotlib.pyplot as plt

gpath = '/home/kale-chen/Documents/Geant4/TPPTsim/temp/gammaPath.txt'
bpath = '/home/kale-chen/Documents/Geant4/TPPTsim/temp/bout.bin'

# Read bout.bin, singles-like file containing hit information as seen by
# simulated detector, this is like what we would see in our real data.

# Final bout np array contains: 
# (ID of scintillator, event ID (used for grouping), energy in MeV)
rows = []
i_scint = None
with open(bpath, "rb") as f:
    while True:
        marker = f.read(1)
        if not marker:
            break
        if marker == b"\xee":
            i_scint = struct.unpack("<i", f.read(4))[0]
        elif marker == b"\xff":
            event_id, time_ns, energy_mev = struct.unpack("<qdd", f.read(24))
            rows.append((i_scint, event_id, time_ns, energy_mev))
        else:
            raise ValueError(f"Unexpected marker {marker!r} at offset {f.tell()-1}")
bout = pd.DataFrame(rows, columns=["iScint", "eventId", "time_ns", "energy_MeV"])
bout.sort_values(by='eventId', inplace=True)
bout = np.array(bout)[:, [0, 1, 3]]
bout[:, 2] = bout[:, 2] * 1e6 # Convert energy to eV and save as int32 for simplicity
bout = bout.astype(np.int32) # 
print(bout.shape)

if not limit2:
    eventids = np.unique(bout[:, 1])
else:
    eventids = []
    curreventid = bout[0, 1]
    currhits = 0
    for row in tqdm(bout):
        if not np.isclose(row[1], curreventid): # End of event, apply selection criteria
            if currhits <= 2:
                eventids.append(curreventid)
            curreventid = row[1]
            currhits = 0
        else:
            currhits += 1
print(len(eventids))

prelength = 0
postlength = 0
kept = pd.DataFrame()
for chunk in tqdm(pd.read_csv(gpath, sep = ' ', chunksize=1000000)):
    prelength += chunk.shape[0]
    kept = pd.concat([kept, chunk[chunk['eventId'].isin(eventids)]])
    postlength = kept.shape[0]
print(f"Prelength: {prelength}, Postlength: {postlength}, Removed: {prelength - postlength}")
kept.to_csv(gpath.replace('.txt', 'filtered.txt'), sep = ' ', index=False)