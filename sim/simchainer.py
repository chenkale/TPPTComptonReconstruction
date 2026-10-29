# For running many simulations and chaining them together

import numpy as np
import os
import subprocess
import pandas as pd
import time
'''
Config
'''
njobs = 1000

geant4_dir = "/home/kale-chen/Documents/Geant4/install/bin/geant4.sh"
working_dir = "/home/kale-chen/Documents/Geant4/TPPTsim/temp"

sim_dir = "/home/kale-chen/Documents/Geant4/TPPTsim/build/sim"
builder_dir = "/home/kale-chen/Documents/Geant4/TPPTbuilder/build/builder"
coinc_dir = "/home/kale-chen/Documents/Geant4/TPPTcoinc/build/coinc"

sim_out_dir = os.path.join(working_dir, "sout.bin")
builder_out_dir = os.path.join(working_dir, "bout.bin")
coinc_out_dir = os.path.join(working_dir, "cout.txt")
gen_pos_dir = os.path.join(working_dir, "positions.txt")

out_dir = os.path.join(working_dir, 'finalcoincs.bin')
gen_hist_dir = os.path.join(working_dir, 'genhist.bin')
'''
Script
'''
start_time = time.time()


# Load Geant4
def load_geant4():
    command = f"""
    source "{geant4_dir}"
    env -0
    """
    result = subprocess.run(
        ["bash", "-c", command],
        capture_output=True,
        check=True
    )
    for entry in result.stdout.split(b"\0"):
        if b"=" in entry:
            key, value = entry.split(b"=", 1)
            os.environ[key.decode()] = value.decode()
load_geant4()

def run_sim(idoffset, seed):
    subprocess.run([sim_dir, "-i", str(idoffset), "-s", str(seed)], check=True)
def run_builder():
    subprocess.run([builder_dir], check=True)
def run_coinc():
    subprocess.run([coinc_dir], check=True)

if os.path.exists(out_dir):
    os.remove(out_dir)
if os.path.exists(gen_hist_dir):
    os.remove(gen_hist_dir)
genimg = np.zeros((320, 320, 110), dtype = np.int16)
genimg.tofile(gen_hist_dir)
for i in range(njobs):
    run_sim(i*1000000, i)
    run_builder()
    run_coinc()
    coincs = pd.read_csv(coinc_out_dir, sep=',', header=None)
    positions = pd.read_csv(gen_pos_dir, sep=' ', header=None, usecols=[1,2,3,4,5,6,7,8,9])
    with open(out_dir, 'ab') as f:
        for row in coincs.iterrows():
            np.asarray(row[1][2:4], dtype = np.int16).tofile(f)
            index = int(row[1][0])
            np.asarray(positions.iloc[index - i*1000000], dtype = np.float32).tofile(f)
    genimg = np.fromfile(gen_hist_dir, dtype=np.int16).reshape((320, 320, 110))
    inds = positions.iloc[:, [0, 1, 2]].to_numpy(dtype=np.float32)
    inds = np.floor(inds + np.array([160, 160, 55], dtype=np.float32)).astype(np.int32)
    valid_x = (inds[:, 0] >= 0) & (inds[:, 0] < 320)
    valid_y = (inds[:, 1] >= 0) & (inds[:, 1] < 320)
    valid_z = (inds[:, 2] >= 0) & (inds[:, 2] < 110)
    valid_mask = valid_x & valid_y & valid_z
    inds_valid = inds[valid_mask]

    if inds_valid.shape[0] > 0:
        flat_inds = np.ravel_multi_index(
            (inds_valid[:, 0], inds_valid[:, 1], inds_valid[:, 2]),
            dims=genimg.shape,
            mode='clip'
        )
        updates = np.bincount(flat_inds, minlength=genimg.size)
        if updates.any():
            genimg += updates.reshape(genimg.shape)
    genimg.tofile(gen_hist_dir)

end_time = time.time()
print(f"Total time: {round(end_time - start_time, 1)} s")