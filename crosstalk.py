import numpy as np
import os
import energycal
from tqdm import tqdm

leakagematrix = np.fromfile(os.path.join(os.path.dirname(__file__), 'leakagematrix.bin'), dtype=np.float32).reshape(64, 64)

def correctenergies(crystal1, crystal2, measuredenergy1, measuredenergy2):
    '''
    Corrects the energies of the two crystals as per the leakagematrix.
    Returns corrected energies.
    '''
    #correctedenergy1 = measuredenergy1 - leakagematrix[crystal1, crystal2] * measuredenergy2
    #correctedenergy2 = measuredenergy2 - leakagematrix[crystal2, crystal1] * measuredenergy1
    denominator = 1 - leakagematrix[crystal1, crystal2] * leakagematrix[crystal2, crystal1]
    correctedenergy1 = (measuredenergy1 - leakagematrix[crystal2, crystal1] * measuredenergy2) / denominator
    correctedenergy2 = (measuredenergy2 - leakagematrix[crystal1, crystal2] * measuredenergy1) / denominator
    return correctedenergy1, correctedenergy2

def obtainleakagematrix():
    path = os.path.join(os.path.dirname(__file__), 'data/testcoinc.bin')
    row_size = 16
    num_rows = os.path.getsize(path) // row_size
    dtypes = np.dtype([('time', np.int64), 
                        ('energy', np.float32), 
                        ('id', np.int32)])
    data = np.memmap(path, dtype=dtypes, mode='r', shape=(num_rows,))
    print(data[0:3])
    print(num_rows)
    photopeaks = np.fromfile(os.path.join(os.path.dirname(__file__), 'photopeaks.bin'))
    escalars = energycal.obtainscalars(photopeaks)
    
    atol = 10.0
    leakagematrix = np.zeros((64, 64)) # Models leakage from crystal at row -> col
    sourcecounts = np.zeros(64)
    for i in tqdm(range(0, num_rows, 1000000)):
        chunk = np.array(data[i:i+1000000])
        chunk = chunk[chunk['id'] > 127]
        chunk = chunk[chunk['id'] < 192]
        chunk['energy'] = energycal.convertenergy(chunk['energy'], escalars[chunk['id']])
        chunk['id'] = chunk['id'] % 64
        lasttime = chunk[0]['time']
        sidx, eidx = 0, 0
        for i in range(len(chunk)):
            if chunk[i]['time'] - lasttime > 10000:
                if np.any(np.isclose(chunk[sidx:eidx]['energy'], 511.0, atol = atol)):
                    centralid = chunk[sidx:eidx]['id'][np.where(np.isclose(chunk[sidx:eidx]['energy'], 511.0, atol = atol))[0][0]]
                    leakagematrix[centralid, chunk[sidx:eidx]['id']] += chunk[sidx:eidx]['energy']
                    sourcecounts[centralid] += 1
                sidx = eidx
                lasttime = chunk[i]['time']
            eidx += 1
    leakagematrix = leakagematrix / sourcecounts
    for i in range(64): # Normalize by central crystal
        #leakagematrix[:, i] = leakagematrix[:, i] / leakagematrix[i, i]
        leakagematrix[i] = leakagematrix[i] / leakagematrix[i, i]
    leakagematrix[np.isnan(leakagematrix)] = 0
    leakagematrix = leakagematrix.astype(np.float32)
    leakagematrix.tofile(os.path.join(os.path.dirname(__file__), 'leakagematrix.bin'))
    return leakagematrix

def plotleakagematrixrow(id):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    indices = {
        0: (4, 0),
        1: (4, 1),
        2: (7, 2),
        3: (5, 0),
        4: (5, 3),
        5: (5, 2),
        6: (4, 3),
        7: (7, 0),
        8: (6, 1),
        9: (7, 3),
        10: (5, 1),
        11: (6, 3),
        12: (4, 2),
        13: (6, 2),
        14: (6, 0),
        15: (7, 1),
        16: (3, 0),
        17: (3, 1),
        18: (2, 0),
        19: (2, 1),
        20: (0, 0),
        21: (1, 0),
        22: (0, 1),
        23: (1, 1),
        24: (3, 2),
        25: (1, 2),
        26: (2, 2),
        27: (4, 4),
        28: (0, 2),
        29: (3, 3),
        30: (0, 3),
        31: (1, 3),
        32: (2, 3),
        33: (3, 4),
        34: (2, 4),
        35: (0, 4),
        36: (1, 4),
        37: (0, 5),
        38: (5, 4),
        39: (1, 5),
        40: (2, 5),
        41: (3, 5),
        42: (1, 6),
        43: (0, 6),
        44: (0, 7),
        45: (3, 6),
        46: (1, 7),
        47: (2, 6),
        48: (3, 7),
        49: (2, 7),
        50: (6, 5),
        51: (6, 6),
        52: (7, 6),
        53: (4, 6),
        54: (5, 6),
        55: (6, 7),
        56: (7, 7),
        57: (7, 5),
        58: (7, 4),
        59: (4, 5),
        60: (5, 7),
        61: (5, 5),
        62: (6, 4),
        63: (4, 7),
        64: (11, 7),
        65: (11, 6),
        66: (10, 4),
        67: (10, 6),
        68: (11, 4),
        69: (9, 7),
        70: (9, 5),
        71: (8, 7),
        72: (9, 6),
        73: (11, 3),
        74: (10, 7),
        75: (10, 3),
        76: (11, 5),
        77: (8, 5),
        78: (10, 5),
        79: (8, 6),
        80: (12, 7),
        81: (14, 7),
        82: (13, 7),
        83: (15, 7),
        84: (13, 6),
        85: (12, 6),
        86: (14, 6),
        87: (15, 6),
        88: (12, 5),
        89: (14, 5),
        90: (13, 5),
        91: (9, 4),
        92: (15, 5),
        93: (15, 4),
        94: (14, 4),
        95: (12, 4),
        96: (13, 4),
        97: (13, 3),
        98: (14, 3),
        99: (12, 3),
        100: (15, 3),
        101: (15, 2),
        102: (8, 4),
        103: (14, 2),
        104: (15, 1),
        105: (13, 2),
        106: (14, 1),
        107: (12, 2),
        108: (15, 0),
        109: (13, 1),
        110: (14, 0),
        111: (12, 1),
        112: (13, 0),
        113: (12, 0),
        114: (8, 2),
        115: (10, 1),
        116: (8, 1),
        117: (11, 1),
        118: (9, 1),
        119: (9, 0),
        120: (8, 0),
        121: (9, 2),
        122: (9, 3),
        123: (11, 2),
        124: (10, 0),
        125: (10, 2),
        126: (8, 3),
        127: (11, 0)
    }
    def get_color(value):
        return plt.cm.jet(value)
    row = leakagematrix[:,id]
    print(row)
    maxval = np.max(row)
    fig, ax = plt.subplots(figsize=(8, 8))
    for i in range(64):
        y, x = indices[i]
        y_ = 15 - y
        rect = patches.Rectangle((x, y_), 1, 1, linewidth=0, facecolor=get_color(row[i] / maxval), alpha=0.6)
        ax.add_patch(rect)
        ax.text(x + 0.5, y_ + 0.5, f'{row[i]:.2f}', fontsize=8, ha='center', va='center')
    ax.axis('off')
    ax.set_xlim(0, 8)
    ax.set_ylim(8, 16)
    ax.set_aspect('equal')
    ax.set_title(f'Leakage matrix row {id}')
    plt.show()

if __name__ == '__main__':
    #leakagematrix = obtainleakagematrix()
    #print(leakagematrix[0:6,0:6])
    plotleakagematrixrow(29)
