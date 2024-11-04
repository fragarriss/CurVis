import numpy as np
from pathlib import Path

l0 = 3
T = 60 # Time of path evolution

pointsNumber = 1000 # Number of points in the path

header = ','.join(["t", "l", "theta", "phi", "fx", "fy", "fz", "upx", "upy", "upz"])

def phifunc(t):
    return 2*np.pi*t/T

ts = np.linspace(0, T, pointsNumber)
phis = [phifunc(t) for t in ts]

# Generating path

path = []

for i in range(pointsNumber):
    t = ts[i]
    l = l0
    theta = np.pi/2
    phi = phis[i]
    fx = -1
    fy = 0
    fz = 0
    upx = 0
    upy = 0
    upz = 1

    path.append(','.join([str(float(x)) for x in [t, l, theta, phi, fx, fy, fz, upx, upy, upz]]))

    # if t > 5: # Break after 5 seconds of path
        # break

# Saving path as a csv file

fileText = '\n'.join([header] + path)

filePath = Path(__file__).parent / 'path_orbit.csv'
with open(filePath, 'w') as f:
    f.write(fileText)

print(filePath)
