import numpy as np

l0 = -4 # Initial coordinate
l1 = 4 # Final coordinate
T = 20 # Time of path evolution

b0 = 3

pointsNumber = 1000 # Number of points in the path

header = ','.join(["t", "l", "theta", "phi", "fx", "fy", "fz", "upx", "upy", "upz"])

def func_l(t):
    return l0 + (l1 - l0) * t / T

def b(l):
    return b0*np.exp(-10*(l/l0)**2)

def alpha(l):
    return np.pi-np.arctan(b(l)/l)

def fx(l):
    return np.sign(l)*np.cos(alpha(l))

def fy(l):
    return np.sign(l)*np.sin(alpha(l))

ls = np.linspace(l0, l1, pointsNumber)
ts = np.linspace(0, T, pointsNumber)
bs = b(ls)
alphas = alpha(ls)
alphas = [a if a > 0 else a + np.pi for a in alphas]

fys = fy(ls)
fxs = fx(ls)

# Generating path

path = []

for i in range(pointsNumber):
    t = ts[i]
    l = ls[i]
    theta = np.pi/2
    phi = 0
    fx = fxs[i]
    fy = fys[i]
    fz = 0
    upx = 0
    upy = 0
    upz = 1

    path.append(','.join([str(float(x)) for x in [t, l, theta, phi, fx, fy, fz, upx, upy, upz]]))

    # if t > 5: # Break after 5 seconds of path
        # break

# Saving path as a csv file

fileText = '\n'.join([header] + path)

with open('path_through.csv', 'w') as f:
    f.write(fileText)
