from matplotlib import pyplot as plt
import pandas as pd

fname1 = "D:\\3dsim\\singleChaosData\\probe.csv"
fname2 = "D:\\3dsim\\singleChaosData\\perturbations.csv"

df1 = pd.read_csv(fname1)
df2 = pd.read_csv(fname2)

# Probe
probeX = df1.loc[:, "x"]
probeY = df1.loc[:, "y"]
probeZ = df1.loc[:, "z"]

# Perturbations
px = df2.loc[:, "x"]
py = df2.loc[:, "y"]
pz = df2.loc[:, "z"]

# Plot
ax = plt.axes(projection="3d")
ax.scatter(probeX, probeY, probeZ)
# ax.scatter(149597870700, 0, 0)
ax.scatter(px, py, pz)

plt.show()