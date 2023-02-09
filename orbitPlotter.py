import matplotlib
import os
import pandas as pd

# Desktop
# fname = "D:\\3dsim\data.csv"

# Laptop
fname = "C:\\Users\\natha\\eclipse-workspace\\3dSim\\data.csv"


df = pd.read_csv(fname)
print(df.head())
earthX = df.loc[:, "ex"]
earthY = df.loc[:, "ey"]
earthZ = df.loc[:, "ez"]
print(earthX)
sunX = df.loc[:, "sx"]
sunY = df.loc[:, "sy"]
sunZ = df.loc[:, "sz"]
jwX = df.loc[:, "jwx"]
jwY = df.loc[:, "jwy"]
jwZ = df.loc[:, "jwz"]

# Sphere Plot
from matplotlib import pyplot as plt

# Orbit Plot
ax = plt.axes(projection="3d")
# ax.scatter(earthX, earthY, earthZ)
ax.scatter(149597870700, 0, 0)
ax.plot(earthX, earthY, earthZ)
ax.plot(jwX, jwY, jwZ)
ax.scatter(earthX.iloc[-1], earthY.iloc[-1], earthZ.iloc[-1])
ax.scatter(jwX.iloc[-1], jwY.iloc[-1], jwZ.iloc[-1])

plt.show()


