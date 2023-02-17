from matplotlib import pyplot as plt
import pandas as pd

# Desktop
fname = "D:\\3dsim\\unitCircle.csv"

# Laptop
# fname = "C:\\Users\\natha\\eclipse-workspace\\3dSim\\unitCircle.csv"

df = pd.read_csv(fname)
# print(df.head())
x = df.iloc[:, 0]
y = df.iloc[:, 1]
z = df.iloc[:, 2]

# Sphere Plot

# Orbit Plot
ax = plt.axes(projection="3d")
ax.scatter(x, y, z)

plt.show()