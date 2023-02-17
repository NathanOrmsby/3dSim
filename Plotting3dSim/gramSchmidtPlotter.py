from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

# Desktop
fname1 = "D:\\3dsim\\gramSchmidtInitial.csv"
fname2 = "D:\\3dsim\\gramSchmidtFinal.csv"

# Laptop
# fname = "C:\\Users\\natha\\eclipse-workspace\\3dSim\\unitCircle.csv"

df1 = pd.read_csv(fname1)
df2 = pd.read_csv(fname2)
# print(df.head())

# Initial File
xi = df1.iloc[:, 0]
zi = df1.iloc[:, 1]

# Final File
xf = df2.iloc[:, 0]
zf = df2.iloc[:, 1]

# Plot initial vectors
plt.rcParams["figure.figsize"] = [12.8, 7.2]
plt.rcParams["figure.autolayout"] = True
soa = np.array([[0, 0, xi[0], zi[0]], [0, 0, xi[1], zi[1]], [0, 0, xi[2], zi[2]]])
X, Y, U, V = zip(*soa)
plt.figure()
ax = plt.gca()
ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1, color=['red', 'green', 'yellow'])
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
plt.draw()
plt.show()

# Plot final vectors
plt.rcParams["figure.figsize"] = [12.8, 7.2]
plt.rcParams["figure.autolayout"] = True
soa = np.array([[0, 0, xf[0], zf[0]], [0, 0, xf[1], zf[1]], [0, 0, xf[2], zf[2]]])
X, Y, U, V = zip(*soa)
plt.figure()
ax = plt.gca()
ax.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=1, color=['red', 'green', 'yellow'])
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
plt.draw()
plt.show()


