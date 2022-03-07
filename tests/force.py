# libraries
from tkinter import W
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import matplotlib as mpl

mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['font.size'] = 16
plt.rcParams['figure.figsize'] = [5.6, 4]
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1
colors = cm.get_cmap('Set1', 9)
 
kb = 1.38064852 * 1e-23
eps0 = 8.854187817e-12
ec = 1.60217662e-19

grid_size = [1e-2]*3
x = np.linspace(0, 10, int((10 - 0) / grid_size[0]))
y = np.linspace(0, 10, int((10 - 0) / grid_size[1]))
z = np.linspace(0, 10, int((10 - 0) / grid_size[2]))
Y, Z = np.meshgrid(y, z)
points = (x, y, z)

def Ey_func(y, z): return (y - 5) * 1e-3 / ec
def Ez_func(y, z): return 1 / (abs(y - 5) + 0.1) * 1e-3 / ec
def Ti_func(y, z): return (z + abs(y - 5)) * 1e-20
def Tn_func(y, z): return (z + abs(y - 5))* 1e-20

plt.contourf(y, z, Ey_func(Y, Z), 30, cmap='jet')
plt.xlabel('$y$')
plt.ylabel('$z$')
plt.title('$E_y$')
plt.colorbar()
plt.tight_layout()
plt.savefig('Ey.png')
plt.show()


plt.contourf(y, z, Ez_func(Y, Z), 30, cmap='jet')
plt.xlabel('$y$')
plt.ylabel('$z$')
plt.title('$E_z$')
plt.colorbar()
plt.tight_layout()
plt.savefig('Ez.png')
plt.show()

plt.contourf(y, z, Ti_func(Y, Z), 30, cmap='jet')
plt.xlabel('$y$')
plt.ylabel('$z$')
plt.title('$T_i$')
plt.colorbar()
plt.tight_layout()
plt.savefig('Ti.png')
plt.show()

plt.contourf(y, z, Tn_func(Y, Z), 30, cmap='jet')
plt.xlabel('$y$')
plt.ylabel('$z$')
plt.title('$T_n$')
plt.colorbar()
plt.tight_layout()
plt.savefig('Tn.png')
plt.show()