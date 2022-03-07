import sys
import matplotlib.pyplot as plt
  
kb = 1.38064852 * 1e-23

# setting path
sys.path.append('../')

# import the module
from CodeA import core
import numpy as np

bd = {'xmin': ['open', 0.],
      'xmax': ['open', 10.],
      'ymin': ['open', 0.],
      'ymax': ['open', 10.],
      'zmin': ['open', 0.],
      'zmax': ['open', 10.]}
grid_size = [0.1, 0.1, 0.1]

#########################################################
# create simulation grid
x = np.linspace(0, 10, int((10 - 0) / grid_size[0]))
y = np.linspace(0, 10, int((10 - 0) / grid_size[1]))
z = np.linspace(0, 10, int((10 - 0) / grid_size[2]))
points = (x, y, z)

# create temperature data
def Ti_func(x, y, z): return 0*x + 0*y + 0*z
def Tn_func(x, y, z): return 0*x + 0*y + 0*z
Ti = Ti_func(*np.meshgrid(*points, indexing='ij'))
Tn = Tn_func(*np.meshgrid(*points, indexing='ij'))

# create E field data
def Ex_func(x, y, z): return 0*x + 0*y + 0*z
def Ey_func(x, y, z): return y - 5
def Ez_func(x, y, z): return 1 / (abs(y - 5) + 0.1)
Ex_field = Ex_func(*np.meshgrid(*points, indexing='ij'))
Ey_field = Ey_func(*np.meshgrid(*points, indexing='ij'))
Ez_field = Ez_func(*np.meshgrid(*points, indexing='ij'))
#########################################################


sim = core.Simulation(total_time=2, dt=1e-2, temperature_i=Ti, temperature_n=Tn,
                      mass_n=1, mass_i=1, density_n=1, density_i=1, debye_length=1,
                      E_field=[Ex_field, Ey_field, Ez_field], conductivity=1,
                      grid_size=grid_size, boundary=bd)

pos, vel = [], []
numpart = 16
for i in np.linspace(1, 9, numpart):
    pos.append([5, i, 9])
    vel.append([0, 0, 0])
sim.particles(pos=pos, vel=vel, r=[1]*numpart, rho=[1]*numpart, charge=[1]*numpart)

sim.run()


# plot the results
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

for i in range(numpart):
    traj_pos, traj_vel = sim.get_trajectory(i)
    plt.plot(traj_pos[:, 1], traj_pos[:, 2], '-')

plt.xlabel('$y$')
plt.ylabel('$z$')
plt.title('$F_g + F_E$')
plt.tight_layout()
plt.savefig('FE.png', dpi=500)
plt.show()