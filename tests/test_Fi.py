import sys
import matplotlib.pyplot as plt
  
kb = 1.38064852 * 1e-23
eps0 = 8.854187817e-12
ec = 1.60217662e-19

# setting path
sys.path.append('../')

# import the module
from CodeA import core
import numpy as np

bd = {'xmin': ['reflect', 0.],
      'xmax': ['reflect', 10.],
      'ymin': ['reflect', 0.],
      'ymax': ['reflect', 10.],
      'zmin': ['open', 0.],
      'zmax': ['reflect', 10.]}
grid_size = [0.1, 0.1, 0.1]

#########################################################
# create simulation grid
x = np.linspace(0, 10, int((10 - 0) / grid_size[0]))
y = np.linspace(0, 10, int((10 - 0) / grid_size[1]))
z = np.linspace(0, 10, int((10 - 0) / grid_size[2]))
points = (x, y, z)

# create temperature data
def Ti_func(x, y, z): return (z + abs(y - 5)) * 1e-20
def Tn_func(x, y, z): return (z + abs(y - 5))* 1e-20
Ti = Ti_func(*np.meshgrid(*points, indexing='ij'))
Tn = Tn_func(*np.meshgrid(*points, indexing='ij'))

# create E field data
def Ex_func(x, y, z): return 0*x + 0*y + 0*z
def Ey_func(x, y, z): return (y - 5) * 1e-3 / ec
def Ez_func(x, y, z): return 1 / (abs(y - 5) + 0.1) * 1e-3 / ec
Ex_field = Ex_func(*np.meshgrid(*points, indexing='ij'))
Ey_field = Ey_func(*np.meshgrid(*points, indexing='ij'))
Ez_field = Ez_func(*np.meshgrid(*points, indexing='ij'))

# create v_thermal_n
def v_thn(x, y, z): return 0*x + 0*y + 0*z + 1
v_thn = v_thn(*np.meshgrid(*points, indexing='ij'))

# create v_neutral
def v_n_funcx(x, y, z): return 0*x + 0*y + 0*z
def v_n_funcy(x, y, z): return 0*x + 0*y + 0*z
def v_n_funcz(x, y, z): return 0*x + 0*y + 0*z
v_nx = v_n_funcx(*np.meshgrid(*points, indexing='ij'))
v_ny = v_n_funcy(*np.meshgrid(*points, indexing='ij'))
v_nz = v_n_funcz(*np.meshgrid(*points, indexing='ij'))
v_n = [v_nx, v_ny, v_nz]

# create v_ion
def v_i_funcx(x, y, z): return 0*x + 0*y + 0*z + 0.e0
def v_i_funcy(x, y, z): return 0*x + 0*y + 0*z + 0.e0
def v_i_funcz(x, y, z): return 0*x + 0*y + 0*z + 1.e0
v_ix = v_i_funcx(*np.meshgrid(*points, indexing='ij'))
v_iy = v_i_funcy(*np.meshgrid(*points, indexing='ij'))
v_iz = v_i_funcz(*np.meshgrid(*points, indexing='ij'))
v_i = [v_ix, v_iy, v_iz]

# create ion and neutral density
def nn_func(x, y, z): return 0*x + 0*y + 0*z + 1
def ni_func(x, y, z): return 0*x + 0*y + 0*z + 1
n_n = nn_func(*np.meshgrid(*points, indexing='ij'))
n_i = ni_func(*np.meshgrid(*points, indexing='ij'))
#########################################################


sim = core.Simulation(total_time=5, dt=1e-2, temperature_i=Ti, temperature_n=Tn,
                      mass_n=1, mass_i=1, density_n=n_n, density_i=n_i, debye_length=1,
                      E_field=[Ex_field, Ey_field, Ez_field], conductivity=1, v_thn=v_thn,
                      v_n=v_n, v_i=v_i, grid_size=grid_size, boundary=bd)

pos, vel = [], []
numpart = 16
for i in np.linspace(1, 9, numpart):
    pos.append([5, i, 9])
    vel.append([0, 0, 0])
sim.particles(pos=pos, vel=vel, r=[1]*numpart, rho=[1]*numpart, charge=[-1000 * ec]*numpart)

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
plt.title('$F_g + F_E + F_T + F_n + F_i$')
plt.tight_layout()
plt.savefig('Fi_refl.png', dpi=500)
plt.show()