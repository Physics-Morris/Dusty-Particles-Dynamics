import sys
import matplotlib.pyplot as plt
  
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
x = np.linspace(0, 10, int((10 - 0) / grid_size[0]))
y = np.linspace(0, 10, int((10 - 0) / grid_size[1]))
z = np.linspace(0, 10, int((10 - 0) / grid_size[2]))
points = (x, y, z)

Ti = np.zeros((10, 10, 10))
Tn = np.zeros((10, 10, 10))

def Ex_func(x, y, z): return 0*x + 0*y + 0*z
def Ey_func(x, y, z): return y - 5
def Ez_func(x, y, z): return 1 / (abs(y - 5) + 0.1)
Ex_field = Ex_func(*np.meshgrid(*points, indexing='ij'))
Ey_field = Ey_func(*np.meshgrid(*points, indexing='ij'))
Ez_field = Ez_func(*np.meshgrid(*points, indexing='ij'))


sim = core.Simulation(total_time=2, dt=1e-2, temperature_i=Ti, temperature_n=Tn,
                      mass_n=1, mass_i=1, density_n=1, density_i=1, debye_length=1,
                      E_field=[Ex_field, Ey_field, Ez_field], grid_size=grid_size, boundary=bd)

pos, vel = [], []
numpart = 16
for i in np.linspace(1, 9, numpart):
    pos.append([5, i, 9])
    vel.append([0, 0, 0])
sim.particles(pos=pos, vel=vel, r=[1]*numpart, rho=[1]*numpart, charge=[1]*numpart)

sim.run()

for i in range(numpart):
    traj_pos, traj_vel = sim.get_trajectory(i)
    plt.plot(traj_pos[:, 1], traj_pos[:, 2], '.-')
plt.show()