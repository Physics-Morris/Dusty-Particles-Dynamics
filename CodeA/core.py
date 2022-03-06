import numpy
import math
import scipy
import scipy.interpolate
import random


class Simulation:

    # initialize the simulation
    def __init__(self, total_time: float, dt: float, temperature_i: float, temperature_n: float,
                 mass_n: float, mass_i: float, density_n: float, density_i: float, debye_length: float,
                 E_field: float, grid_size: float, boundary: dict) -> None:
        self.total_time = total_time
        self.dt = dt
        self.temperature_i = temperature_i
        self.temperature_n = temperature_n
        self.mass_n = mass_n
        self.density_n = density_n
        self.density_i = density_i
        self.mass_i = mass_i
        self.particles_pos = []
        self.particles_vel = []
        self.particles_r = []
        self.particles_rho = []
        self.particles_charge = []
        self.mass = []
        self.debye_length = debye_length
        self.E_field = E_field
        self.grid_size = grid_size
        self.boundary = boundary

    # add non-interacting particles
    def particles(self, pos: float, vel: float, r: float, rho: float, charge: float) -> None:
        for i in range(len(pos)):
            self.particles_pos.append(pos[i])
            self.particles_vel.append(vel[i])
            self.particles_r.append(r[i])
            self.particles_rho.append(rho[i])
            self.particles_charge.append(charge[i])

    def get_trajectory(self, i):
        return self.traj_pos[:, i, :], self.traj_vel[:, i, :]

    def get_trajectory_all(self):
        return self.traj_pos[:, :, :], self.traj_vel[:, :, :]

    # run the simulation
    def run(self) -> None:
        # initialize the variables
        self.particles_pos = numpy.array(self.particles_pos)
        self.particles_vel = numpy.array(self.particles_vel)
        self.particles_r = numpy.array(self.particles_r)
        num_part = len(self.particles_pos)
        self.traj_pos = numpy.zeros((int(self.total_time / self.dt), num_part, 3))
        self.traj_vel = numpy.zeros((int(self.total_time / self.dt), num_part, 3))

        # calculate mass of dusty particles
        for i in range(num_part):
            self.mass.append(self.particles_rho[i] * self.particles_r[i]**3 * 4/3 * math.pi)

        # run the simulation
        for i in range(int(self.total_time / self.dt)):
            self.move(num_part, step=i)

    # move particles
    def move(self, num_part, step):
        self.traj_pos[step, :] = self.particles_pos
        self.traj_vel[step, :] = self.particles_vel
        new_x = numpy.empty(self.particles_pos.shape)
        new_v = numpy.empty(self.particles_vel.shape)
        for i in range(num_part):
            x1, v1 = self.rk4(i, self.dt, self.particles_pos[i, :],
                              self.particles_vel[i, :])
            new_x[i, :] = x1
            new_v[i, :] = v1
        self.particles_pos = new_x
        self.particles_vel = new_v

    # check if the particle is in the boundary
    def boundary_check(self, i, x0):
        msg = 'ok'
        # xmin
        if x0[0] <= self.boundary['xmin'][1]:
            if self.boundary['xmin'][0] == 'reflect': msg = 'reflect_xmin'
            elif self.boundary['xmin'][0] == 'open': msg = 'open_xmin'
        # xmax
        if x0[0] >= self.boundary['xmax'][1]:
            if self.boundary['xmax'][0] == 'reflect': msg = 'reflect_xmax'
            elif self.boundary['xmax'][0] == 'open': msg = 'open_xmax'
        # ymin
        if x0[1] <= self.boundary['ymin'][1]:
            if self.boundary['ymin'][0] == 'reflect': msg = 'reflect_ymin'
            elif self.boundary['ymin'][0] == 'open': msg = 'open_ymax'
        # ymax
        if x0[1] >= self.boundary['ymax'][1]:
            if self.boundary['ymax'][0] == 'reflect': msg = 'reflect_ymax'
            elif self.boundary['ymax'][0] == 'open': msg = 'open_ymax'
        # zmin
        if x0[2] <= self.boundary['zmin'][1]:
            if self.boundary['zmin'][0] == 'reflect': msg = 'reflect_zmin'
            elif self.boundary['zmin'][0] == 'open': msg = 'open_zmin'
        # zmax
        if x0[2] >= self.boundary['zmax'][1]:
            if self.boundary['zmax'][0] == 'reflect': msg = 'reflect_zmax'
            elif self.boundary['zmax'][0] == 'open': msg = 'open_zmax'
        return msg

            
    # Ruggle-Kutta method
    def rk4(self, i, dt: float, x0: float, v0: float):
        msg = self.boundary_check(i, x0) 
        if msg == 'ok':
            k1 = dt * self.total_force(i, x0) / self.mass[i]
            g1 = v0 * dt
            msg = self.boundary_check(i, x0 + g1 / 2.0) 
            if msg == 'ok':
                k2 = dt * self.total_force(i, x0 + g1 / 2.0) / self.mass[i]
                g2 = dt * (v0 + k1 / 2.0)
                msg = self.boundary_check(i, x0 + g2 / 2.0) 
                if msg == 'ok':
                    k3 = dt * self.total_force(i, x0 + g2 / 2.0) / self.mass[i]
                    g3 = dt * (v0 + k2 / 2.0)
                    msg = self.boundary_check(i, x0 + g3)
                    if msg == 'ok':
                        k4 = dt * self.total_force(i, x0 + g3) / self.mass[i]
                        g4 = dt * (v0 + k3)
                        v1 = v0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
                        x1 = x0 + (g1 + 2.0 * g2 + 2.0 * g3 + g4) / 6.0
                        return x1, v1

        if msg == 'reflect_xmin' or 'reflect_xmax': 
            x1 = x0 + dt * v0 * numpy.array([-1, 1, 1])
            v1 = v0 * numpy.array([-1, 1, 1])
        if msg == 'reflect_ymin' or 'reflect_ymax': 
            x1 = x0 + dt * v0 * numpy.array([1, -1, 1])
            v1 = v0 * numpy.array([1, -1, 1])
        if msg == 'reflect_zmin' or 'reflect_zmax': 
            x1 = x0 + dt * v0 * numpy.array([1, 1, -1])
            v1 = v0 * numpy.array([1, 1, -1])
        if msg == 'open_xmin' or 'open_xmax': 
            x1 = x0
            v1 = numpy.array([0, 0, 0])
        if msg == 'open_ymin' or 'open_ymax': 
            x1 = x0
            v1 = numpy.array([0, 0, 0])
        if msg == 'open_zmin' or 'open_zmax': 
            x1 = x0
            v1 = numpy.array([0, 0, 0])
        return x1, v1

    # calculate the total force 
    def total_force(self, i, x0):
        # create meshgrid
        xmin, xmax = self.boundary['xmin'][1], self.boundary['xmax'][1]
        ymin, ymax = self.boundary['ymin'][1], self.boundary['ymax'][1]
        zmin, zmax = self.boundary['zmin'][1], self.boundary['zmax'][1]
        x = numpy.linspace(xmin, xmax, int((xmax - xmin) / self.grid_size[0]))
        y = numpy.linspace(ymin, ymax, int((ymax - ymin) / self.grid_size[1]))
        z = numpy.linspace(zmin, zmax, int((zmax - zmin) / self.grid_size[2]))
        points = (x, y, z)
        total = numpy.zeros(3)
        total += self.gravity(i)
        total += self.electrostatic(i, x0, points)
        return total

    # calculate gravity force Fg = 4/3 * pi * r^3 * rho * g
    def gravity(self, i):
        Fg = numpy.array([0, 0, -4/3 * math.pi * self.particles_r[i]**3 * self.particles_rho[i] * 9.8])
        return Fg

    # calculate electrostatic force Fe = Qd * E * ( 1 + (rd / ld)^2 / 3*(1 + rd / ld) )
    def electrostatic(self, i, x0, points):
        Qd = self.particles_charge[i]
        rd = self.particles_r[i]
        ld = self.debye_length
        # interpolate the electric field
        Ex_field = self.E_field[0]
        Ey_field = self.E_field[1]
        Ez_field = self.E_field[2]
        Ex = scipy.interpolate.interpn(points, Ex_field, x0, method='linear')
        Ey = scipy.interpolate.interpn(points, Ey_field, x0, method='linear')
        Ez = scipy.interpolate.interpn(points, Ez_field, x0, method='linear')
        E_interp = numpy.array([Ex[0], Ey[0], Ez[0]])
        Fe = Qd * E_interp * (1 + (rd / ld)**2 / 3 * (1 + rd / ld))
        return Fe

    # calculate thermophoresis force FT = -32/15 * (pi * m_n / 8 / kb / Tn)**0.5 * r_d^2 * k_tr * del_Tn
    def thermophoresis(self, particle: list) -> list:
        # calculate gradient T
        def gradientT(self):
            pass
        pass

    # calculate neutral drag force Fn = - 4/3 * pi * r_d^2 * n_n * m_n * v_thn * (v_d - v_n) 
    def neutral_drag(self, particle: list) -> list:
        pass

    # calculate ion drag force Fi = F_i-coll + F_i-or = ni * mi * vi * vs * pi * bc^2 +
    #                                                   ni * mi * vi * vs * 4 * pi * bn_pi/2^2 * gamma
    def ion_drag(self, particle: list) -> list:
        pass