import numpy
import math
import random


class Simulation:

    # initialize the simulation
    def __init__(self, total_time: float, dt: float, temperature_i: float, temperature_n: float,
                 mass_n: float, density_n: float, density_i: float, mass_i: float) -> None:
        self.total_time = total_time
        self.dt = dt
        self.temperature_i = temperature_i
        self.temperature_n = temperature_n
        self.mass_n = mass_n
        self.density_n = density_n
        self.density_i = density_i
        self.mass_i = mass_i

    # calculate gradient T
    def gradientT(self):
        pass

    # add density of the simulation
    def density(self, density: float) -> None:
        pass

    # add non-interacting particles
    def particles(self, particles: list) -> None:
        pass

    # run the simulation
    def run(self) -> None:
        pass

    # move particles
    def move(self, particle: list) -> list:
        pass

    # calculate gravity force Fg = 4/3 * pi * r^3 * rho * g
    def gravity(self, particle: list) -> list:
        pass

    # calculate electrostatic force Fe = Qd * E * ( 1 + (rd / ld)^2 / 3*(1 + rd / ld) )
    def electrostatic(self, particle: list) -> list:
        pass

    # calculate thermophoresis force FT = -32/15 * (pi * m_n / 8 / kb / Tn)**0.5 * r_d^2 * k_tr * del_Tn
    def thermophoresis(self, particle: list) -> list:
        pass

    # calculate neutral drag force Fn = - 4/3 * pi * r_d^2 * n_n * m_n * v_thn * (v_d - v_n) 
    def neutral_drag(self, particle: list) -> list:
        pass

    # calculate ion drag force Fi = F_i-coll + F_i-or = ni * mi * vi * vs * pi * bc^2 +
    #                                                   ni * mi * vi * vs * 4 * pi * bn_pi/2^2 * gamma
    def ion_drag(self, particle: list) -> list:
        pass