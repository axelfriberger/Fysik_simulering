import matplotlib.pyplot as plt
import math
import numpy as np


class Rocket:
    def __init__(self, mass, rocket_thrust, air_resistance: bool, 
                drag_coefficient, area, change_in_mass_and_grav_acc: bool):
        self.mass_0 = mass
        self.thrust = rocket_thrust
        self.air_resistance: bool = air_resistance
        self.drag_coefficient = drag_coefficient
        self.area = area
        self.change_in_mass_and_grav_acc: bool = change_in_mass_and_grav_acc

        self.velocity = 0
        self.height: int = 0
        self.g = 9.82
        self.acc = (self.thrust - self.mass_0*self.g)/self.mass_0
        self.time = 0
        self.empty_rocket_weight = 140_000 


    def mass(self, time: float) -> float:
        K = 4989
        return max(self.mass_0 - K*time, self.empty_rocket_weight) \
            if self.change_in_mass_and_grav_acc else self.mass_0

    def grav_acc(self, height: float) -> float:
        G = 6.6743*10**(-11)
        M = 5.97*10**24
        R = 6371*10**3
        return (G*M)/((R+height)**2) \
            if self.change_in_mass_and_grav_acc else self.g 

    def force(self, time, height):
        return self.thrust - self.mass(time)*self.grav_acc(height) \
        if self.mass(time) >= self.rocket_weight else -self.rocket_weight * self.grav_acc(height)

    def acceleration(self, time, height) -> float:
        a = (self.thrust - self.mass(time)*self.grav_acc(height))/self.mass(time)
        if self.air_resistance:
            a -= self.drag_acc(self.velocity, self.time)
        return a

    def air_density(self, height):
        rho0 = 1.225
        H = 8000.0 
        return rho0 * np.exp(-height/H)\
            if height < 100_000 else 0

    def drag_acc(self, velocity, time):
        return (self.drag_coefficient*self.air_density(self.height)*self.area*velocity**2)/(2*self.mass(time))

    def update(self, dt):
        self.height += self.velocity * dt + 0.5 * self.acceleration(self.time, self.height) * dt**2
        self.velocity += self.acceleration(self.time, self.height) * dt
        self.time += dt

    

class App:
    def __init__(self, *rockets):
        self.rockets = rockets
        self.time_values = [[] for _ in rockets]
        self.height_values = [[] for _ in rockets]

    def run(self, height_max, dt):
        for i, rocket in enumerate(self.rockets):
            while rocket.height <= height_max+50_000:
                rocket.update(dt)
                self.time_values[i].append(rocket.time)
                self.height_values[i].append(rocket.height)

    def plot(self):
        plt.figure()
        for i, rocket in enumerate(self.rockets):
            plt.plot(self.time_values[i], self.height_values[i], 
            label=f"Air Resistance={rocket.air_resistance}, Change in mass and grav acc={rocket.change_in_mass_and_grav_acc}")
            idx = next(x[0] for x in enumerate(self.height_values[i]) if x[1] >= height_max)
            plt.plot(self.time_values[i][idx], self.height_values[i][idx], 'bo') 
            plt.vlines(self.time_values[i][idx], 0, self.height_values[i][idx], colors='r', linestyles='--', color="black")

        plt.plot([x for x in range(int(self.time_values[0][-1]))], [height_max for y in range(int(self.time_values[0][-1]))], 
                    label="Hight = 408_000m",linestyle="--", color="black")

        plt.title("Rocket Height over Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Height (m)")
        #plt.legend()
        plt.show()


if __name__ == "__main__":
    mass = 3_000_000
    thrust = 34_500_000
    area = (5**2)*3.14

    bool_combinations = [(False, False), (True, False), (False, True), (True, True)]
    rockets = [Rocket(mass=mass, rocket_thrust=thrust, air_resistance=air_resistance,
                change_in_mass_and_grav_acc=change_in_mass_and_grav_acc, drag_coefficient=0.95, area=area) 
                for air_resistance, change_in_mass_and_grav_acc in bool_combinations]
    dt = 0.1
    height_max = 408_000*3
    app = App(*rockets)
    app.run(height_max, dt)
    app.plot()
    


