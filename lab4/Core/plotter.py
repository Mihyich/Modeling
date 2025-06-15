import os
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt

class Plotter:
    folderUup = None
    folderEnergySource = None
    folderEnergyStrength = None

    def __init__(self, folderPath):
        self.folderUup = folderPath / Path("u_up")
        self.folderEnergySource = folderPath / Path("EnergySource")
        self.folderEnergyStrength = folderPath / Path("EnergyStrenght")

        if not folderPath.exists():
            folderPath.mkdir(parents=True, exist_ok=True)

        if not self.folderUup.exists():
            self.folderUup.mkdir(parents=True, exist_ok=True)

        if not self.folderEnergySource.exists():
            self.folderEnergySource.mkdir(parents=True, exist_ok=True)

        if not self.folderEnergyStrength.exists():
            self.folderEnergyStrength.mkdir(parents=True, exist_ok=True)


    def plot_T_evolution(self, r: list[float], Thistory: list[list[float]], tau: float):
        matplotlib.use('TkAgg')

        plt.figure(figsize=(10, 6))
        for i, T in enumerate(Thistory):
            if i % 8 == 0 or i == len(Thistory) - 1:
                plt.plot(r, T, label=f"t = {i * tau:.2e} с")
        plt.title("Температурное поле во времени")
        plt.xlabel("r (см)")
        plt.ylabel("T(r) (К)")
        plt.legend()
        plt.grid(True)
        plt.show()


    def plot_u_and_up(self, r, u, u_p, i: int):
        matplotlib.use('Agg')
        fname = self.folderUup / Path(f"u_up{i}.svg")

        plt.figure(figsize=(10, 5))
        plt.plot(r, u, label="u(r)", color='blue')
        # plt.plot(r, u_p, '--', label=r"$u_p(r)$", color='orange')
        plt.title("Электромагнитное поле и функция Планка")
        plt.xlabel("r (см)")
        plt.ylabel("Значение")
        plt.legend()
        plt.grid(True)
        plt.savefig(fname=fname, format='svg', bbox_inches='tight', dpi=600)
        plt.close()


    def plot_q(self, r, q, i: int):
        matplotlib.use('Agg')
        fname = self.folderEnergySource / Path(f"energySource{i}.svg")

        plt.figure(figsize=(10, 5))
        plt.plot(r, q, color='green')
        plt.title("Источник тепла q(r)")
        plt.xlabel("r (см)")
        plt.ylabel("q(r) (Вт/см³)")
        plt.grid(True)
        plt.savefig(fname=fname, format='svg', bbox_inches='tight', dpi=600)
        plt.close()


    def plot_E_over_time(self, times, E_values):
        matplotlib.use('Agg')
        fname = self.folderEnergyStrength / Path(f"energyStrenght.svg")

        plt.figure(figsize=(10, 5))
        plt.plot(times, E_values, color='purple')
        plt.title("Напряжённость электрического поля E(t)")
        plt.xlabel("t (с)")
        plt.ylabel("E (В/см)")
        plt.grid(True)
        plt.savefig(fname=fname, format='svg', bbox_inches='tight', dpi=600)
        plt.close()


    def plot_current(self, times, I_values):
        plt.figure(figsize=(10, 5))
        plt.plot(times, I_values, color='red')
        plt.title("Ток во времени")
        plt.xlabel("t (с)")
        plt.ylabel("I(t) (А)")
        plt.grid(True)
        plt.show()


    def plot_material_functions(self, T_table, sigma_table, lambda_table, c_table):
        plt.figure(figsize=(12, 8))

        plt.subplot(3, 1, 1)
        plt.plot(T_table, sigma_table, marker='o', linestyle='--')
        plt.title("Электропроводность σ(T)")
        plt.grid(True)

        plt.subplot(3, 1, 2)
        plt.plot(T_table, lambda_table, marker='o', linestyle='--', color='green')
        plt.title("Теплопроводность λ(T)")
        plt.grid(True)

        plt.subplot(3, 1, 3)
        plt.plot(T_table, c_table, marker='o', linestyle='--', color='magenta')
        plt.title("Теплоёмкость c_T(T)")
        plt.grid(True)

        plt.tight_layout()
        plt.show()


    def plot_dT_dr(self, r, T):
        dr = r[1] - r[0]
        dT_dr = [(T[i+1] - T[i]) / dr for i in range(len(T)-1)]
        plt.figure(figsize=(10, 5))
        plt.plot(r[:-1], dT_dr, color='brown')
        plt.title("Градиент температуры ∂T/∂r")
        plt.xlabel("r (см)")
        plt.ylabel("dT/dr")
        plt.grid(True)
        plt.show()


    def plot_max_temperature_over_time(self, Thistory, times):
        max_temps = [max(T) for T in Thistory]
        plt.figure(figsize=(10, 5))
        plt.plot(times, max_temps, color='darkred')
        plt.title("Максимальная температура во времени")
        plt.xlabel("t (с)")
        plt.ylabel("T_max (К)")
        plt.grid(True)
        plt.show()


    # def plot_2d_heatmap(self, r, Thistory, times):
    #     T_array = np.array(Thistory)
    #     plt.figure(figsize=(12, 6))
    #     plt.imshow(T_array, aspect='auto', cmap='hot', extent=[0, R, t_max, 0])
    #     plt.colorbar(label="Температура, К")
    #     plt.xlabel("r (см)")
    #     plt.ylabel("t (с)")
    #     plt.title("Температурное поле T(r,t)")
    #     plt.show()