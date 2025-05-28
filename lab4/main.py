from Core.defines import filePathODE, filePathPDE, variant, N, R, T0, Tw, p, tau, I_max, t_max, S
from ODE.absorption import initAbsorptionCoefTable as initAbsorpCoefODE
from PDE.absorption import initAbsorptionCoefTable as initAbsorpCoefPDE
from PDE.solution import solutePDE

import matplotlib.pyplot as plt


def main():
    T_tableODE, K_table = initAbsorpCoefODE(filePathODE, variant)
    T_tablePDE, sigma_table, lambda_table, cT_table = initAbsorpCoefPDE(filePathPDE)

    r, Thistory = solutePDE(N, R, T0, Tw, p, tau, I_max, t_max, S, T_tableODE, K_table, T_tablePDE, sigma_table, lambda_table, cT_table)

    plt.figure(figsize=(10, 6))
    for i, T in enumerate(Thistory):
        if i % 10 == 0 or i == len(Thistory) - 1:
            plt.plot(r, T, label=f"t = {i * tau:.2e} с")
    plt.xlabel("r, см")
    plt.ylabel("T(r), К")
    plt.title("Температурное поле во времени")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()