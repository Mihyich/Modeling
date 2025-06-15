from Core.defines import filePathODE, filePathPDE, folderPathPlot, variant, N, R, T0, Tw, p, tau, I_max, t_max, S
from Core.plotter import Plotter
from ODE.absorption import initAbsorptionCoefTable as initAbsorpCoefODE
from PDE.absorption import initAbsorptionCoefTable as initAbsorpCoefPDE
from PDE.solution import solutePDE


def main():
    plotter = Plotter(folderPathPlot)
    T_tableODE, K_table = initAbsorpCoefODE(filePathODE, variant)
    T_tablePDE, sigma_table, lambda_table, cT_table = initAbsorpCoefPDE(filePathPDE)
    r, Thistory = solutePDE(N, R, T0, Tw, p, tau, I_max, t_max, S, T_tableODE, K_table, T_tablePDE, sigma_table, lambda_table, cT_table, plotter)

    plotter.plot_T_evolution(r, Thistory, tau)


if __name__ == "__main__":
    main()