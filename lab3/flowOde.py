import matplotlib.pyplot as plt
from RadiationTransfer.absoption import initAbsorptionCoefTable
from RadiationTransfer.solution import solute
from RadiationTransfer.defines import N, T0, Tw, R, p, filePath, variant

def main():
    T, k = initAbsorptionCoefTable(filePath, variant)
    r, u = solute(N, T0, Tw, R, p, T, k)

    plt.plot(r, u, "r", linewidth=1, label=f'u(r)')
    plt.title(f"Решение ОДУ")
    plt.xlabel("r")
    plt.ylabel("u(r)")
    plt.grid(True)
    plt.legend()
    plt.gcf().canvas.manager.set_window_title("График решение ОДУ")
    plt.show()


if __name__ == "__main__":
    main()