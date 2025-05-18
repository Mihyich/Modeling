import matplotlib.pyplot as plt
from RadiationTransfer.absoption import initAbsorptionCoefTable
from RadiationTransfer.solution import soluteODE, soluteFluxDifferential, soluteFluxIntegral
from RadiationTransfer.defines import N, T0, Tw, R, p, filePath, variant


def main():
    T, k = initAbsorptionCoefTable(filePath, variant)
    r, kinterp, u, u_p = soluteODE(N, T0, Tw, R, p, T, k)
    Fdif = soluteFluxDifferential(r, kinterp, u)
    Fint = soluteFluxIntegral(r, kinterp, u, u_p)

    plt.subplot(2, 2, 1)
    plt.plot(r, u_p, "g-", linewidth=1, label=f'u_p(r): Функия Планка ')
    plt.plot(r, u, "b-", linewidth=1, label=f'u(r): Решение')
    plt.title(f"Решение ОДУ")
    plt.xlabel("r")
    plt.ylabel("u(r), u_p(r)")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(r, Fdif, "r-", linewidth=1, label=f'F(z): Поток')
    plt.xlabel("z")
    plt.ylabel("F(z)")
    plt.title(f"Определение потока численным дифференцированием 2-ого порядка")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.plot(r, Fdif, "r-", linewidth=1, label=f'F(z): Поток')
    plt.xlabel("z")
    plt.ylabel("F(z)")
    plt.title(f"Определение потока интегрированием")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(r, kinterp, "m-", linewidth=1, label=f'k(r)')
    plt.xlabel("r")
    plt.ylabel("k(r)")
    plt.title(f"Интерполяция коэффициентов поглощения")
    plt.grid(True)
    plt.legend()

    plt.gcf().canvas.manager.set_window_title("График решение ОДУ")
    plt.show()


if __name__ == "__main__":
    main()