import matplotlib.pyplot as plt
import Galerkin.galerkin as galerkin
import Galerkin.galerkinNumpy as galerkinNumpy
import Galerkin.dragrace as dragrace
from Galerkin.defines import koefCount, NexusCount, h

def main():
    Xaxis_g, Yaxis_g, rssq_g = galerkin.solute(koefCount, h)
    Xaxis_gnp, Yaxis_gnp, rssq_gnp = galerkinNumpy.solute(koefCount, h)
    Xaxis_dr, Yaxis_dr = dragrace.solute(NexusCount)

    print(f"Сумма квадратов невязок в МНК: {rssq_g}")
    print(f"Сумма квадратов невязок в Numpy: {rssq_gnp}")

    plt.plot(Xaxis_g, Yaxis_g, "r", linewidth=1, label=f'МНК для решения СЛАУ с кол-вом коэффициентов: {koefCount}')
    plt.plot(Xaxis_gnp, Yaxis_gnp, "g", linewidth=1, label=f'Numpy для решения СЛАУ с кол-вом коэффициентов: {koefCount}')
    plt.plot(Xaxis_dr, Yaxis_dr, "b", linewidth=1, label=f'Прогонка с кол-вом точек дискретизации сетки: {NexusCount}')
    plt.title(f"Решение u'' - 2xu' + 2u = x")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.grid(True)
    plt.legend()
    plt.gcf().canvas.manager.set_window_title("График сравнения решений дифференциального уравнения")
    plt.show()

if __name__ == "__main__":
    main()