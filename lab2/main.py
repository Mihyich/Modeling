import numpy as np
import matplotlib.pyplot as plt

import main2

c = 3e10
R = 0.35
Tw = 2000
T0 = 1e4
p = 4
last_up = 0
mas = []
def T(r):
    return (Tw - T0) * (r / R)**p + T0

def up(r):
    return (3.084e-4) / (np.exp(4.799e4 / T(r)) - 1)

def k(T):
    
    T_values = [2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
    k_values = [1.600e0, 5.400e0, 1.280e1, 2.500e1, 4.320e1, 6.860e1, 1.024e2, 1.458e2, 2.000e2]
    # k_values = [8.200e-3, 2.768e-2, 6.560e-2, 1.281e-1, 2.214e-1, 3.516e-1, 5.248e-1, 7.472e-1, 1.025e0]
    ln_T_values = np.log(T_values)
    ln_k_values = np.log(k_values)
    
    ln_T = np.log(T)
    ln_k = np.interp(ln_T, ln_T_values, ln_k_values)
    
    return np.exp(ln_k)

def system(r, y):
    u, F = y
    k_r = k(T(r))
    up_r = up(r)
    dudr = -3 * k_r / c * F
    
    if r == 0:
        dFdr = -c * k_r * (u - up_r) / 2
    else:
        dFdr = -F / r - c * k_r * (u - up_r)
    
    return np.array([dudr, dFdr])

def runge_kutta_4(f, r, y, h):
    k1 = f(r, y)
    k2 = f(r + h/2, y + h/2 * k1)
    k3 = f(r + h/2, y + h/2 * k2)
    k4 = f(r + h, y + h * k3)
    return y + h/6 * (k1 + 2*k2 + 2*k3 + k4)

def solve_ode(f, r0, y0, r_end, h):
    global mas
    r_values = np.arange(r0, r_end + h, h)
    y_values = []
    y = y0
    for r in r_values:
        y_values.append(y)
        y_new = runge_kutta_4(f, r, y, h)
        
        u_new = y_new[0]
        while True:
            if abs((up(r) - u_new) / up(r)) < 10**-6:
                h /= 2
                up_r = up(r)
                y_new[0] = up_r
                mas.append(h)
            else:
                h = 0.0001
                break
        
        y = y_new
    
    return r_values, np.array(y_values)

def compute_error(lambda_):
    r0 = 0
    u0 = lambda_ * up(r0)
    F0 = 0
    y0 = np.array([u0, F0])
    
    r_end = R
    h = 0.001
    r_values, y_values = solve_ode(system, r0, y0, r_end, h)
    
    u_R = y_values[-1, 0]
    F_R = y_values[-1, 1]
    dudr_R = -3 * k(T(R)) / c * F_R
    
    left_side = -1 / (3 * k(T(R))) * dudr_R
    right_side = 0.39 * u_R
    
    return abs(left_side - right_side)

def golden_section_search(left, right, precision):
    phi = (1 + np.sqrt(5)) / 2 
    while abs(right - left) > precision:
        mid1 = right - (right - left) / phi
        mid2 = left + (right - left) / phi
        
        error1 = compute_error(mid1)
        error2 = compute_error(mid2)
        
        if error1 < error2:
            right = mid2
        else:
            left = mid1
    
    return (left + right) / 2

lambda_opt = golden_section_search(0, 1, 1e-6)
print(f"Оптимальное значение кси: {lambda_opt}")

r0 = 0
u0 = lambda_opt * up(r0)
F0 = 0
y0 = np.array([u0, F0])
r_end = R
h = 0.01
r_values, y_values = solve_ode(system, r0, y0, r_end, h)

u_values = y_values[:, 0]
F_values = y_values[:, 1]
up_values = [up(r) for r in r_values]
k_values = [k(T(r)) for r in r_values]

plt.figure(figsize=(12, 8))

# Первый подграфик: Поток излучения F(r)
plt.subplot(2, 2, 1)
plt.plot(r_values, F_values, label=f'F(r)')
plt.xlabel('r')
plt.ylabel('F(r)')
plt.title('Поток излучения F(r)')
plt.legend()

r_vals, u_vals = main2.trapezoidal_method(R, main2.dr, main2.u0_opt)

# Второй подграфик: Объемная плотность энергии u(r) и функция Планка up(r)
plt.subplot(2, 2, 2)
plt.plot(r_values, u_values, label=f'u(r)', color='blue')
#plt.plot(r_values, up_values, label=f'up(r)', color='red', linestyle='--')
#plt.plot(r_values, [v for v in u_vals[:, 0]], color="black")
plt.xlabel('r')
plt.ylabel('Значения')
plt.title('Объемная плотность энергии u(r) и функция Планка up(r)')
plt.legend()

# Третий подграфик: Коэффициент поглощения k(r)
plt.subplot(2, 2, 3)
plt.plot(r_values, k_values, label=f'k(r)')
plt.xlabel('r')
plt.ylabel('k(r)')
plt.title('Коэффициент поглощения k(r)')
plt.legend()

plt.tight_layout()
plt.show()