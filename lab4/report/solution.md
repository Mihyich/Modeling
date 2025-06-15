# Разностная схема

$
\left(\frac{\tau}{2 c_T r_n h^2}r_{n-\frac{1}{2}}\lambda_{n-\frac{1}{2}}^{s+1}\right)T_{n-1}^{s+1}+\left(1+\frac{\tau}{2 c_T r_n h^2} \left[r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s+1}+r_{n-\frac{1}{2}} \lambda_{n-\frac{1}{2}}^{s+1}\right]\right)T_n^{s+1}+\left(-\frac{\tau}{2 c_T r_n h^2}r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s+1}\right)T_{n+1}^{s+1}=\left(\frac{\tau}{2 c_T r_n h^2}r_{n-\frac{1}{2}}\lambda_{n-\frac{1}{2}}^{s}\right)T_{n-1}^s+\left(1-\frac{\tau}{2 c_T r_n h^2} \left[r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s}+r_{n-\frac{1}{2}} \lambda_{n-\frac{1}{2}}^{s}\right]\right)T_n^s+\left(\frac{\tau}{2 c_T r_n h^2}r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s}\right)T_{n+1}^s+\frac{r}{c_T}\left[\sigma(T_n^s)E^2-q_n^s\right]
$

## Коэффициенты

$$
A_n^{s+1}=\frac{\tau}{2 c_T r_n h^2}r_{n-\frac{1}{2}}\lambda_{n-\frac{1}{2}}^{s+1}
$$

$$
B_n^{s+1}=1+\frac{\tau}{2 c_T r_n h^2} \left[r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s+1}+r_{n-\frac{1}{2}} \lambda_{n-\frac{1}{2}}^{s+1}\right]
$$

$$
C_n^{s+1}=-\frac{\tau}{2 c_T r_n h^2}r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s+1}
$$

$$
F_n^s=\left(\frac{\tau}{2 c_T r_n h^2}r_{n-\frac{1}{2}}\lambda_{n-\frac{1}{2}}^{s}\right)T_{n-1}^s+\left(1-\frac{\tau}{2 c_T r_n h^2} \left[r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s}+r_{n-\frac{1}{2}} \lambda_{n-\frac{1}{2}}^{s}\right]\right)T_n^s+\left(\frac{\tau}{2 c_T r_n h^2}r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s}\right)T_{n+1}^s+\frac{r}{c_T}\left[\sigma(T_n^s)E^2-q_n^s\right]
$$

## Вспомогательные вычисления

$$
r_n=nh, \space n=0,\dots,N
$$

$$
r_{n+\frac{1}{2}}=\frac{r_n+r_{n+1}}{2}=\left(n+\frac{1}{2}\right)h
$$

$$
r_{n-\frac{1}{2}}=\frac{r_{n-1}-r_n}{2}=\left(n-\frac{1}{2}\right)h
$$

$$
\lambda_{n+\frac{1}{2}}^s=\frac{\lambda(T_n^s)+\lambda(T_{n+1}^s)}{2}
$$

$$
\lambda_{n-\frac{1}{2}}^s=\frac{\lambda(T_{n-1}^s)+\lambda(T_n^s)}{2}
$$

$$
\lambda_{n+\frac{1}{2}}^{s+1}=\frac{\lambda(T_n^{s+1})+\lambda(T_{n+1}^{s+1})}{2}
$$

$$
\lambda_{n-\frac{1}{2}}^{s+1}=\frac{\lambda(T_{n-1}^{s+1})+\lambda(T_n^{s+1})}{2}
$$

$$
E = \frac{1}{2 \pi \int_0^R \sigma(T(r))rdr}
$$

$$
q(r) = ck(r)(u_p - u)
$$

$$
u_p(r) = \frac{3.084 \cdot 10^{-4}}{\exp{\left( \frac{4.799 \cdot 10^4}{T(r)} \right)} - 1}
$$

$$
T(r) = (T_w - T_0) \left( \frac{r}{R} \right)^p + T_0
$$

> $$
> \sigma(T), \space \lambda(T), \space c(T)
> $$
>
> вычисляются с помощью логарифмической интерполяции по таблице

> Поиск $u(r)$ - решенная задача (лабораторная №2)
>
> Решение зависит только от переданной температуры $T$
>
> То есть на каждом слое $s$, необходимо вызвать решение для $u(r)$ при новой $Т$

# Переход между слоями

Общий вид полученной разностной схемы с найденными коэффициентами:

$$
A_n^{s+1}T_{n-1}^{s+1}+B_n^{}T_n^{s+1}+C_n^{s+1}T_{n+1}^{s+1}=F_n^s
$$

Рассмотрим $s$-й слой:

$$
A_n^{s+1,(k)}T_{n-1}^{s+1,(k)}+B_n^{s+1}T_n^{s+1,(k)}+C_n^{s+1,(k)}T_{n+1}^{s+1,(k)}=F_n^s
$$

Для нового слоя приближенно принимается:

$$
T_n^{s+1,(0)}=T_n^s
$$

Пересчитываются коэффициенты на основе $T_n^{s+1,(0)}$:

$$
\lambda_{n+\frac{1}{2}}^{s+1,(0)}=\frac{\lambda(T_n^{s+1,(0)})+\lambda(T_{n+1}^{s+1,(0)})}{2}
$$

$$
\lambda_{n-\frac{1}{2}}^{s+1,(0)}=\frac{\lambda(T_{n-1}^{s+1,(0)})+\lambda(T_n^{s+1,(0)})}{2}
$$

Вычисляются коэффициенты:

$$
A_n^{s+1,(0)}=\frac{\tau}{2 c_T r_n h^2}r_{n-\frac{1}{2}}\lambda_{n-\frac{1}{2}}^{s+1,(0)}
$$

$$
B_n^{s+1,(0)}=1+\frac{\tau}{2 c_T r_n h^2} \left[r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s+1,(0)}+r_{n-\frac{1}{2}} \lambda_{n-\frac{1}{2}}^{s+1,(0)}\right]
$$

$$
C_n^{s+1,(0)}=-\frac{\tau}{2 c_T r_n h^2}r_{n+\frac{1}{2}}\lambda_{n+\frac{1}{2}}^{s+1,(0)}
$$

Составляется матрица, $(k=0)$:

$$
A_n^{s+1,(k)}T_{n-1}^{s+1,(k)}+B_n^{s+1}T_n^{s+1,(k)}+C_n^{s+1,(k)}T_{n+1}^{s+1,(k)}=F_n^s
$$

Она трехдиагональная, решается методом прогонки с учетом краевых условий:

$$
\begin{cases}
    t=0, T(x, 0) = T_0 + (T_w - T_0) \left( \frac{r}{R} \right)^p \\
    r=0, \frac{\partial T}{\partial t}=0, \frac{du}{dr}=0 \\
    r=R, T(R)=T_w, -\frac{1}{3k(R)} \frac{du}{dr} = 0.39 \cdot u(R)
\end{cases}
$$

Откуда находятся новые $T_n^{s+1,(k+1)}$, то есть произошел переход на итерацию $k=1$. Условием остановки итераций служит:

$$
\max_n\left|\frac{T_n^{s+1,(k+1)} - T_n^{s+1,(k)}}{T_n^{s+1,(k+1)}}\right| > \varepsilon
$$