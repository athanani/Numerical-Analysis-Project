import numpy as np


def Simpson(a, b, N):
    c = (b - a) / N
    x = np.zeros(N + 1)
    y = np.zeros(N + 1)
    for i in range(N + 1):
        x[i] = i*c
        y[i] = np.sin(x[i])
    integral = y[0] + y[N]
    sum1 = 0
    for i in range(1, int(N/2)):
        sum1 += y[2*i]
    integral += 2 * sum1
    sum1 = 0
    for i in range(1, int(N/2) + 1):
        sum1 += y[(2*i) - 1]
    integral += 4 * sum1
    integral = integral * ((b - a) / (3 * N))
    M = y[0]
    for i in range (1, N + 1):
        if (abs(y[i]) > abs(M)):
            M = abs(y[i])
    error = (((b - a)**5)/(180*N**4))*M
    return integral, error


def Trapezoid(a, b, N):
    c = (b - a) / N
    x = np.zeros(N + 1)
    y = np.zeros(N + 1)
    for i in range(N + 1):
        x[i] = i*c
        y[i] = np.sin(x[i])
    integral = y[0] + y[N]
    sum1 = 0
    for i in range(1, N):
        sum1 += y[i]
    integral += 2 * sum1
    integral = integral * ((b - a) / (2 * N))
    M = abs(y[0])
    for i in range(1, N + 1):
        if (abs(-y[i]) > M):
            M = abs(-y[i])
    error = (((b - a) ** 3) / (12 * N ** 2)) * M
    return integral, error


a = 0
b = np.pi/2
N = 10
x, error = Simpson(a, b, N)
x1, error1 = Trapezoid(a, b, N)
print("Integral value of sin(x) in [0, π/2] with Simpson method: ", x)
print("Theoretical error of Simpson method: ", "%.16f" % error)
print("Arithmetic error of Simpson method: ", "%.16f" % abs(1 - x))
print("\nIntegral value of sin(x) in [0, π/2] with Trapezoid method: ", x1)
print("Theoretical error of Trapezoid method: ", "%.16f" % error1)
print("Arithmetic error of Simpson method: ", "%.16f" % (1 - x1))