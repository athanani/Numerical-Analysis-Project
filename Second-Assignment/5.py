import numpy as np
import copy
import matplotlib.pyplot as plt

def LU(A):
    n = len(A)
    U = copy.deepcopy(A)
    P = [[0.0 for c in range(n)] for i in range(n)]
    L = [[0.0 for c in range(n)] for i in range(n)]
    for i in range(n):
        P[i][i] = 1.0
        L[i][i] = 1.0
        for j in range(n):
            U[i][j] *= 1.0
    for c in range(n - 1):
        max = c
        for i in range(c+1, n):
            if abs(U[max][c]) < abs(U[i][c]):
                max = i
        if max != c:
            for j in range(n):
                met = U[c][j]
                U[c][j] = U[max, j]
                U[max][j] = met
                met = P[c][j]
                P[c][j] = P[max][j]
                P[max][j] = met
            for i in range(c, 0, -1):
                met = L[c][i-1]
                L[c][i-1] = L[max][i-1]
                L[max][i-1] = met
        for k in range(c+1, n):
            so = U[k][c]
            L[k][c] = U[k][c] / U[c][c]
            for i in range(n):
                U[k][i] -= so * U[c][i] / U[c][c]
    return L, U, P


def gauss(A, b):
    n = len(A)
    L, U, P = LU(A)
    x = np.matmul(P, b)
    for c in range(n - 1):
        for k in range(c, n - 1):
            so = L[k + 1][c] / L[c][c]
            x[k + 1] -= so * x[c]
            for i in range(n):
                L[k + 1][i] -= so * L[c][i]
    for c in range(n):
        x[c] /= L[c][c]
    for c in range(n - 1, 0, -1):
        for k in range(c, 0, -1):
            so = U[k - 1][c] / U[c][c]
            x[k - 1] -= so * x[c]
            for i in range(n):
                U[k - 1][i] -= so * U[c][i]
    for c in range(n):
        x[c] /= U[c][c]
    return x


def order(x, y):
    n = len(x)
    for c in range(n - 1):
        max1 = c
        for i in range(c + 1, n):
            if x[max1] > x[i]:
                max1 = i
        if max1 != c:
            x[c], x[max1] = x[max1], x[c]
            y[c], y[max1] = y[max1], y[c]
    return np.array(x), np.array(y)


def coefficients(x, y):
    n = len(x)
    A = np.zeros((n, n))
    dd = np.zeros(n)
    for c in range(n):
        A[c, 0] = y[c]
    for c in range(1, n):
        for j in range(n - c):
            A[j, c] = (A[j + 1, c - 1] - A[j, c - 1]) / (x[c + j] - x[j])
    for c in range(0, n):
        dd[c] = A[0, c]
    return np.array(dd)


def NewtonInterpolation(x, y, x0):
    n = len(x)
    dd = coefficients(x, y)
    sum1 = dd[0]
    for c in range(1, n):
        product = dd[c]
        for i in range(c):
            product *= (x0 - x[i])
        sum1 += product
    return sum1


def Splines(x, y, x0):
    n = len(x)
    index = 0
    delta = np.zeros(n)
    Delta = np.zeros(n)
    A = np.zeros((n, n))
    u = np.zeros(n)
    d = np.zeros(n)
    b = np.zeros(n)
    A[0, 0] = 1
    A[n - 1, n - 1] = 1
    a = copy.deepcopy(y)
    for c in range(n):
        if x[c] <= x0:
            index = c
    for c in range(n - 1):
        delta[c] = x[c + 1] - x[c]
        Delta[c] = y[c + 1] - y[c]
    for c in range(1, n - 1):
        A[c, c - 1] = delta[c - 1]
        A[c, c] = 2 * delta[c - 1] + 2 * delta[c]
        A[c, c + 1] = delta[c]
        u[c] = 3 * (Delta[c] / delta[c] - Delta[c - 1] / delta[c - 1])
    l = gauss(A, u)
    for c in range(n - 1):
        d[c] = (l[c + 1] - l[c]) / 3 * delta[c]
        b[c] = Delta[c] / delta[c] - (delta[c] / 3) * (2 * l[c] + l[c + 1])
    return a[index] + b[index] * (x0 - x[index]) + l[index] * pow(x0 - x[index], 2) + d[index] * pow(x0 - x[index], 3)


def LeastSquares(x, y, x0, degree):
    n = len(x)
    A = np.zeros((n, degree + 1))
    b = copy.deepcopy(y)
    for c in range(n):
        for i in range(degree + 1):
            A[c,i] = pow(x[c], i)
    ATA = np.matmul(A.T, A)
    ATb = np.matmul(A.T, b)
    xnew = gauss(ATA, ATb)
    sum1 = 0
    for i in range(degree + 1):
        sum1 += xnew[i] * pow(x0, i)
    return sum1


def plot(x, y):
    figure = plt.figure(figsize=(7, 4))
    x1 = np.linspace(-np.pi, np.pi, 200)
    y1 = np.sin(x1)
    err1 = np.zeros(len(x1))
    err2 = np.zeros(len(x1))
    err3 = np.zeros(len(x1))
    y2 = [NewtonInterpolation(x, y, x1[i]) for i in range(len(x1))]
    plt.title("Newton Polynomial Interpolation Error")
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.grid(True)
    for i in range (len(x1)):
        err1[i] = y1[i] - y2[i]
    plt.plot(x1, err1, markersize=1, color='r')
    plt.savefig('NPIE.png', dpi=300)
    plt.show()
    y2 = [Splines(x, y, x1[i]) for i in range(len(x1))]
    plt.title("Cubic Spline Interpolation Error")
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.grid(True)
    for i in range (len(x1)):
        err2[i] = y1[i] - y2[i]
    plt.plot(x1, err2, markersize=1, color='g')
    plt.savefig('CSIE.png', dpi=300)
    plt.show()
    y2 = [LeastSquares(x, y, x1[i], 4) for i in range(len(x1))]
    plt.title("Least Squares Interpolation Error")
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.grid(True)
    for i in range (len(x1)):
        err3[i] = y1[i] - y2[i]
    plt.plot(x1, err3, markersize=1, color='b')
    plt.savefig('LSIE.png', dpi=300)
    plt.show()
    return err1, err2, err3


x = np.array([-3.14, -1.16, 3.14, 2.59, 0, -0.67, -1.69, 0.95, -2.56, 1.75])
y = np.sin(x)
for i in range (len(y)):
    y[i] = round(y[i], 8)
x, y = order(x, y)
y1 = np.zeros(len(x))
y2 = np.zeros(len(x))
y3 = np.zeros(len(x))
y4 = np.zeros(len(x))
for i in range(len(x)):
    y1[i] = NewtonInterpolation(x, y, x[i])
    y2[i] = Splines(x, y, x[i])
    y3[i] = LeastSquares(x, y, x[i], 4)

sum1 = 0
print("Value approach with Newton Polynomial Interpolation:")
for i in range(len(y1)):
    print(y1[i], "for x =", x[i])
    sum1 += y[i] - y1[i]
print("Average error :", "%.19f" % (sum1/len(y)))
sum1 = 0
print("\nValue approach with Cubic Spline Interpolation:")
for i in range(len(y2)):
    print(y2[i], "for x =", x[i])
    sum1 += y[i] - y2[i]
print("Average error :", (sum1/len(y)))
sum1 = 0
print("\nValue approach with Least Squares Interpolation:")
for i in range(len(y3)):
    print(y3[i], "for x =", x[i])
    sum1 += y[i] - y3[i]
print("Average error :", "%.19f" % (sum1/len(y)))


errorN, errorS, errorLS = plot(x, y)
sum1 = 0
acc = 0
for i in range (len(errorN)):
    sum1 += errorN[i]
sum1 /= len(errorN)
for i in range (19):
    if (0.5 * (10 ** (-i)) < abs(sum1)):
        break
    acc = i
print("\nAverage error for Newton Polynomial Interpolation:", "%.10f" % sum1, "and ", acc, "digits of accuracy.")
sum1 = 0
acc = 0
for i in range (len(errorS)):
    sum1 += errorS[i]
sum1 /= len(errorS)
for i in range (19):
    if (0.5 * (10 ** (-i)) < abs(sum1)):
        break
    acc = i
print("Average error for Cubic Splines Interpolation:", "%.10f" % sum1, "and ", acc, "digits of accuracy.")
sum1 = 0
acc = 0
for i in range (len(errorLS)):
    sum1 += errorLS[i]
sum1 /= len(errorLS)
for i in range (19):
    if (0.5 * (10 ** (-i)) < abs(sum1)):
        break
    acc = i
print("Average error for Least Squares Interpolation:", "%.10f" % sum1, "and ", acc, "digits of accuracy.")