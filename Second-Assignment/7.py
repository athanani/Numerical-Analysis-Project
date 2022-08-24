import numpy as np
import copy


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
    return round(sum1, 2)


x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y = np.array([11.06, 10.64, 10.7, 10.64, 10.64, 11.02, 10.8, 10.66, 10.56, 10.4])
y1 = np.zeros(len(x))
y2 = np.zeros(4)
print("Closing price prediction of TITC for one session later (10/7/2020):")
xn = LeastSquares(x, y, 12, 2)
print(xn, ", with 2nd degree polynomial.")
xn = LeastSquares(x, y, 12, 3)
print(xn, ", with 3rd degree polynomial.")
xn = LeastSquares(x, y, 12, 4)
print(xn, ", with 4rd degree polynomial.\n")
for c in range (2, 5):
    for i in range (len(x)):
        y1[i]= LeastSquares(x, y, i, c)
    print("The first 10 closing prices of TITC with polynomial of degree", c, ":\n", y1, "\n")

print("Closing price prediction of TITC for five sessions later (16/7/2020):")
for c in range (2, 5):
    for i in range (13, 17):
        y2[i-13] = LeastSquares(x, y, i, c)
    print("Prediction of closing prices for TITC for 5 days later (16/7/2020) with polynomial of degree", c, ":\n", y2, "\n")

y = np.array([12.94, 12.88, 12.66, 12.64, 12.68, 12.8, 12.8, 12.88, 12.88, 12.86])
print("Closing price prediction of ΕΛΛ for one session later (10/7/2020):")
xn = LeastSquares(x, y, 12, 2)
print(xn, ", with 2nd degree polynomial.")
xn = LeastSquares(x, y, 12, 3)
print(xn, ", with 3rd degree polynomial.")
xn = LeastSquares(x, y, 12, 4)
print(xn, ", with 4rd degree polynomial.\n")
for c in range (2, 5):
    for i in range (len(x)):
        y1[i]= LeastSquares(x, y, i, c)
    print("The first 10 closing prices of ΕΛΛ with polynomial of degree", c, ":\n", y1, "\n")
for c in range (2, 5):
    for i in range (13, 17):
        y2[i-13] = LeastSquares(x, y, i, c)
    print("Prediction of closing prices for ΕΛΛ for 5 days later (16/7/2020) with polynomial of degree", c, ":\n", y2, "\n")