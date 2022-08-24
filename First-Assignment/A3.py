import numpy as np
import copy
import math

def LU(A):
    n = len(A)
    U = copy.deepcopy(A)
    P = [[0.0 for c in range(n)] for i in range(n)]
    L = [[0.0 for c in range(n)] for i in range(n)]
    for i in range (n):
        P[i][i] = 1.0
        L[i][i] = 1.0
        for j in range (n):
            U[i][j] *= 1.0
    max = 0
    so = 0
    for c in range (n-1):
        for i in range(c, n - 1):
            if abs(U[i][c]) < abs(U[i + 1][c]):
                max = i + 1
        if max != 0:
            U[c], U[max] = U[max], U[c]
            P[c], P[max] = P[max], P[c]
            for i in range (c, 0, -1):
                L[c][i-1], L[max][i-i] = L[max][i-i], L[c][i-1]
        for k in range(c, n - 1):
            so = U[k + 1][c] / U[c][c]
            L[k + 1][c] = so
            for i in range(n):
                U[k + 1][i] -= so * U[c][i]
    return L, U, P

def gauss(A,b):
    n = len(A)
    L, U, P = LU(A)
    x =  np.matmul(P, b)
    for c in range (n-1):
        for k in range(c, n - 1):
            so = L[k + 1][c] / L[c][c]
            x[k + 1] -= so * x[c]
            for i in range(n):
                L[k + 1][i] -= so * L[c][i]
    for c in range(n):
        x[c] /= L[c][c]
    for c in range(n-1, 0, -1):
        for k in range(c, 0, -1):
            so = U[k - 1][c] / U[c][c]
            x[k - 1] -= so * x[c]
            for i in range(n):
                U[k - 1][i] -= so * U[c][i]
    for c in range(n):
        x[c] /= U[c][c]
    return x

def cholesky(A):
    n = len(A)
    L = [[0.0 for c in range(n)] for i in range(n)]
    for c in range(n):
        for i in range(c+1):
            sum = 0
            if(i == c):
                for k in range (i):
                    sum += pow(L[i][k], 2)
                L[i][i] = math.sqrt(A[i][i]-sum)
            else:
                for k in range (i):
                    sum += L[c][k] * L[i][k]
                L[c][i] = (A[c][i] - sum) / L[i][i]
    return L

def GaussSiedel(n, tol):
    A = [[0 for c in range(n)] for i in range(n)]
    b = [1 for c in range(n)]
    b[0] = 3
    b[n-1] = 3
    for c in range (n):
        for i in range (c, n):
            if(c == i):
                A[c][c] = 5
            if (i == c + 1):
                A[c][i] = -2
                A[i][c] = -2
    x = [0 for c in range(n)]
    bol = True
    while(bol):
        x1 = copy.deepcopy(x)
        for c in range (n):
            sum = b[c]
            for i in range (n):
                if (c != i):
                    sum -= x[i] * A[c][i]
            x[c] = sum / A[c][c]
            max = 0
        for c in range (n):
            if(abs(x[c]-x1[c]))>max:
                max = abs(x[c]-x1[c])
        if max < tol:
            bol = False
    return x

A = [[1, -2, 1],
     [2, 1, -3],
     [4, -7, 1]]
b = [0, 5, -1]
c = gauss(A, b)
print("1) Vector x:\n", np.array([c]).T)

A = [[4, -2, 2],
     [-2, 2, -4],
     [2, -4, 11]]
L = cholesky(A)
print("2) Matrix L:\n", np.array(L))

A = [[8,1,1],[1,8,1],[1,1,8]]
b = [10,10,10]
x = GaussSiedel(10, 0.00005)
print("3) Vector x:\n", np.array([x]).T)