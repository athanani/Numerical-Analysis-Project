import math
import random

def f(x):
    return 94*(math.cos(x)**3)-24*math.cos(x)+177*(math.sin(x)**2)-108*(math.sin(x)**4)-72*(math.cos(x)**3)*(math.sin(x)**2)-65

def df(x):
    return (216*math.cos(x)**2-432*math.cos(x))*math.sin(x)**3+(-144*math.cos(x)**4-282*math.cos(x)**2+354*math.cos(x)+24)*math.sin(x)

def d2f(x):
    return (432-432*math.cos(x))*math.sin(x)**4+(1224*math.cos(x)**3-1296*math.cos(x)**2+564*math.cos(x)-354)*math.sin(x)**2-144*math.cos(x)**5-282*math.cos(x)**3+354*math.cos(x)**2+24*math.cos(x)

def modifiedNewtonRaphson(x, tol, i):
    if df(x) == 0:
        print("Divide by zero!")
    else:
        while abs((f(x)/df(x)) + ((f(x)**2)*d2f(x))/(2*(df(x)**3))) >= tol:
            i+=1
            x = x - (f(x)/df(x)) - ((f(x)**2)*d2f(x))/(2*(df(x)**3))
        return round(x, 5), i

def modifiedBisection(a, b, tol, i):
    c = 0
    if f(a) * f(b) > 0:
        print("No root found!")
    else:
        while (b - a) / 2.0 > tol:
            i += 1
            c = random.uniform(a,b)
            if f(c) == 0:
                return round(c, 5), i
            elif f(a) * f(c) < 0:
                b = c
            else:
                a = c
        return round(c, 5), i

def modifiedSecant(a, b, c, tol, i):
    if f(a) * f(b) * f(c) > 0:
        print("No root found!")
    else:
        while abs(f(c)) > tol:
            i+=1
            r = f(c) / f(b)
            q = f(a) / f(b)
            s = f(c) / f(a)
            d = c - (r*(r-q)*(c-b)+(1-r)*s*(c-a))/((q-1)*(r-1)*(s-1))
            if f(d) == 0:
                return round(d, 5), i
            else:
                a = b
                b = c
                c = d
        return round(d, 5), i

def Bisection(a, b, tol, i):
    if f(a) * f(b) > 0:
        print("No root found!")
    else:
        while (b - a) / 2.0 > tol:
            i += 1
            c = (a + b) / 2.0
            if f(c) == 0:
                return round(c, 5), i
            elif f(a) * f(c) < 0:
                b = c
            else:
                a = c
        return round(c, 5), i

def NewtonRaphson(x, tol, i):
    if df(x) == 0:
        print("Divide by zero!")
    else:
        while abs(f(x)/df(x)) >= tol:
            i+=1
            x = x - f(x)/df(x)
        return round(x, 5), i

def Secant(a, b, tol, i):
    if f(a) * f(b) > 0:
        print("No root found!")
    else:
        while abs(f(b)) > tol:
            i+=1
            c = b - (f(b)*(b-a))/(f(b)-f(a))
            if f(c) == 0:
                return round(c, 5), i
            else:
                a = b
                b = c
        return round(c, 5), i

tol = 0.000005

roota1, ia1 = modifiedNewtonRaphson(0.85, tol, 0)
roota2, ia2 = modifiedNewtonRaphson(1, tol, 0)
roota3, ia3 = modifiedNewtonRaphson(2.5, tol, 0)

rootb1, ib1 = modifiedBisection(0, 0.9, tol, 0)
rootb2, ib2 = modifiedBisection(0.9, 1.5, tol, 0)
rootb3, ib3 = modifiedBisection(1.5, 2.5, tol, 0)

rootc1, ic1 = modifiedSecant(0.9, 0.8, 0.7, tol, 0)
rootc2, ic2 = modifiedSecant(0.9, 1.5, 2.1, tol, 0)
rootc3, ic3 = modifiedSecant(2.5, 2.6, 3, tol, 0)

rootca1, cia1 = NewtonRaphson(0.85, tol, 0)
rootca2, cia2 = NewtonRaphson(1, tol, 0)
rootca3, cia3 = NewtonRaphson(2.5, tol, 0)

rootcb1, cib1 = Bisection(0, 0.9, tol, 0)
rootcb2, cib2 = Bisection(0.9, 1.5, tol, 0)
rootcb3, cib3 = Bisection(1.5, 2.5, tol, 0)

rootcc1, cic1 = Secant(0.9, 0.8, tol, 0)
rootcc2, cic2 = Secant(1, 1.5, tol, 0)
rootcc3, cic3 = Secant(2.3, 2.9, tol, 0)

print("1) Roots of the modified Newton Raphson Method:\n", roota1, ", number of iterations needed:", ia1, "\n",
      roota2,   ", number of iterations needed:", ia2, "\n", roota3, ", number of iterations needed:", ia3)

print("Roots of the modified Bisection Method:\n", rootb1, ", number of iterations needed:", ib1, "\n",
      rootb2,   ", number of iterations needed:", ib2, "\n", rootb3, ", number of iterations needed:", ib3)

print("Roots of the modified Secant Method:\n", rootc1, ", number of iterations needed:", ic1, "\n",
      rootc2,   ", number of iterations needed:", ic2, "\n", rootc3, ", number of iterations needed:", ic3)

print("Roots of the classic Newton Raphson Method:\n", rootca1, ", number of iterations needed:", cia1, "\n",
      rootca2,   ", number of iterations needed:", cia2, "\n", rootca3, ", number of iterations needed:", cia3)

print("Roots of the classic Bisection Method:\n", rootb1, ", number of iterations needed:", ib1, "\n",
      rootcb2,   ", number of iterations needed:", cib2, "\n", rootcb3, ", number of iterations needed:", cib3)

print("Roots of the classic Secant Method:\n", rootc1, ", number of iterations needed:", ic1, "\n",
      rootcc2,   ", number of iterations needed:", cic2, "\n", rootcc3, ", number of iterations needed:", cic3)

it = [[0 for c in range(3)] for i in range(10)]
for c in range (10):
    rb1, i1 = modifiedBisection(0, 0.9, tol, 0)
    rb2, i2 = modifiedBisection(0.9, 1.5, tol, 0)
    rb3, i3 = modifiedBisection(1.5, 2.5, tol, 0)
    it[c][0] = i1
    it[c][1] = i2
    it[c][2] = i3
print("2) Number of iterations of the modified Bisection Method\nfor each of the roots after 10 runs of the method:\n",
      rootb1, "(1st column),", rootb2, "(2nd column),", rootb3, "(3rd column).")
for c in it:  
    print(c)