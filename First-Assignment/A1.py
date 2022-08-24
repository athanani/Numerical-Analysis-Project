import math

def f(x):
    return math.exp(math.sin(x)**3) + x**6 - 2*x**4 - x**3 - 1

def df(x):
    return 3*math.exp(math.sin(x)**3)*math.cos(x)*math.sin(x)**2 + 6*x**5 - 8*x**3 - 3*x**2

def conv(x, tol):
    if df(x) == 0:
        return False
    else:
        if abs(df(x)) < tol:
            return False
        else:
            return True

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
roota1, ia1 = Bisection(-1.5, 1.5, tol, 0)
roota2, ia2 = Bisection(1.5, 2, tol, 0)
roota3, ia3 = Bisection(-2, -1, tol, 0)

rootb1, ib1 = NewtonRaphson(1, tol, 0)
rootb2, ib2 = NewtonRaphson(1.6, tol, 0)
rootb3, ib3 = NewtonRaphson(-1.5, tol, 0)

rootc1, ic1 = Secant(-1.2, 0.2, tol, 0)
rootc2, ic2 = Secant(1.3, 2, tol, 0)
rootc3, ic3 = Secant(-1.3, -1.15, tol, 0)

print("Roots of the Bisection Method:\n1)", roota1, ", number of iterations needed:", ia1, "\n2)",
      roota2,   ", number of iterations needed:", ia2, "\n3)", roota3, ", number of iterations needed:", ia3)

print("Roots of the Newton Raphson Method:\n1)","%.5f" %rootb1, ", number of iterations needed:", ib1, "\n2)",
      rootb2,   ", number of iterations needed:", ib2, "\n3)", rootb3, ", number of iterations needed:", ib3)

print("Roots of the Secant Method:\n1)", rootc1, ", number of iterations needed:", ic1, "\n2)",
      rootc2,   ", number of iterations needed:", ic2, "\n3)", rootc3, ", number of iterations needed:", ic3)

if(conv(rootb1, tol)):
    print("Root", "%.5f" %rootb1, "converges quadratically.")
else:
    print("Root", "%.5f" % rootb1, "converges linearly.")

if(conv(rootb2, tol)):
    print("Root", rootb2, "converges quadratically.")
else:
    print("Root", rootb2, "converges linearly.")

if(conv(rootb3, tol)):
    print("Root", rootb3, "converges quadratically.")
else:
    print("Root", rootb3, "converges linearly.")