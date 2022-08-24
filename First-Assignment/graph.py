import matplotlib.pyplot as plt
import numpy as np
import math

def f(x):
    return math.exp(math.sin(x)**3) + x**6 - 2*x**4 - x**3 - 1

def graph():
    x = np.linspace(-2, 2, 100)
    y = np.e**((np.sin(x))**3) + x**6 - 2*x**4 - x**3 - 1
    figure = plt.figure(figsize=(7, 4))
    plt.plot(x, y)
    plt.xlim([-3, 3])
    plt.ylim([-3, 3])
    plt.xticks([-2, -1, 0, 1, 2])
    plt.yticks([-2, -1, 0, 1, 2])
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.legend(["f(x)"])
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(True)
    plt.savefig('graph.png', dpi=300)
    plt.show()

graph()