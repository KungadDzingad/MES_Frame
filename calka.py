import math


def f(x, y):
    return -2*x*x*y + 2*x*y + 4


def f2(x,y):
    return -5*math.pow(x, 2) *y


def calculate(func, a):
    if a == 2:
        w = [1, 1]
        p = [-1/math.sqrt(3), 1/math.sqrt(3)]
    elif a == 3:
        w = [5 / 9, 8 / 9, 5 / 9]
        p = [-math.sqrt(3/5), 0, math.sqrt(3/5)]
    else:
        raise Exception("Wrong number of points")

    suma = 0
    for i in range(len(w)):
        for j in range(len(w)):
            suma += w[i] * w[j] * func(p[i], p[j])

    return suma


print(calculate(f,2))





