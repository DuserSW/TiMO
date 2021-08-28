import sympy
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

x1, x2 = sympy.symbols("x1 x2")
expr1 = x1 + 1 / x1 ** 2
expr2 = (x1 - 2) ** 2 + (x2 - 1) ** 2


def get_value_from_function(expr, expr_args):
    _x1, _x2, _x3, _x4, _x5 = sympy.symbols("x1 x2 x3 x4 x5")
    length = len(expr_args)
    if length == 1:
        return expr.subs({_x1: expr_args[0]}).n()
    elif length == 2:
        return expr.subs({_x1: expr_args[0], _x2: expr_args[1]}).n()


def create_point(ksi, x, val):
    point = []

    for i in range(len(x)):
        if ksi[i] == 1 or ksi[i] == -1:
            point.append(val)
        else:
            point.append(x[i])

    return point


def golden_split_search(a, b, epsilon, max_iteration, expr, x, ksi):
    i = 0
    alpha = (sympy.sqrt(5) - 1).n() / 2
    c = b - alpha * (b - a)
    d = a + alpha * (b - a)

    while abs(b - a) > epsilon:
        point_in_c = create_point(ksi, x, c)
        func_value_of_c = get_value_from_function(expr, point_in_c)

        point_in_d = create_point(ksi, x, d)
        func_value_of_d = get_value_from_function(expr, point_in_d)

        if func_value_of_c < func_value_of_d:
            a = a
            b = d
            d = c
            c = b - alpha * (b - a)
        else:
            a = c
            b = b
            c = d
            d = a + alpha * (b - a)

        i = i + 1

        if i > max_iteration:
            (a + b) / 2

    return (a + b) / 2


def gauss_siedel(expr, x, ksi, e_0, epsilon_j, L):
    steps_x1 = []
    steps_x2 = []
    steps_x1.append(x[0])
    steps_x2.append(x[1])

    for i in range(L):

        x_old = []
        for j in range(len(x)):
            x_old.append(x[j])

        for j in range(len(x)):
            first_move_x = golden_split_search(x[j], x[j] + e_0, epsilon_j, 100, expr, x, ksi[j])
            second_move_x = golden_split_search(x[j], x[j] - e_0, epsilon_j, 100, expr, x, ksi[j])

            first_move_point = create_point(ksi[j], x, first_move_x)
            second_move_point = create_point(ksi[j], x, second_move_x)

            if get_value_from_function(expr, first_move_point) < get_value_from_function(expr, second_move_point):
                x[j] = first_move_x
            else:
                x[j] = second_move_x

            if j == 0:
                steps_x1.append(x[0])
                steps_x2.append(x_old[1])
            else:
                steps_x1.append(x[0])
                steps_x2.append(x[1])
        print(steps_x1, steps_x2)
        sum = 0
        for j in range(len(x)):
            sum += math.sqrt((x[j] - x_old[j]) ** 2)

        if sum <= epsilon_j:
            print("first stop condition passed!")
            return steps_x1, steps_x2

        if abs(get_value_from_function(expr, x) - get_value_from_function(expr, x_old)) <= epsilon_j:
            print("second stop condition passed!")
            return steps_x1, steps_x2

    print("third stop condition passed!")
    return steps_x1, steps_x2
