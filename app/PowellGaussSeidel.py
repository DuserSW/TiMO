import sympy
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import copy


def create_ksi(N):
    ksi = []
    for i in range(N):
        tmp = []
        for j in range(N):
            if i == j:
                tmp.append(1)
            else:
                tmp.append(0)

        ksi.append(tmp)

    return ksi


def get_value_from_function(expr, expr_args):
    _x1, _x2, _x3, _x4, _x5 = sympy.symbols("x1 x2 x3 x4 x5")
    length = len(expr_args)
    if length == 1:
        return expr.subs({_x1: expr_args[0]}).n()
    elif length == 2:
        return expr.subs({_x1: expr_args[0], _x2: expr_args[1]}).n()
    elif length == 3:
        return expr.subs({_x1: expr_args[0], _x2: expr_args[1], _x3: expr_args[2]}).n()
    elif length == 4:
        return expr.subs({_x1: expr_args[0], _x2: expr_args[1], _x3: expr_args[2], _x4: expr_args[3]}).n()
    elif length == 5:
        return expr.subs({_x1: expr_args[0], _x2: expr_args[1], _x3: expr_args[2], _x4: expr_args[3], _x5: expr_args[4]}).n()


def create_point(ksi, x, val):
    point = []

    for i in range(len(x)):
        if ksi[i] == 1 or ksi[i] == -1:
            point.append(val)
        else:
            point.append(x[i])

    return point


def heaviside(expr, expr_args, theta_i):
    if get_value_from_function(expr, expr_args) + theta_i > 0:
        return 1
    elif get_value_from_function(expr, expr_args) + theta_i <= 0:
        return 0


def powell_penalty_func(expr, expr_args, sigma, theta, limitations):
    expr_val = get_value_from_function(expr, expr_args)

    penalty_val = 0
    for i in range(len(limitations)):
        heaviside_val = heaviside(limitations[i], expr_args, theta[i])

        if heaviside_val == 1:
            penalty_val += (sigma[i] * pow(get_value_from_function(limitations[i], expr_args) + theta[i], 2))

    return expr_val + penalty_val


def golden_split_search(a, b, epsilon, max_iteration, expr, x, ksi, sigma, theta, limitations):
    i = 0
    alpha = (sympy.sqrt(5) - 1).n() / 2
    c = b - alpha * (b - a)
    d = a + alpha * (b - a)

    while abs(b - a) > epsilon:
        point_in_c = create_point(ksi, x, c)
        func_value_of_c = powell_penalty_func(expr, point_in_c, sigma, theta, limitations)

        point_in_d = create_point(ksi, x, d)
        func_value_of_d = powell_penalty_func(expr, point_in_d, sigma, theta, limitations)

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


def gauss_siedel(expr, x, sigma, theta, limitations, ksi, e_0, epsilon_j, L):
    steps_x1 = []
    steps_x2 = []
    steps_x1.append(x[0])
    steps_x2.append(x[1])

    for i in range(L):

        x_old = []
        for j in range(len(x)):
            x_old.append(x[j])

        for j in range(len(x)):
            first_move_x = golden_split_search(x[j], x[j] + e_0, epsilon_j, 100, expr, x, ksi[j], sigma, theta, limitations)
            second_move_x = golden_split_search(x[j], x[j] - e_0, epsilon_j, 100, expr, x, ksi[j], sigma, theta, limitations)

            first_move_point = create_point(ksi[j], x, first_move_x)
            second_move_point = create_point(ksi[j], x, second_move_x)

            if powell_penalty_func(expr, first_move_point, sigma, theta, limitations) < powell_penalty_func(expr, second_move_point, sigma, theta, limitations):
                x[j] = first_move_x
            else:
                x[j] = second_move_x

            if j == 0:
                steps_x1.append(x[0])
                steps_x2.append(x_old[1])
            else:
                steps_x1.append(x[0])
                steps_x2.append(x[1])

        sum = 0
        for j in range(len(x)):
            sum += math.sqrt((x[j] - x_old[j]) ** 2)

        if sum <= epsilon_j:
            return steps_x1, steps_x2

        if abs(powell_penalty_func(expr, x, sigma, theta, limitations) - powell_penalty_func(expr, x_old, sigma, theta, limitations)) <= epsilon_j:
            return steps_x1, steps_x2

    return steps_x1, steps_x2


def powell_calc_new_c_and_point(f, g, x, e_0, epsilon_j, sigma, theta, L, ksi):
    steps_x1, steps_x2 = gauss_siedel(f, x, sigma, theta, g, ksi, e_0, epsilon_j, L)
    new_point = [steps_x1[len(steps_x1) - 1], steps_x2[len(steps_x2) - 1]]

    step3_boolean = []
    for i in range(len(g)):
        if get_value_from_function(g[i], new_point) + theta[i] > 0:
            step3_boolean.append(True)
        else:
            step3_boolean.append(False)

    step3_values = []
    for i in range(len(g)):
        if step3_boolean[i]:
            step3_values.append(abs(get_value_from_function(g[i], new_point)))

    c = max(step3_values)

    return c, new_point


def powell(f, g, x, c, c_min, e_0, epsilon_j, L, ksi):
    list_of_steps = [copy.copy(x)]
    list_of_c = [copy.copy(c)]
    theta = [0, 0]
    sigma = [1, 1]
    last_exec_step6 = False
    m1 = 0.25
    m2 = 10

    c_0 = c
    c, new_point = powell_calc_new_c_and_point(f, g, x, e_0, epsilon_j, sigma, theta, L, ksi)
    list_of_steps.append(copy.copy(new_point))
    list_of_c.append(copy.copy(c))
    print("POWELL c = {0}, point = {1},{2}".format(c, new_point[0], new_point[1]))

    # krok4
    while c_min < c:
        # krok 5
        if c < c_0:
            # krok8
            if last_exec_step6:
                # krok8 (i)
                for i in range(len(g)):
                    theta[i] = max(get_value_from_function(g[i], new_point) + theta[i], 0)

                # krok7
                x = new_point
                last_exec_step6 = False
            else:
                # krok8 (ii)
                if c <= m1 * c_0:
                    # krok8 (i)
                    for i in range(len(g)):
                        theta[i] = max(get_value_from_function(g[i], new_point) + theta[i], 0)

                    # krok7
                    x = new_point
                    last_exec_step6 = False
                else:
                    # krok6
                    last_exec_step6 = True

                    for i in range(len(g)):
                        if abs(get_value_from_function(g[i], new_point)) > m1 * c_0 and get_value_from_function(g[i], new_point) + theta[i] > 0:
                            sigma[i] = m2 * sigma[i]
                            theta[i] = theta[i] / m2
                            # krok 7
                            x = new_point
        # False z kroku 5
        else:
            c = c_0

            # krok6
            last_exec_step6 = True

            for i in range(len(g)):
                if abs(get_value_from_function(g[i], new_point)) > m1 * c_0 and get_value_from_function(g[i], new_point) + theta[i] > 0:
                    sigma[i] = m2 * sigma[i]
                    theta[i] = theta[i] / m2
                    # krok 7
                    x = new_point

        c_0 = c
        c, new_point = powell_calc_new_c_and_point(f, g, x, e_0, epsilon_j, sigma, theta, L, ksi)
        list_of_steps.append(copy.copy(new_point))
        list_of_c.append(copy.copy(c))
        print("POWELL c = {0}, point = {1},{2}".format(c, new_point[0], new_point[1]))

    return list_of_steps, list_of_c