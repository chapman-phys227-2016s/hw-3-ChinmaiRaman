#! /usr/bin/env python

"""
File: root_finder_examples.py
Copyright (c) 2016 Chinmai Raman
License: MIT
Course: PHYS227
Assignment: A.11
Date: Feb 24, 2016
Email: raman105@mail.chapman.edu
Name: Chinmai Raman
Description: Tests different methods for finding roots of a nonlinear function.
"""

import numpy as np

def Newton(f, fprime, x0, eps = 1e-7, N = 100):
    """
    Uses Newton's method to estimate the root(s) of a function
    """
    n = 1
    x = np.zeros(N + 1)
    x[0] = x0
    while abs(f(x[n - 1])) > eps and n <= N:
        if abs(fprime(float(x[n - 1]))) < 1e-14:
            raise ValueError("Error. Diverges due to small value of denominator")
        x[n] = x[n - 1] - f(float(x[n - 1])) / fprime(float(x[n - 1]))
        n += 1
    return x[:n]

def bisect(f, a, b, eps = 1e-7):
    """
    Uses the bisection method to estimate the root(s) of a function
    """
    f_a = f(a)
    if f_a * f(b) > 0:
        return None, 0

    i = 0
    m_list = []
    while b - a > eps:
        i += 1
        m = (b + a) / float(2)
        f_m = f(m)
        if f_a * f_m <= 0:
            b = m
        else:
            a = m
            f_a = f_m
        m_list.append(m)
    return m_list

def secant(f, x0, x1, eps = 1e-7, N = 100):
    """
    Uses the secant method to estimate the root(s) of a function
    """
    n = 2
    x = np.zeros(N + 1)
    x[0] = x0
    x[1] = x1
    while abs(f(x[n - 1]) * (x[n - 1] - x[n - 2])) > eps and n <= N:
        if abs((f(float(x[n-1])) - f(float(x[n - 2])))) < 1e-14:
            raise ValueError("Error. Diverges due to small value of denominator")
        x[n] = x[n - 1] - (f(x[n - 1]) * (x[n - 1] - x[n - 2])) / (f(float(x[n-1])) - f(float(x[n - 2])))
        n += 1
    return x[:n]

def test_Newton():
    def f(x):
        return 1 - x**2

    def fprime(x):
        return -2 * x

    assert(abs(f(Newton(f, fprime, -1)[-1]) - 0) < 1e-6), 'Failure'

def test_bisect():
    def f(x):
        return 1 - x**2

    assert(abs(f(bisect(f, -2, 0)[0]) - 0) < 1e-6), 'Failure'

def test_secant():
    def f(x):
        return 1 - x**2

    assert(abs(f(secant(f, 5, 3)[-1]) - 0) < 1e-4), 'Failure'