#! /usr/bin/env python

"""
File: arclength.py
Copyright (c) 2016 Chinmai Raman
License: MIT
Course: PHYS227
Assignment: A.13
Date: Feb 25, 2016
Email: raman105@mail.chapman.edu
Name: Chinmai Raman
Description: Computes the arc length of a curve
"""

import numpy as np
from sympy import *

def integral(f, a, x, N = 20):
    index_set = range(N + 1)
    x = np.linspace(a, x, N + 1)
    f_ = np.zeros_like(x)
    s = np.zeros_like(x)
    f_[0] = f(x[0])
    s[0] = 0

    for n in index_set[1:]:
        f_[n] = f(x[n])
        s[n] = s[n - 1] + 0.5 * (x[n] - x[n - 1]) * (f_[n - 1] + f_[n])
    return x, s

def derivative(f):
    x = Symbol('x')
    f = f(x)
    return f.diff(x)

def f(x):
    f = x**2

def arclength(f, a, b, n):
    integral(np.sqrt(1 + ((derivative(f, a)) ** 2)), a, x)

def test_derivative():
    def f(x):
        return x**2
    print derivative(f, 2)
    assert(abs(derivative(f, 2) - 4.0) < 1e-6)

def test_integral():
    def f(x):
        return 2 * x

    assert(abs(integral(f, 0, 4)[1][-1] - 16.0) < 1e-6)