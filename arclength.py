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
import sympy as sp

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

def derivative(f='sin(x)'):
    x = sp.Symbol('x')
    df = sp.lambdify([x], sp.diff(f))
    return df

def helper(g):
    x = sp.Symbol('x')
    dg = sp.lambdify([x], sp)

def arclength(f, a, b, n):
    df = derivative(f)
    g = lambda x: np.sqrt(1 + (df(x) ** 2))
    return integral(g, a, b, n)

def test_derivative():
    def f(x):
        return x**2
    print derivative(f, 2)
    assert(abs(derivative(f, 2) - 4.0) < 1e-6)

def test_integral():
    def f(x):
        return 2 * x

    assert(abs(integral(f, 0, 4)[1][-1] - 16.0) < 1e-6)

def test_arclength():
    assert(abs(arclength('x', 0, 2, 100)[1][-1] - 2.8284271247461) < 1e-7)