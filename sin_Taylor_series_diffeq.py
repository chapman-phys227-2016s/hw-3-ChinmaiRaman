#! /usr/bin/env python

"""
File: sin_Taylor_series_diffeq.py
Copyright (c) 2016 Chinmai Raman
License: MIT
Course: PHYS227
Assignment: A.14
Date: Feb 25, 2016
Email: raman105@mail.chapman.edu
Name: Chinmai Raman
Description: Implements difference equations for computing a Taylor polynomial approximation to sin(x)
"""

import numpy as np

def sin_Taylor(x, N):
    n = 1
    an_prev = x
    sn_prev = 0
    while n <= N + 1:
        sn = sn_prev + an_prev
        an = - x**2 / ((2 * n + 1) * 2 * n) * an_prev
        sn_prev = sn
        an_prev = an
        n += 1
    return sn, abs(an)

def test_Taylor():
    assert(abs(sin_Taylor(np.pi / 2.0, 20)[0] - 1.0) < 1e-6), 'Failure'

def test_partc():
    assert(abs(sin_Taylor(np.pi / 2.0, 2)[0] - 1.0) < 1e-2), 'Failure'