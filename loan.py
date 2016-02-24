#! /usr/bin/env python

"""
File: loan.py
Copyright (c) 2016 Chinmai Raman
License: MIT
Course: PHYS227
Assignment: A.4
Date: Feb 23, 2016
Email: raman105@mail.chapman.edu
Name: Chinmai Raman
Description: Computes the development of a loan over time
"""

import numpy as np

def loan(p, L, N):
    """
    Calculates the amount paid each month and the amount of the loan remaining at each month.
    """
    index_set = range(N + 1)
    x = np.zeros(len(index_set))
    x[0] = L
    y = np.zeros_like(x)
    for n in index_set[1:]:
        y[n] = p / (12.0 * 100.0) * x[n - 1] + L / float(N)
        x[n] = x[n - 1] + p / (12.0 * 100.0) * x[n - 1] - y[n]
    return y, x

def test_loan():
    assert(abs(loan(6, 10000, 12)[1][12] - 0.0) < 1e-10)