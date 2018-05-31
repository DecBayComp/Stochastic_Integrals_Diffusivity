"""
Small functions used throughout the package
"""


def p(n, n_pi, dim):
    p = dim * (n + n_pi - 1.0) / 2.0 - 1.0
    return p


def n_pi(dim):
    return 5 - dim
