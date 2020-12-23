# This file was generated, do not modify it. # hide
using FunctionIntegrator
a = abs(chebyshev_quadrature(x -> cos(x), 1000, 1, 0, pi/2) - 1);
show(a)