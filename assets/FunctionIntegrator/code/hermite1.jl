# This file was generated, do not modify it. # hide
using FunctionIntegrator
a = abs(hermite_quadrature(x -> x^2*exp(-x^2), 100, 1)-sqrt(pi)/2)
show(a)