# This file was generated, do not modify it. # hide
using FunctionIntegrator
a = abs(laguerre_quadrature(x -> x*exp(-x), 100, 1)-1)
show(a)