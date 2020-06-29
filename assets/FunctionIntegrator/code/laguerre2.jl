# This file was generated, do not modify it. # hide
using FunctionIntegrator
a = abs(laguerre_quadrature(x -> x, 100, 2)-1)
show(a)