# This file was generated, do not modify it. # hide
using FunctionIntegrator
a = abs(simpsons_rule(x -> cos(x), 1000, 0, pi/2)-1)
show(a)