# This file was generated, do not modify it. # hide
using FunctionIntegrator
a = abs(simpsons38_rule(x -> cos(x), 1002, 0, pi/2)-1)
show(a)