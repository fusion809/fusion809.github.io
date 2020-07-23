# This file was generated, do not modify it. # hide
using FunctionIntegrator
a = abs(rombergs_method(x -> cos(x), 10, 0, pi/2)-1)
show(a)