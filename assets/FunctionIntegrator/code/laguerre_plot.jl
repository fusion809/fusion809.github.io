# This file was generated, do not modify it. # hide
using PyPlot
x = LinRange(0,15,1001);
clf()
y = x.*exp.(-x);
PyPlot.plot(x,y)
PyPlot.savefig(joinpath(@OUTPUT, "laguerre_plot.png"), dpi=80)