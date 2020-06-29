# This file was generated, do not modify it. # hide
using PyPlot
x = LinRange(-10,10,1001);
clf()
y = (x.^2).*exp.(-x.^2);
PyPlot.plot(x,y)
PyPlot.savefig(joinpath(@OUTPUT, "hermite_plot.png"), dpi=50)