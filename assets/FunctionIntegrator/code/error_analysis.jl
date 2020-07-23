# This file was generated, do not modify it. # hide
using FunctionIntegrator, PyPlot
nc = round.(exp.(1/4*log(10)*(1:24)));
chebyshev1_error         = zeros(length(nc));
chebyshev2_error         = zeros(length(nc));
chebyshev3_error         = zeros(length(nc));
chebyshev4_error         = zeros(length(nc));
jacobi_error             = zeros(length(nc));
laguerre_error           = zeros(length(nc));
legendre_error           = zeros(length(nc));
lobatto_error            = zeros(length(nc));
radau_error              = zeros(length(nc));
rectangle_left_error     = zeros(length(nc));
rectangle_midpoint_error = zeros(length(nc));
rectangle_right_error    = zeros(length(nc));
rombergs_error           = zeros(30);
simpsons_error           = zeros(length(nc));
simpsons38_error         = zeros(length(nc));
trapezoidal_error        = zeros(length(nc));

for i=1:length(nc)
    chebyshev1_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 1, 0, 100)-pi/2);
    chebyshev2_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 2, 0, 100)-pi/2);
    chebyshev3_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 3, 0, 100)-pi/2);
    chebyshev4_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 4, 0, 100)-pi/2);
    jacobi_error[i]             = abs(jacobi_quadrature(x -> sech(x), nc[i], 1, 1, 0, 100)-pi/2);
    laguerre_error[i]           = abs(laguerre_quadrature(x -> sech(x), nc[i], 1)-pi/2);
    legendre_error[i]           = abs(legendre_quadrature(x -> sech(x), nc[i], 0, 100)-pi/2);
    lobatto_error[i]            = abs(lobatto_quadrature(x -> sech(x), nc[i], 0, 100)-pi/2);
    radau_error[i]              = abs(radau_quadrature(x -> sech(x), nc[i], 0, 100)-pi/2);
    rectangle_left_error[i]     = abs(rectangle_rule_left(x -> sech(x), nc[i], 0, 100)-pi/2);
    rectangle_midpoint_error[i] = abs(rectangle_rule_midpoint(x -> sech(x), nc[i], 0, 100)-pi/2);
    rectangle_right_error[i]    = abs(rectangle_rule_right(x -> sech(x), nc[i], 0, 100)-pi/2);
    simpsons_error[i]           = abs(simpsons_rule(x -> sech(x), 2*round(nc[i]/2), 0, 100)-pi/2);
    simpsons38_error[i]         = abs(simpsons38_rule(x -> sech(x), 3*round(nc[i]/3), 0, 100)-pi/2);
    trapezoidal_error[i]        = abs(trapezoidal_rule(x -> sech(x), nc[i], 0, 100)-pi/2);
end

for i=1:30
    rombergs_error[i]           = abs(rombergs_method(x -> sech(x), i, 0, 100)-pi/2);
end

PyPlot.figure(1)
PyPlot.clf()
PyPlot.plot(nc, chebyshev1_error);
PyPlot.title("Chebyshev k=1")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev1_error_plot.png"), dpi=80)
PyPlot.figure(2)
PyPlot.clf()
PyPlot.plot(nc, chebyshev2_error); 
PyPlot.title("Chebyshev k=2")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev2_error_plot.png"), dpi=80)
PyPlot.figure(3)
PyPlot.clf()
PyPlot.plot(nc, chebyshev3_error); 
PyPlot.title("Chebyshev k=3")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev3_error_plot.png"), dpi=80)
PyPlot.figure(4)
PyPlot.clf()
PyPlot.plot(nc, chebyshev4_error)
PyPlot.title("Chebyshev k=4")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev4_error_plot.png"), dpi=80)
PyPlot.figure(5)
PyPlot.clf()
PyPlot.plot(nc, jacobi_error)
PyPlot.title("Jacobi")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "jacobi_error_plot.png"), dpi=80)
PyPlot.figure(6)
PyPlot.clf()
PyPlot.plot(nc, laguerre_error)
PyPlot.title("Laguerre")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "laguerre_error_plot.png"), dpi=80)
PyPlot.figure(7)
PyPlot.clf()
PyPlot.plot(nc, legendre_error)
PyPlot.title("Legendre")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "legendre_error_plot.png"), dpi=80)
PyPlot.figure(8)
PyPlot.clf()
PyPlot.plot(nc, lobatto_error)
PyPlot.title("Lobatto")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "lobatto_error_plot.png"), dpi=80)
PyPlot.figure(9)
PyPlot.clf()
PyPlot.plot(nc, radau_error)
PyPlot.title("Radau")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "radau_error_plot.png"), dpi=80)
PyPlot.figure(10)
PyPlot.clf()
PyPlot.plot(nc, rectangle_left_error)
PyPlot.title("Rectangle rule (left)")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rectangle_left_error_plot.png"), dpi=80)
PyPlot.figure(11)
PyPlot.clf()
PyPlot.plot(nc, rectangle_midpoint_error)
PyPlot.title("Rectangle rule (midpoint)")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rectangle_midpoint_error_plot.png"), dpi=80)
PyPlot.figure(12)
PyPlot.clf()
PyPlot.plot(nc, rectangle_right_error)
PyPlot.title("Rectangle rule (right)")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rectangle_right_error_plot.png"), dpi=80)
PyPlot.figure(13)
PyPlot.clf()
PyPlot.plot(1:30, rombergs_error)
PyPlot.title("Romberg's method")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rombergs_error_plot.png"), dpi=80)
PyPlot.figure(14)
PyPlot.clf()
PyPlot.plot(2*round.(nc/2), simpsons_error)
PyPlot.title("Simpson's rule")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "simpsons_error_plot.png"), dpi=80)
PyPlot.figure(15)
PyPlot.clf()
PyPlot.plot(3*round.(nc/3), simpsons38_error)
PyPlot.title("Simpson's 3/8 rule")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "simpsons38_error_plot.png"), dpi=80)
PyPlot.figure(16)
PyPlot.clf()
PyPlot.plot(nc, trapezoidal_error)
PyPlot.title("Trapezoidal rule")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "trapezoidal_error_plot.png"), dpi=80)