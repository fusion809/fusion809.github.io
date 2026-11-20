using DifferentialEquations, Plots, LaTeXStrings
ENV["GKSwstype"] = "100"
g, m, l = 9.81, 1.0, 1.0
rho, Cd, A = 2.0, 0.05, 5.0

function pendulum!(du, u, p, t)
    du[1] = u[2]/(m*l^2)
    du[2] = -m*g*l*cos(u[1]) - rho*Cd*A*abs(u[2])*u[2]/(2*m^2*l)
end

tspan = (0.0, 100.0)
theta0 = pi/2-0.01
ptheta0 = 0

plot()
solut = solve(ODEProblem(pendulum!, [theta0, ptheta0], tspan))
plot!(solut, idxs=(1, 2), label=false)

xlabel!(L"\theta (\mathrm{rad})")
ylabel!(L"p_{\theta}(\mathrm{kg}\cdot\mathrm{m}^2\cdot\mathrm{rad}\cdot\mathrm{s}^{-1})")
title!("Dampened pendulum phase plot")
savefig("Phase_plot_dampened.svg")
