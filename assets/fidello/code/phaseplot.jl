# This file was generated, do not modify it. # hide
using DifferentialEquations, Plots, LaTeXStrings
ENV["GKSwstype"] = "100"
g, m, l = 9.81, 1.0, 1.0

function pendulum!(du, u, p, t)
    du[1] = u[2]/(m*l^2)
    du[2] = -m*g*l*cos(u[1])
end

tspan = (0.0, 10.0)
ics = [(θ₀, pθ₀) for θ₀ in range(-π/2+0.1, π/2-0.1, length=9), pθ₀ in range(-0.1, 0.1, length=3)]

plot()
for u₀ in ics
    local sol = solve(ODEProblem(pendulum!, collect(u₀), tspan))
    plot!(sol, idxs=(1, 2), label=false)
end

xlabel!(L"\theta (\mathrm{rad})")
ylabel!(L"p_{\theta}(\mathrm{kg}\cdot\mathrm{m}^2\cdot\mathrm{rad}\cdot\mathrm{s}^{-1})")
title!("Simple pendulum phase plot")
savefig("Phase_plot.svg")