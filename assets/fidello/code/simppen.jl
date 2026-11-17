# This file was generated, do not modify it. # hide
using DifferentialEquations, Plots, LaTeXStrings

rho, Cd, A, g, m, l = 2.0, 0.01, 5.0, 9.81, 1.0, 1.0

function pendulum!(du, u, p, t)
    du[1] = u[2]
    du[2] = -(g/l)*cos(u[1])-rho*Cd*A*l*abs(u[2])*u[2]/(2*m)
end

tspan = (0.0, 10.0)
u0 = [0.0, 0.0]  # θ₀ = 0, θ̇₀ = 0

prob = ODEProblem(pendulum!, u0, tspan)
sol = solve(prob, saveat=0.01)

plot(sol.t, sol[1, :], label=false)
xlabel!(L"t\;(\mathrm{seconds})")
ylabel!(L"\theta\;(\mathrm{radians})")
title!("Dampened pendulum")
savefig("Dampened_pendulum_theta.svg")