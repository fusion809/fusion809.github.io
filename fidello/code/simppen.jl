using DifferentialEquations, Plots, LaTeXStrings
ENV["GKSwstype"] = "100"
g, l = 9.81, 1.0
fps = 30

function pendulum!(du, u, p, t)
    du[1] = u[2]
    du[2] = -(g/l)*cos(u[1])
end

tspan = (0.0, 10.0)
u0 = [0.0, 0.0]  # θ₀ = 0, θ̇₀ = 0

prob = ODEProblem(pendulum!, u0, tspan)
sol = solve(prob, saveat=1/fps)

plot(sol.t, sol[1, :], label=false)
xlabel!(L"t\;(\mathrm{seconds})")
ylabel!(L"\theta\;(\mathrm{radians})")
title!("Simple pendulum plot")
savefig("Simple_pendulum_theta.svg")
# Animate with θ measured counterclockwise from the positive x-axis
x = l .* cos.(sol[1, :])
y = l .* sin.(sol[1, :])

anim = @animate for i in eachindex(sol.t)
    plot(
        [0, x[i]], [0, y[i]];
        xlims=(-1.1l, 1.1l),
        ylims=(-1.1l, 1.1l),
        aspect_ratio=:equal,
        label=false,
        xlabel=L"x\ (\mathrm{metres})",
        ylabel=L"y\ (\mathrm{metres})",
        title="Simple pendulum",
    )
    scatter!([x[i]], [y[i]]; markersize=8, label=false)
end

gif(anim, "Simple_pendulum.gif"; fps=fps)