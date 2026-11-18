# This file was generated, do not modify it. # hide
using DifferentialEquations
using Plots
using LaTeXStrings

# Constants
g = 9.81
m = 1.0
ρ = 1.225
C_D = 0.47
A = 0.01
fps = 30

# Equations of motion
function projectile!(du, u, p, t)
    x, y, px, py = u

    v = sqrt(px^2 + py^2)

    du[1] = px / m
    du[2] = py / m
    du[3] = -(ρ * C_D * A * px * v) / (2m^2)
    du[4] = -m*g - (ρ * C_D * A * py * v) / (2m^2)
end

# Ground condition
function ground_condition(u, t, integrator)
    u[2]
end

# Terminate when the projectile crosses y = 0 while descending
function ground_affect!(integrator)
    terminate!(integrator)
end

ground_callback = ContinuousCallback(
    ground_condition,
    nothing,           # Do nothing when crossing y = 0 upwards
    ground_affect!    # Terminate when crossing y = 0 downwards
)

# Initial conditions
x₀ = 0.0
y₀ = 0.0
vx₀ = 30.0
vy₀ = 30.0

u₀ = [x₀, y₀, m*vx₀, m*vy₀]

# Integrate
tspan = (0.0, 10.0)

prob = ODEProblem(projectile!, u₀, tspan)

sol = solve(
    prob,
    callback=ground_callback,
    saveat=1/fps
)

# Animation limits
xmin = minimum(sol[1, :])
xmax = maximum(sol[1, :])
ymax = maximum(sol[2, :]) * 1.1

# Animation
anim = @animate for i in eachindex(sol.t)

    plot(
        sol[1, 1:i],
        sol[2, 1:i],
        xlabel=L"x\;(\mathrm{metres})",
        ylabel=L"y\;(\mathrm{metres})",
        title="Projectile motion",
        label=false,
        aspect_ratio=:equal,
        xlims=(xmin, xmax),
        ylims=(0, ymax)
    )

    scatter!(
        [sol[1, i]],
        [sol[2, i]],
        label=false,
        markersize=6
    )
end

gif(anim, "projectile_motion.gif", fps=fps)