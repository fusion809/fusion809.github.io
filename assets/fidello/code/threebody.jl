# This file was generated, do not modify it. # hide
using DifferentialEquations
using Plots, Measures, Colors
using Printf, LaTeXStrings
ENV["GKSwstype"] = "100"
# ============================================================
# Three-body gravitational system
#
# Body M is fixed at the origin.
# Bodies m₁ and m₂ orbit M while also perturbing one another.
#
# Unlike a strongly hierarchical system, the masses here are
# comparable enough for the mutual interaction to produce
# strongly non-Keplerian and potentially chaotic motion.
#
# Units are normalised:
#   G = 1
#   mass = arbitrary mass unit
#   distance = arbitrary distance unit
#   time = corresponding dynamical time unit
# ============================================================

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

G              = 1.0
M              = 3.0
m₁             = 1.5
m₂             = 1.5
a₁             = 4.0
e₁             = 0.15
a₂             = 5.0
e₂             = 0.10
tspan          = (0.0, 280.0)
abstol         = 1e-11
reltol         = 1e-11
fps            = 60
marksize       = 10
stepsize       = 25/fps
legendfontsize = 24
guidefontsize  = 24
tickfontsize   = 20
titlefontsize  = 30

# ------------------------------------------------------------
# Initial approximate Keplerian orbits
# ------------------------------------------------------------

# Periapsis distances

r₁₀ = a₁ * (1 - e₁)
r₂₀ = a₂ * (1 - e₂)

# ------------------------------------------------------------
# Initial angular positions
# ------------------------------------------------------------

# Start on opposite sides of the central body. A small
# asymmetry is introduced so that the system does not remain
# artificially symmetric.

θ₁₀ = 0.0
θ₂₀ = π + 0.15

# ------------------------------------------------------------
# Initial Keplerian periapsis velocities
# ------------------------------------------------------------

v₁₀ = sqrt(G * M * (1 + e₁) / r₁₀)
v₂₀ = sqrt(G * M * (1 + e₂) / r₂₀)

# Radial velocity is zero at periapsis.

ṙ₁₀ = 0.0
ṙ₂₀ = 0.0

# Since v = r θ̇:

θ̇₁₀ = v₁₀ / r₁₀
θ̇₂₀ = v₂₀ / r₂₀

# ------------------------------------------------------------
# Initial state vector
#
# u = [
#     r₁,
#     r₂,
#     θ₁,
#     θ₂,
#     ṙ₁,
#     ṙ₂,
#     θ̇₁,
#     θ̇₂
# ]
# ------------------------------------------------------------

u₀ = [
    r₁₀,
    r₂₀,
    θ₁₀,
    θ₂₀,
    ṙ₁₀,
    ṙ₂₀,
    θ̇₁₀,
    θ̇₂₀
]

# ------------------------------------------------------------
# Equations of motion
# ------------------------------------------------------------

function three_body!(du, u, p, t)

    r₁, r₂, θ₁, θ₂, ṙ₁, ṙ₂, θ̇₁, θ̇₂ = u

    Δθ = θ₂ - θ₁

    r₁₂ = sqrt(
        r₁^2 +
        r₂^2 -
        2*r₁*r₂*cos(Δθ)
    )

    # --------------------------------------------------------
    # Radial accelerations
    # --------------------------------------------------------

    r̈₁ =
        r₁*θ̇₁^2 -
        G*M/r₁^2 -
        G*m₂*(r₁ - r₂*cos(Δθ))/r₁₂^3

    r̈₂ =
        r₂*θ̇₂^2 -
        G*M/r₂^2 -
        G*m₁*(r₂ - r₁*cos(Δθ))/r₁₂^3

    # --------------------------------------------------------
    # Angular accelerations
    # --------------------------------------------------------

    θ̈₁ =
        -2*ṙ₁*θ̇₁/r₁ +
        G*m₂*r₂*sin(Δθ)/(r₁*r₁₂^3)

    θ̈₂ =
        -2*ṙ₂*θ̇₂/r₂ -
        G*m₁*r₁*sin(Δθ)/(r₂*r₁₂^3)

    # --------------------------------------------------------
    # First-order system
    # --------------------------------------------------------

    du[1] = ṙ₁
    du[2] = ṙ₂
    du[3] = θ̇₁
    du[4] = θ̇₂

    du[5] = r̈₁
    du[6] = r̈₂
    du[7] = θ̈₁
    du[8] = θ̈₂
end

# ------------------------------------------------------------
# Integrate
#
# 100 time units gives the system ample time to undergo
# several orbital periods and develop strongly perturbed
# trajectories.
# ------------------------------------------------------------

prob = ODEProblem(
    three_body!,
    u₀,
    tspan
)

solut = solve(
    prob,
    Vern8(),
    saveat=stepsize,
    abstol=abstol,
    reltol=reltol
)

# ------------------------------------------------------------
# Convert polar coordinates to Cartesian coordinates
# ------------------------------------------------------------

x₁ = solut[1, :] .* cos.(solut[3, :])
y₁ = solut[1, :] .* sin.(solut[3, :])

x₂ = solut[2, :] .* cos.(solut[4, :])
y₂ = solut[2, :] .* sin.(solut[4, :])

# ------------------------------------------------------------
# Calculate body-body separation
# ------------------------------------------------------------

separation = sqrt.(
    (x₂ .- x₁).^2 +
    (y₂ .- y₁).^2
)
# ------------------------------------------------------------
# Plot limits
# ------------------------------------------------------------

x_all = vcat(0.0, x₁, x₂)
y_all = vcat(0.0, y₁, y₂)

xmin = minimum(x_all)
xmax = maximum(x_all)
ymin = minimum(y_all)
ymax = maximum(y_all)

span = max(
    xmax - xmin,
    ymax - ymin
)

margin = 0.08 * span

xc = (xmin + xmax) / 2
yc = (ymin + ymax) / 2

xmin = xc - span/2 - margin
xmax = xc + span/2 + margin

ymin = yc - span/2 - margin
ymax = yc + span/2 + margin

# ------------------------------------------------------------
# Animation
# ------------------------------------------------------------

m1col = RGB24(0, 0, 155/255)
m1outcol = RGB24(190/255, 190/255, 255/255)
m2col = RGB24(0, 155/255, 155/255)
m2outcol = RGB24(190/255, 255/255, 255/255)
Mcol  = RGB24(65/255, 130/255, 65/255)

anim = @animate for i in eachindex(solut.t)

    p = plot(
        x₁[1:i],
        y₁[1:i],
        linewidth=1,
        color=m1col,
        label=false,
        xlims=(xmin, xmax),
        ylims=(ymin, ymax),
        aspect_ratio=:equal,
        xlabel=L"x",
        ylabel=L"y",
        title="Three-body system — t = $(@sprintf("%.2f", solut.t[i]))",
        guidefontsize=guidefontsize,
        tickfontsize=tickfontsize,
        titlefontsize=titlefontsize,
        legend=:topleft,
        legendfontsize=legendfontsize,
        left_margin=2Plots.mm,
        right_margin=2Plots.mm,
        top_margin=2Plots.mm,
        bottom_margin=2Plots.mm,
        size=(1600, 1600)
    )

    plot!(
        p,
        x₂[1:i],
        y₂[1:i],
        linewidth=1,
        color=m2col,
        label=false
    )

    # Central body
    scatter!(
        p,
        [0.0],
        [0.0],
        markershape=:circle,
        markersize=marksize,
        color=Mcol,
        markerstrokecolor=Mcol,
        label="M"
    )

    # Body 1
    scatter!(
        p,
        [x₁[i]],
        [y₁[i]],
        markershape=:circle,
        markersize=round(marksize/2),
        color=m1col,
        markerstrokecolor=m1outcol,
        label="m₁"
    )

    # Body 2
    scatter!(
        p,
        [x₂[i]],
        [y₂[i]],
        markershape=:circle,
        markersize=round(marksize/2),
        color=m2col,
        markerstrokecolor=m2outcol,
        label="m₂"
    )
end

# ------------------------------------------------------------
# Save GIF
#
# saveat = 1/fps
# fps    = fps
#
# Therefore one simulated second corresponds to one second
# of animation playback.
# ------------------------------------------------------------

gif(
    anim,
    "three_body_system.gif",
    fps=fps
)