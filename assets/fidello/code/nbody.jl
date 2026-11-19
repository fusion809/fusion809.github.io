# This file was generated, do not modify it. # hide
using DifferentialEquations
using Plots, Measures, Colors
using Printf, LaTeXStrings
ENV["GKSwstype"] = "100"

# ============================================================
# N-body gravitational system with a fixed central mass
#
# M is fixed at the origin.
# N orbiting bodies interact gravitationally with M and with
# one another.
#
# State vector:
# u = [r₁,...,rₙ, θ₁,...,θₙ, ṙ₁,...,ṙₙ, θ̇₁,...,θ̇ₙ]
#
# All quantities use normalized units.
# ============================================================

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
# Add or remove entries to change the number of orbiting bodies.
# All parameter vectors must have the same length.
G              = 1.0
M              = 3.0
masses         = [1.5, 2.0, 0.01]
N              = length(masses)
# Initial conditions
## Semi-major axes
a              = [5.0, 5.0, 15.0]
## Eccentricities
e              = [0.20, 0.05, 0.01]
## Initial angular positions
θ₀             = [
    0.0
    π - 0.15
    π/3 + 0.06
]
# Frames per second of animation
fps            = 60
# Solution params
tspan          = (0.0, 60.0)
abstol         = 1e-11
reltol         = 1e-11
stepsize       = 8/fps
# Animation params
marksize       = 10
legendfontsize = 24
guidefontsize  = 24
tickfontsize   = 20
titlefontsize  = 30

# ------------------------------------------------------------
# Validate parameters
# ------------------------------------------------------------
@assert N >= 1 "At least one orbiting body is required."
@assert length(a) == N "The semimajor axis vector must have N entries."
@assert length(e) == N "The eccentricity vector must have N entries."
@assert all(>(0), masses) "All orbiting masses must be positive."
@assert all(>(0), a) "All semimajor axes must be positive."
@assert all(x -> 0 <= x < 1, e) "Eccentricities must satisfy 0 ≤ e < 1."

# ------------------------------------------------------------
# Initial approximate Keplerian orbits
# ------------------------------------------------------------

# Periapsis distances

r₀ = a .* (1 .- e)

# Initial Keplerian periapsis velocities

v₀ = sqrt.(G .* M .* (1 .+ e) ./ r₀)

# Radial velocities vanish at periapsis.

ṙ₀ = zeros(N)

# Since v = r θ̇:

θ̇₀ = v₀ ./ r₀

# ------------------------------------------------------------
# Initial state vector
# ------------------------------------------------------------

u₀ = vcat(
    r₀,
    θ₀,
    ṙ₀,
    θ̇₀
)

# ------------------------------------------------------------
# Equations of motion
# ------------------------------------------------------------

function n_body!(du, u, p, t)

    G, M, masses, N = p

    # Obtain views of the state vector.

    r  = @view u[1:N]
    θ  = @view u[N+1:2N]
    ṙ = @view u[2N+1:3N]
    θ̇ = @view u[3N+1:4N]

    # Obtain views of the derivative vector.

    dr  = @view du[1:N]
    dθ  = @view du[N+1:2N]
    dṙ = @view du[2N+1:3N]
    dθ̇ = @view du[3N+1:4N]

    # First-order kinematic equations.

    dr .= ṙ
    dθ .= θ̇

    # --------------------------------------------------------
    # Gravitational accelerations
    # --------------------------------------------------------

    for k in 1:N

        # Angular differences relative to body k.

        Δθ = θ .- θ[k]

        # Distances between body k and every other body.

        rᵢₖ = sqrt.(
            r[k]^2 .+
            r.^2 .-
            2 .* r[k] .* r .* cos.(Δθ)
        )

        # Exclude self-interaction.

        rᵢₖ[k] = Inf

        # ----------------------------------------------------
        # Radial acceleration
        # ----------------------------------------------------

        radial_interaction = sum(
            masses .* (
                r[k] .- r .* cos.(Δθ)
            ) ./ rᵢₖ.^3
        )

        r̈ₖ =
            r[k]*θ̇[k]^2 -
            G*M/r[k]^2 -
            G*radial_interaction

        # ----------------------------------------------------
        # Angular acceleration
        # ----------------------------------------------------

        angular_interaction = sum(
            masses .* r .* sin.(Δθ) ./ rᵢₖ.^3
        )

        θ̈ₖ =
            -2*ṙ[k]*θ̇[k]/r[k] +
            G*angular_interaction/r[k]

        dṙ[k] = r̈ₖ
        dθ̇[k] = θ̈ₖ

    end
end

# ------------------------------------------------------------
# Integrate
# ------------------------------------------------------------

p = (G, M, masses, N)

prob = ODEProblem(
    n_body!,
    u₀,
    tspan,
    p
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

x = [
    solut[k, :] .* cos.(solut[N+k, :])
    for k in 1:N
]

y = [
    solut[k, :] .* sin.(solut[N+k, :])
    for k in 1:N
]

# ------------------------------------------------------------
# Calculate pairwise separations
# ------------------------------------------------------------

separations = Dict{Tuple{Int,Int},Vector{Float64}}()

for k in 1:N-1
    for i in k+1:N

        separations[(k,i)] = sqrt.(
            (x[k] .- x[i]).^2 .+
            (y[k] .- y[i]).^2
        )

    end
end

# ------------------------------------------------------------
# Plot limits
# ------------------------------------------------------------

x_all = vcat([0.0], x...)
y_all = vcat([0.0], y...)

xmin = minimum(x_all)
xmax = maximum(x_all)
ymin = minimum(y_all)
ymax = maximum(y_all)

span = max(xmax - xmin, ymax - ymin)

margin = 0.08 * span

xc = (xmin + xmax)/2
yc = (ymin + ymax)/2

xmin = xc - span/2 - margin
xmax = xc + span/2 + margin

ymin = yc - span/2 - margin
ymax = yc + span/2 + margin

# ------------------------------------------------------------
# Colours
# ------------------------------------------------------------

palette = [
    RGB(0.0, 0.0, 0.65),   # Dark blue
    RGB(0.0, 0.55, 0.55),  # Dark teal
    RGB(0.65, 0.0, 0.0),   # Dark red
    RGB(0.45, 0.0, 0.65),  # Purple
    RGB(0.0, 0.45, 0.15),  # Dark green
    RGB(0.65, 0.25, 0.0),  # Dark orange
    RGB(0.25, 0.25, 0.25), # Dark grey
    RGB(0.0, 0.4, 0.65)    # Steel blue
]

bodycolors = [palette[mod1(k, length(palette))] for k in 1:N]

centralcolor = RGB24(65/255, 130/255, 65/255)

# ------------------------------------------------------------
# Animation
# ------------------------------------------------------------

anim = @animate for frame in eachindex(solut.t)

    p = plot(
        xlims=(xmin, xmax),
        ylims=(ymin, ymax),
        aspect_ratio=:equal,
        xlabel=L"x",
        ylabel=L"y",
        title="N-body system — t = $(@sprintf("%.2f", solut.t[frame]))",
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

    # Plot trajectories and current positions.

    for k in 1:N

        plot!(
            p,
            x[k][1:frame],
            y[k][1:frame],
            linewidth=1,
            color=bodycolors[k],
            label=false
        )

        scatter!(
            p,
            [x[k][frame]],
            [y[k][frame]],
            markershape=:circle,
            markersize=round(marksize/2),
            color=bodycolors[k],
            markerstrokecolor=bodycolors[k],
            label=latexstring("m_{$k}")
        )

    end

    # Central body.

    scatter!(
        p,
        [0.0],
        [0.0],
        markershape=:circle,
        markersize=marksize,
        color=centralcolor,
        markerstrokecolor=centralcolor,
        label="M"
    )

end

# ------------------------------------------------------------
# Save GIF
#
# saveat = 1/fps
# fps    = fps
#
# One simulated second corresponds to one second of playback.
# ------------------------------------------------------------

gif(
    anim,
    "n_body_system.gif",
    fps=fps
)