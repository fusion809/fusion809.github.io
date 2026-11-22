# This file was generated, do not modify it. # hide
#!/usr/bin/env julia

using Pkg

# Install missing packages automatically.
for package in ("FFTW", "GLMakie")
    if !haskey(Pkg.project().dependencies, package)
        Pkg.add(package)
    end
end

using FFTW
import GLMakie
using Printf

# ============================================================
# Configuration
# ============================================================

const T_ROOM = 20.0
const T_BALL = 100.0

const ALPHA = 1e-5
const RADIUS = 0.10

const LX = 1.0
const LY = 1.0
const LZ = 1.0

const NX = 64
const NY = 64
const NZ = 64

const SIMULATION_TIME = 1800.0

const NUMBER_OF_FRAMES = 600
const FRAME_RATE = 20

const OUTPUT_FILE = "heat_equation_3d.gif"

# Smooth transition between the hot region and room temperature.
const TRANSITION_CELLS = 6.0

# RKF45 tolerances.
const RELTOL = 1e-7
const ABSTOL = 1e-6

const INITIAL_DT = 0.05
const MAX_DT = 2.0
const MIN_DT = 1e-9

# Conservative explicit-RK stability restriction.
const STABILITY_FACTOR = 0.50

# ============================================================
# Spatial grid
# ============================================================

dx = LX / NX
dy = LY / NY
dz = LZ / NZ

x = collect((0:NX-1) .* dx .- LX / 2)
y = collect((0:NY-1) .* dy .- LY / 2)
z = collect((0:NZ-1) .* dz .- LZ / 2)

# ============================================================
# Smooth initial hot sphere
# ============================================================

const TRANSITION_WIDTH = TRANSITION_CELLS * max(dx, dy, dz)

function smootherstep(s)
    return s^3 * (10.0 + s * (-15.0 + 6.0 * s))
end

function initial_temperature(r)
    inner_radius = RADIUS - TRANSITION_WIDTH / 2
    outer_radius = RADIUS + TRANSITION_WIDTH / 2

    if r <= inner_radius
        return T_BALL
    elseif r >= outer_radius
        return T_ROOM
    else
        s = (r - inner_radius) / TRANSITION_WIDTH
        weight = 1.0 - smootherstep(s)

        return T_ROOM + (T_BALL - T_ROOM) * weight
    end
end

T = Array{Float64}(undef, NX, NY, NZ)

for k in 1:NZ
    for j in 1:NY
        for i in 1:NX
            r = sqrt(x[i]^2 + y[j]^2 + z[k]^2)
            T[i, j, k] = initial_temperature(r)
        end
    end
end

# ============================================================
# Fourier spectral Laplacian
# ============================================================

kx = 2π .* FFTW.fftfreq(NX, 1 / dx)
ky = 2π .* FFTW.fftfreq(NY, 1 / dy)
kz = 2π .* FFTW.fftfreq(NZ, 1 / dz)

KX = reshape(kx, NX, 1, 1)
KY = reshape(ky, 1, NY, 1)
KZ = reshape(kz, 1, 1, NZ)

lambda = -ALPHA .* (
    KX.^2 .+
    KY.^2 .+
    KZ.^2
)

function heat_rhs(u, lambda)
    return real.(ifft(lambda .* fft(u)))
end

# ============================================================
# RKF45
# ============================================================

function rkf45_step(u, h, lambda)

    k1 = heat_rhs(u, lambda)

    k2 = heat_rhs(
        u .+ h .* (1 / 4) .* k1,
        lambda
    )

    k3 = heat_rhs(
        u .+ h .* (
            (3 / 32) .* k1 .+
            (9 / 32) .* k2
        ),
        lambda
    )

    k4 = heat_rhs(
        u .+ h .* (
            (1932 / 2197) .* k1 .-
            (7200 / 2197) .* k2 .+
            (7296 / 2197) .* k3
        ),
        lambda
    )

    k5 = heat_rhs(
        u .+ h .* (
            (439 / 216) .* k1 .-
            8 .* k2 .+
            (3680 / 513) .* k3 .-
            (845 / 4104) .* k4
        ),
        lambda
    )

    k6 = heat_rhs(
        u .+ h .* (
            -(8 / 27) .* k1 .+
            2 .* k2 .-
            (3544 / 2565) .* k3 .+
            (1859 / 4104) .* k4 .-
            (11 / 40) .* k5
        ),
        lambda
    )

    # Fourth-order estimate.
    u4 = u .+ h .* (
        (25 / 216) .* k1 .+
        (1408 / 2565) .* k3 .+
        (2197 / 4104) .* k4 .-
        (1 / 5) .* k5
    )

    # Fifth-order estimate.
    u5 = u .+ h .* (
        (16 / 135) .* k1 .+
        (6656 / 12825) .* k3 .+
        (28561 / 56430) .* k4 .-
        (9 / 50) .* k5 .+
        (2 / 55) .* k6
    )

    return u5, u4
end

# ============================================================
# Stable timestep restriction
# ============================================================

lambda_max_abs = maximum(abs.(lambda))

dt_stability = STABILITY_FACTOR / lambda_max_abs

@printf(
    "Maximum |lambda| = %.8f s^-1\n",
    lambda_max_abs
)

@printf(
    "RKF45 stability timestep limit = %.8f s\n",
    dt_stability
)

# ============================================================
# Adaptive integration
# ============================================================

function integrate_to!(u, t, target, h, lambda)

    while t < target

        h = min(
            h,
            MAX_DT,
            dt_stability,
            target - t
        )

        if h < MIN_DT
            error(
                "RKF45 timestep fell below MIN_DT at t = $t"
            )
        end

        u5, u4 = rkf45_step(u, h, lambda)

        scale =
            ABSTOL .+
            RELTOL .* max.(abs.(u), abs.(u5))

        error_norm =
            maximum(
                abs.(u5 .- u4) ./ scale
            )

        if error_norm <= 1.0
            u .= u5
            t += h
        end

        factor =
            if error_norm == 0.0
                5.0
            else
                clamp(
                    0.9 * error_norm^(-0.2),
                    0.2,
                    5.0
                )
            end

        h = min(
            MAX_DT,
            dt_stability,
            h * factor
        )
    end

    return t, h
end

# ============================================================
# GLMakie visualisation
# ============================================================

GLMakie.activate!()

fig = GLMakie.Figure(
    size = (1000, 800)
)

ax = GLMakie.Axis3(
    fig[1, 1],
    aspect = (1, 1, 1),
    xlabel = "",
    ylabel = "",
    zlabel = "",
    title = ""
)

# ============================================================
# Temperature colour map
#
# 20 °C is opaque white, matching the background.
# ============================================================

temperature_colormap = [
    GLMakie.RGBAf(1.0, 1.0, 1.0, 1.0),  # 20 °C
    GLMakie.RGBAf(0.05, 0.10, 0.80, 1.0),
    GLMakie.RGBAf(0.0, 0.65, 1.0, 1.0),
    GLMakie.RGBAf(0.0, 0.90, 0.45, 1.0),
    GLMakie.RGBAf(1.0, 0.90, 0.0, 1.0),
    GLMakie.RGBAf(1.0, 0.35, 0.0, 1.0),
    GLMakie.RGBAf(0.85, 0.0, 0.0, 1.0)
]

# Makie expects volume dimensions as two endpoints.
volume_plot = GLMakie.volume!(
    ax,
    (x[1], x[end]),
    (y[1], y[end]),
    (z[1], z[end]),
    T;
    algorithm = :mip,
    colormap = temperature_colormap,
    colorrange = (T_ROOM, T_BALL)
)

# ============================================================
# Colourbar
# ============================================================

GLMakie.Colorbar(
    fig[1, 2],
    volume_plot;
    label = "Temperature (°C)",
    ticks = (
        collect(20.0:10.0:100.0),
        string.(collect(20.0:10.0:100.0))
    )
)

# ============================================================
# Simulation time display
# ============================================================

time_label = GLMakie.Label(
    fig[2, 1],
    "Simulation time: 0.00 s";
    tellwidth = false,
    fontsize = 20
)

# ============================================================
# Animation
# ============================================================

frame_times = collect(
    range(
        0.0,
        SIMULATION_TIME;
        length = NUMBER_OF_FRAMES
    )
)

u = copy(T)

# Keep mutable animation state local to the record callback.
let t = 0.0,
    dt = min(INITIAL_DT, dt_stability)

    GLMakie.record(
        fig,
        OUTPUT_FILE,
        enumerate(frame_times);
        framerate = FRAME_RATE
    ) do (frame_number, target_time)

        if frame_number > 1

            t, dt = integrate_to!(
                u,
                t,
                target_time,
                dt,
                lambda
            )

        end

        # Update the volume with the actual temperature field.
        volume_plot[4] = u

        # Update simulation time.
        elapsed_minutes = floor(Int, target_time / 60)
        remaining_seconds = target_time - 60 * elapsed_minutes
        time_label.text =
            @sprintf(
                "Simulation time: %d min %.1f s",
                elapsed_minutes,
                remaining_seconds
            )

        # Console diagnostics.
        Tmin = minimum(u)
        Tmax = maximum(u)

        @printf(
            "Frame %4d/%d | t = %d min %.1f s | dt = %.5f s | T range = %.8f–%.8f °C\n",
            frame_number,
            NUMBER_OF_FRAMES,
            elapsed_minutes,
            remaining_seconds,
            dt,
            Tmin,
            Tmax
        )
    end
end

println()
println("Animation written to: $OUTPUT_FILE")