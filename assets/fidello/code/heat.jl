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

const ALPHA = 1.0e-5
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

const GIF_FILE = "heat_equation_3d.gif"

# Only display regions at least this many degrees
# above room temperature.
#
# This affects visualisation only.
const VISUAL_THRESHOLD = 1.0

# Width of the smooth transition at the edge of the
# initially hot sphere, measured in grid cells.
#
# A wider transition substantially reduces high-frequency
# Fourier components and therefore reduces Gibbs oscillations.
const TRANSITION_CELLS = 6.0

# Adaptive RKF45 tolerances.
const RELTOL = 1.0e-7
const ABSTOL = 1.0e-6

const INITIAL_DT = 0.05
const MAX_DT = 2.0
const MIN_DT = 1.0e-9

# Conservative fraction of the RKF45 negative-real-axis
# stability limit used for the diffusion timestep.
#
# This is deliberately conservative because the heat equation
# is a stiff system when represented spectrally.
const STABILITY_FACTOR = 0.50

# ============================================================
# Spatial grid
# ============================================================

dx = LX / NX
dy = LY / NY
dz = LZ / NZ

x = collect((0:NX-1) .* dx .- LX/2)
y = collect((0:NY-1) .* dy .- LY/2)
z = collect((0:NZ-1) .* dz .- LZ/2)

const TRANSITION_WIDTH =
    TRANSITION_CELLS * max(dx, dy, dz)

# ============================================================
# Fourier spectral Laplacian
# ============================================================

kx = 2π .* FFTW.fftfreq(NX, 1/dx)
ky = 2π .* FFTW.fftfreq(NY, 1/dy)
kz = 2π .* FFTW.fftfreq(NZ, 1/dz)

KX = reshape(kx, NX, 1, 1)
KY = reshape(ky, 1, NY, 1)
KZ = reshape(kz, 1, 1, NZ)

lambda = -ALPHA .* (
    KX.^2 .+
    KY.^2 .+
    KZ.^2
)

# Largest magnitude eigenvalue of the discrete Fourier
# representation of the Laplacian.
lambda_max_abs = maximum(abs.(lambda))

# Conservative timestep based on the fastest Fourier mode.
#
# RKF45 has a finite stability region on the negative real axis.
# Using a conservative factor here prevents the adaptive error
# controller from taking an apparently accurate but unstable
# timestep for rapidly decaying high-frequency modes.
dt_stability = STABILITY_FACTOR / lambda_max_abs

# ============================================================
# Smooth initial spherical temperature distribution
# ============================================================

# Quintic smootherstep:
#
# f(s) = 6s^5 - 15s^4 + 10s^3
#
# f(0) = 0
# f(1) = 1
#
# The first and second derivatives vanish at both endpoints.

function smootherstep(s)

    return s^3 * (
        10.0 +
        s * (
            -15.0 +
            6.0*s
        )
    )
end

function initial_temperature(r)

    inner_radius =
        RADIUS - TRANSITION_WIDTH / 2

    outer_radius =
        RADIUS + TRANSITION_WIDTH / 2

    if r <= inner_radius

        return T_BALL

    elseif r >= outer_radius

        return T_ROOM

    else

        s = (
            r - inner_radius
        ) / TRANSITION_WIDTH

        weight = 1.0 - smootherstep(s)

        return T_ROOM +
               (T_BALL - T_ROOM) * weight
    end
end

T = Array{Float64}(undef, NX, NY, NZ)

for k in 1:NZ, j in 1:NY, i in 1:NX

    r = sqrt(
        x[i]^2 +
        y[j]^2 +
        z[k]^2
    )

    T[i,j,k] = initial_temperature(r)
end

# ============================================================
# Heat equation right-hand side
# ============================================================

function heat_rhs(u, lambda)

    return real.(
        ifft(
            lambda .* fft(u)
        )
    )
end

# ============================================================
# Runge–Kutta–Fehlberg 4(5)
# ============================================================

function rkf45_step(u, h, lambda)

    k1 = heat_rhs(
        u,
        lambda
    )

    k2 = heat_rhs(
        u .+
        h .* (1/4) .* k1,
        lambda
    )

    k3 = heat_rhs(
        u .+
        h .* (
            (3/32) .* k1 .+
            (9/32) .* k2
        ),
        lambda
    )

    k4 = heat_rhs(
        u .+
        h .* (
            (1932/2197) .* k1 .-
            (7200/2197) .* k2 .+
            (7296/2197) .* k3
        ),
        lambda
    )

    k5 = heat_rhs(
        u .+
        h .* (
            (439/216) .* k1 .-
            8 .* k2 .+
            (3680/513) .* k3 .-
            (845/4104) .* k4
        ),
        lambda
    )

    k6 = heat_rhs(
        u .+
        h .* (
            -(8/27) .* k1 .+
            2 .* k2 .-
            (3544/2565) .* k3 .+
            (1859/4104) .* k4 .-
            (11/40) .* k5
        ),
        lambda
    )

    u4 = u .+
        h .* (
            (25/216) .* k1 .+
            (1408/2565) .* k3 .+
            (2197/4104) .* k4 .-
            (1/5) .* k5
        )

    u5 = u .+
        h .* (
            (16/135) .* k1 .+
            (6656/12825) .* k3 .+
            (28561/56430) .* k4 .-
            (9/50) .* k5 .+
            (2/55) .* k6
        )

    return u5, u4
end

# ============================================================
# Adaptive integration
# ============================================================

function integrate_to!(
    u,
    t,
    target,
    h,
    lambda,
    dt_stability
)

    while t < target

        # Apply both the requested maximum timestep and
        # the diffusion stability restriction.
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

        u5, u4 = rkf45_step(
            u,
            h,
            lambda
        )

        scale =
            ABSTOL .+
            RELTOL .* max.(
                abs.(u),
                abs.(u5)
            )

        error_norm = maximum(
            abs.(u5 .- u4) ./ scale
        )

        if error_norm <= 1.0

            u .= u5
            t += h
        end

        factor = if error_norm == 0.0

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
# Three-dimensional visualisation
# ============================================================

GLMakie.activate!()

# Temperature excess above room temperature.
#
# The threshold is applied only to the rendered field.
# T itself is never modified.

deltaT = max.(
    T .- T_ROOM,
    0.0
)

deltaT .= ifelse.(
    deltaT .>= VISUAL_THRESHOLD,
    deltaT,
    0.0
)

deltaT_obs =
    GLMakie.Observable(
        copy(deltaT)
    )

# Original colour map.
heat_colormap = [

    GLMakie.RGBAf(
        0.0, 0.0, 0.0, 0.0
    ),

    GLMakie.RGBAf(
        0.05, 0.10, 0.80, 1.0
    ),

    GLMakie.RGBAf(
        0.0, 0.65, 1.0, 1.0
    ),

    GLMakie.RGBAf(
        0.0, 0.90, 0.45, 1.0
    ),

    GLMakie.RGBAf(
        1.0, 0.90, 0.0, 1.0
    ),

    GLMakie.RGBAf(
        1.0, 0.35, 0.0, 1.0
    ),

    GLMakie.RGBAf(
        0.85, 0.0, 0.0, 1.0
    )
]

fig = GLMakie.Figure(
    size = (1100, 900)
)

ax = GLMakie.Axis3(
    fig[1,1],
    xlabel = "x (m)",
    ylabel = "y (m)",
    zlabel = "z (m)",
    aspect = :data
)

GLMakie.volume!(
    ax,
    (-LX/2, LX/2),
    (-LY/2, LY/2),
    (-LZ/2, LZ/2),
    deltaT_obs;
    algorithm = :mip,
    colormap = heat_colormap,
    colorrange = (
        0.0,
        T_BALL - T_ROOM
    )
)

# Display actual temperatures on the colourbar.
GLMakie.Colorbar(
    fig[1,2],
    colormap = heat_colormap,
    limits = (
        0.0,
        T_BALL - T_ROOM
    ),
    ticks = (
        collect(
            0.0:10.0:80.0
        ),
        string.(
            collect(
                20.0:10.0:100.0
            )
        )
    ),
    label = "Temperature (°C)",
    height = GLMakie.Relative(0.75)
)

time_label = GLMakie.Label(
    fig[2,1:2],
    "Simulation time: 0.0 s | Maximum temperature: 100.0 °C",
    fontsize = 20
)

# ============================================================
# Prepare simulation state
# ============================================================

println()

println(
    "Starting heat equation simulation..."
)

println(
    "Grid: $(NX) × $(NY) × $(NZ)"
)

println(
    "Simulation duration: $(SIMULATION_TIME) seconds"
)

println(
    "Output frames: $(NUMBER_OF_FRAMES)"
)

println(
    "GIF playback duration: " *
    "$(NUMBER_OF_FRAMES / FRAME_RATE) seconds"
)

println(
    "Visualisation threshold: " *
    "$(VISUAL_THRESHOLD) °C"
)

println(
    "Smooth transition width: " *
    "$(TRANSITION_WIDTH) m " *
    "($(TRANSITION_CELLS) grid cells)"
)

println(
    @sprintf(
        "Maximum Fourier eigenvalue magnitude: %.8f s^-1",
        lambda_max_abs
    )
)

println(
    @sprintf(
        "Diffusion stability timestep: %.8f s",
        dt_stability
    )
)

println()

println(
    @sprintf(
        "Initial temperature range: %.8f–%.8f °C",
        minimum(T),
        maximum(T)
    )
)

println(
    @sprintf(
        "Initial mean temperature: %.8f °C",
        sum(T) / length(T)
    )
)

# Persistent integration state.
t_state = Ref(0.0)
h_state = Ref(INITIAL_DT)

output_times = range(
    0.0,
    SIMULATION_TIME,
    length = NUMBER_OF_FRAMES
)

# ============================================================
# Generate the animated GIF
# ============================================================

GLMakie.record(
    fig,
    GIF_FILE,
    1:NUMBER_OF_FRAMES;
    framerate = FRAME_RATE
) do frame

    # Target time for this frame.
    target_time =
        output_times[frame]

    # Advance the physical temperature field.
    new_t, new_h =
        integrate_to!(
            T,
            t_state[],
            target_time,
            h_state[],
            lambda,
            dt_stability
        )

    t_state[] = new_t
    h_state[] = new_h

    # Calculate temperature excess for visualisation.
    #
    # T itself remains untouched.
    @. deltaT =
        max(
            T - T_ROOM,
            0.0
        )

    @. deltaT =
        ifelse(
            deltaT >= VISUAL_THRESHOLD,
            deltaT,
            0.0
        )

    # Update the volume-rendering observable.
    copyto!(
        deltaT_obs[],
        deltaT
    )

    GLMakie.notify(
        deltaT_obs
    )

    # Calculate maximum physical temperature.
    maximum_temperature =
        maximum(T)

    minimum_temperature =
        minimum(T)

    # Update the timer and temperature display.
    time_label.text[] =
        @sprintf(
            "Simulation time: %.1f s (%.1f min) | Maximum temperature: %.4f °C",
            t_state[],
            t_state[] / 60,
            maximum_temperature
        )

    # Print progress and diagnostics.
    println(
        @sprintf(
            "Frame %3d/%3d | t = %8.2f s | dt = %.5f s | T range = %.8f–%.8f °C",
            frame,
            NUMBER_OF_FRAMES,
            t_state[],
            h_state[],
            minimum_temperature,
            maximum_temperature
        )
    )
end

println()
println(
    "GIF successfully generated: $GIF_FILE"
)