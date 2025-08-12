# This file was generated, do not modify it. # hide
using FFTW, StaticArrays, Plots, Interpolations

"""
    RKF45(f::Function, params, t0::Number, tf::Number, conds::Vector, epsilon, dtInitial, tolType="absolute", dtMin=(tf-t0)/1e8)

f should be a function that returns the RHS of the ODE being solved expressed as a system of first-order ODEs
in an array of the form `[dx1/dt, dx2/dt, dx3/dt, dx4/dt, ..., dxn/dt]`. Its arguments should be: params 
(an object containing problem parameters), t (a Float64) and `vars::Vector` (a column vector of the form 
`[element1; element2; element3; ...; elementn]`).

`params` should be a named tuple containing parameter values. For the simple
pendulum problem with pendulum length 1 metre, for example, it can be written
as:

`params = (g = 9.81, l = 1.0)`.

`t0` is the value of t (the independent variable) at the beginning of the integration.

`tf` is the value of t at the end of the integration.

`conds` is an SVector containing initial conditions. For the simple pendulum
problem, for example, the following code may be used (where theta0 and
thetaDot0) have been defined elsewhere:

`conds = @SVector [theta0, thetaDot0]`.

`epsilon` is the error tolerance for the problem.

`dtInitial` is the initial guess for the step size.

`tolType` is the type of tolerance used. Either "absolute" or "relative".

`dtMin` is the minimum step size allowed. Default is (tf-t0)/1e8.
"""
function RKF45(f::Function, params::NamedTuple, t0::Float64, tf::Float64, conds::SVector, epsilon::Float64, dtInitial::Float64, tolType::String = "absolute", dtMin::Float64 = (tf-t0)/1e8)
    # Initialize relevant variables
    dt = dtInitial;
    t = Float64[t0];
    vars = [conds];
    i = 1;
    ti = t0;

    # Loop over t under the solution for tf has been found
    while ( ti < tf )
        varsi =  vars[i];
        dt = minimum((dt, tf-ti));

        # RKF45 approximators
        K1 = dt*f(params, ti, varsi);
        K2 = dt*f(params, ti + dt/4, varsi + K1/4);
        K3 = dt*f(params, ti + 3*dt/8, varsi + 3*K1/32 + 9*K2/32);
        K4 = dt*f(params, ti + 12*dt/13, varsi + 1932*K1/2197 - 7200*K2/2197 + 7296*K3/2197);
        K5 = dt*f(params, ti + dt, varsi + 439*K1/216 - 8*K2 + 3680*K3/513 - 845*K4/4104);
        K6 = dt*f(params, ti + dt/2, varsi - 8*K1/27 + 2*K2 - 3544*K3/2565 + 1859*K4/4104 - 11*K5/40);

        # 4/5th order approximations to next step value
        vars1 = varsi + 25*K1/216 + 1408*K3/2565 + 2197*K4/4104 - K5/5;
        vars2 = varsi + 16*K1/135 + 6656*K3/12825 + 28561*K4/56430 - 9*K5/50 + 2*K6/55;

        # Determine if error is small enough to move on to next step
        if (tolType in ["relative", "rel", "R", "r", "Rel", "Relative"])
            R = maximum(abs.(vars2 - vars1)./abs.(vars1))/dt;
        elseif (tolType in ["absolute", "abs", "A", "a", "Abs", "Absolute"])
            R = maximum(abs.(vars2 - vars1))/dt;
        else
            error("tolType is set to an invalid value ($tolType), so exiting...")
        end
        s = (epsilon/(2*R))^(0.25);
        if (R <= epsilon)
            Base.push!(t, ti+dt);
            StaticArrays.push!(vars, vars1);
            i += 1;
            ti = t[i];
        end
        dt *= s;
        if (dt < dtMin)
            @warn("dt has reached $dt at t=$ti which is less than dtMin=$dtMin")
            if (tolType in ["absolute", "abs", "A", "a", "Abs", "Absolute"])
                tolType = "relative";
                @warn("As you are using an absolute tolerance type, we will switch to relative tolerance to see if this fixes the problem...")
            elseif (tolType in ["relative", "rel", "R", "r", "Rel", "Relative"])
                @warn("Breaking out of loop as tolerance type is already set to relative.")
                break
            else
                error("tolType is set to an invalid value ($tolType), so exiting...")
            end
        end
    end

    # Transpose and enter into NamedTuple
    vars = transpose(reduce(hcat, vars));
    return t, vars;
end
function heat(params, t, u::SVector)
    α = params.α

    # Wavenumbers
    k = params.k

    # FFT of u
    u_hat = fft(collect(u))  # convert to Vector for FFTW

    # Second derivative in Fourier space: -(k^2) * u_hat
    uxx_hat = - (k .^ 2) .* u_hat

    # Back to real space
    u_xx = real(ifft(uxx_hat))

    # Return du/dt as SVector for RKF45
    return SVector{length(u_xx)}(α .* u_xx)
end

N = 128                # number of grid points
T = 2π                 # period
L = T                  # domain length
dx = L / N
α = 0.01               # thermal diffusivity
# Wavenumbers for FFT (assumes N even)
k = vcat(0:N÷2-1, 0, -N÷2+1:-1) .* (2π/L)

# Initial condition
x = dx .* (0:N-1)
u0 = 10 * (2 .- cos.(2*pi/T * x))

params = (α=α, L=L, k=k)
conds = SVector{length(u0)}(u0)
t0 = 0.0
tf = 300.0
epsilon = 1e-9
dtInitial = 1e-3

t_vals, u_vals = RKF45(heat, params, t0, tf, conds, epsilon, dtInitial)
Umat = Matrix(u_vals)';
plotlyjs()  # use PlotlyJS for 3D plotting

# Create grids for x and t matching Umat dimensions
X = repeat(x, 1, length(t_vals))           # N × nt
T = repeat(t_vals', length(x), 1)          # N × nt (t_vals' transposed to row vector)

# Surface plot
surface(X, T, Umat, xlabel="x", ylabel="t", zlabel="u(x,t)",
        title="Heat equation", legend=false)
savefig(joinpath(@OUTPUT, "HeatEquation.svg"))

# Animate
nt_uniform = Int64(round((tf-t0)*30));
t_uniform = range(t_vals[1], t_vals[end], length=nt_uniform)  # e.g., 300 frames

U_interp = zeros(size(Umat, 1), length(t_uniform))

for i in 1:size(Umat, 1)
    itp = LinearInterpolation(t_vals, Umat[i, :])
    U_interp[i, :] = itp.(t_uniform)
end

function fixed_width_time(t; digits=3, width=7)
    # Format time with fixed decimal places
    s = string(round(t, digits=digits))
    # Pad with spaces on left so total length is 'width'
    return lpad(s, width)
end

ymin, ymax = extrema(U_interp)
anim = @animate for i in 1:nt_uniform
    plot(x, U_interp[:, i],
         ylim = (ymin, ymax),
         xlabel = "x",
         ylabel = "u(x, t)",
         title = "t = $(fixed_width_time(t_uniform[i]))",
         legend = false)
end

Deltat = t_uniform[2]-t_uniform[1];
mp4(anim, joinpath(@OUTPUT, "heat.mp4"), fps = 1/Deltat)