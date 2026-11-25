
using LinearAlgebra
using DataFrames
using Plots, QuadGK, LaTeXStrings

# Parameters
N = 40
gravity = 9.81
l = 1.0

α = gravity / l
# Integrand
f(θ) = 1 / sqrt(- sin(θ))

# Numerical integration
int, error = quadgk(f, -π, 0)

# Period
T = 2 * sqrt(l / (2gravity)) * int
tf = T;

# Newton iteration parameters
tol = 1e-12
maxiter = 50

# Chebyshev-Lobatto extrema grid (ascending order)
ξ = -cos.(π .* (0:N) ./ N)

# Map to the time interval [0, T]
t = (tf / 2) .* (1 .+ ξ)

# Chebyshev differentiation matrix
c = [2.0; ones(N-1); 2.0] .* (-1.0).^(0:N)

X = repeat(ξ, 1, N+1)
dX = X - X'

D = (c * (1 ./ c)') ./ (dX + LinearAlgebra.I)
D -= Diagonal(vec(sum(D, dims=2)))

# Differentiation matrices with respect to time
Dt = (2 / tf) .* D
D2t = Dt * Dt

# Initial approximation
θ = -(π / 2) .* (1 .- cos.(2π .* t ./ T) )

# Newton-Kantorovich iteration
converged = false

for n in 1:maxiter

    # Residual of the nonlinear differential equation
    F = D2t * θ + α .* cos.(θ)

    # Jacobian of the nonlinear system
    J = D2t - Diagonal(α .* sin.(θ))

    # Right-hand side of Newton correction equation
    rhs = -F

    # Initial condition: Δ(0) = -θ_n(0)
    J[1, :] .= 0.0
    J[1, 1] = 1.0
    rhs[1] = -θ[1]

    # Initial condition: Δ'(0) = -θ_n'(0)
    J[2, :] = Dt[1, :]
    rhs[2] = -(Dt * θ)[1]

    # Solve for the correction
    Δ = J \ rhs

    # Update approximation
    θ .+= Δ

    # Check convergence
    error = norm(Δ, Inf)

    println("Iteration $n: error = $error")

    if error < tol
        converged = true
        println("Newton-Kantorovich iteration converged.")
        break
    end
end

if !converged
    @warn "Newton-Kantorovich iteration did not converge."
end

# Calculate angular velocity and acceleration
θdot = Dt * θ
θddot = D2t * θ

# Display numerical results
df = DataFrame(
    t = t,
    θ = θ,
    θdot = θdot,
    θddot = θddot
)

println("\nNumerical solution:")
show(df, allrows=true, allcols=true)
println()

# Plot angular displacement
p1 = plot(
    t, θ,
    xlabel = L"t\ (\mathrm{seconds})",
    ylabel = L"\theta",
    legend=false,
    linewidth = 2
)

savefig(p1, "pendulum_theta.svg")