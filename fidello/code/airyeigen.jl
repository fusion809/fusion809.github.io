
using LinearAlgebra
using DataFrames
using SpecialFunctions, Plots

# Parameters
N = 80
L = 10.0
nev = 6

# Ascending Chebyshev-Lobatto grid
ξ = -cos.(π .* (0:N) ./ N)

# Chebyshev differentiation matrix
c = [2.0; ones(N-1); 2.0] .* (-1.0).^(0:N)

X = repeat(ξ, 1, N+1)
dX = X - X'

D = (c * (1 ./ c)') ./ (dX + I)
D -= Diagonal(vec(sum(D, dims=2)))

# Map to the semi-infinite interval
x = L .* (1 .+ ξ) ./ (1 .- ξ)
x[1] = 0.0
x[end] = Inf

# Differentiation matrices with respect to x
q = (1 .- ξ).^2 ./ (2L)

Dx = Diagonal(q) * D
D2 = Dx * Dx

# Construct the generalised eigenvalue problem
# u'' - xu = λu

A = zeros(N+1, N+1)
B = zeros(N+1, N+1)

# Neumann boundary condition: u'(0) = 0
A[1, :] = Dx[1, :]

# Differential equation at interior grid points
A[2:N, :] = D2[2:N, :]

A[N+1, :] = Dx[N+1, :]

A[2:N, 2:N] -= diagm(x[2:N])

# Eigenvalue terms at interior grid points
B[2:N, 2:N] = Matrix{Float64}(I, N-1, N-1)

# Solve the generalised eigenvalue problem
E = eigen(A, B)

# Retain finite, real, negative eigenvalues
idx = findall(k ->
    isfinite(real(E.values[k])) &&
    abs(imag(E.values[k])) < 1e-8 &&
    real(E.values[k]) < 0,
    eachindex(E.values)
)

# Sort eigenvalues from largest to smallest
idx = sort(idx, by=k -> real(E.values[k]), rev=true)
idx = idx[1:min(nev, length(idx))]

λ = real.(E.values[idx])

# Extract and normalise eigenfunctions
U = zeros(N+1, length(idx))

for (j, k) in enumerate(idx)
    U[:, j] = real.(E.vectors[:, k])
    U[:, j] ./= maximum(abs.(U[:, j]))
end

# Display eigenvalues
df = DataFrame(
    Eigenvalue = 1:length(λ),
    λ = λ
)

println(df)

# Interpolate first eigenfunction onto a linear grid [0, 5]
x_linear = collect(range(0.0, 10.0, length=500))

# Transform physical coordinates to Chebyshev coordinates
ξ_linear = (x_linear .- L) ./ (x_linear .+ L)

# Barycentric interpolation weights
w = 1.0 ./ c

# Barycentric interpolation
u_linear = zeros(length(x_linear))

for i in eachindex(x_linear)
    k = findfirst(j -> isapprox(ξ_linear[i], ξ[j];
                                atol=8eps(Float64), rtol=0.0),
                  eachindex(ξ))

    if !isnothing(k)
        u_linear[i] = U[k, 1]
    else
        terms = w ./ (ξ_linear[i] .- ξ)
        u_linear[i] = sum(terms .* U[:, 1]) / sum(terms)
    end
end

# Plot interpolated first eigenfunction
p = plot(
    x_linear,
    u_linear,
    xlabel = "x",
    ylabel = "u₁(x)",
    title = "First eigenfunction",
    label = "λ = $(round(λ[1], digits=10))",
    linewidth = 2,
    legend = :topright,
    size = (900, 600)
)

savefig(p, "airyeigen.svg")
