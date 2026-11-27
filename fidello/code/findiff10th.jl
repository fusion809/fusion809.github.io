
using LinearAlgebra
using SparseArrays
using DataFrames

# Generate finite-difference weights using polynomial moments.
function fdweights(offsets, m)
    s = length(offsets)

    V = BigFloat[
        BigFloat(offsets[j])^(k-1)
        for k in 1:s, j in 1:s
    ]

    rhs = zeros(BigFloat, s)
    rhs[m+1] = factorial(big(m))

    return Float64.(V \ rhs)
end

function solve_problem(N)
    h = 2π / N

    A = spzeros(N+1, N+1)
    b = zeros(N+1)

    # Dirichlet boundary conditions
    A[1, 1] = 1
    A[end, end] = 1

    b[1] = 1
    b[end] = 1

    # Interior finite-difference equations:
    # y'' + y = 0
    for i in 1:N-1

        # Select 12 grid points.
        # The stencil is shifted near the boundaries.
        start = clamp(i-5, 0, N-11)

        indices = collect(start:start+11)
        offsets = indices .- i

        weights = fdweights(offsets, 2)

        # Approximate h² y'' + h² y = 0
        A[i+1, start+1:start+12] .= weights
        A[i+1, i+1] += h^2
    end

    x = collect(LinRange(0, 2π, N+1))

    # Solve the discretised system
    y = A \ b

    return x, y, A
end

# Calculate a numerical solution
N = 20
x, y, A = solve_problem(N)

# Compare with cos(x), one possible exact solution
y_exact = cos.(x)

df = DataFrame(
    x = x,
    Approximation = y,
    Exact = y_exact,
    Error = abs.(y .- y_exact)
)

show(stdout, MIME("text/plain"), df)

println("\nMaximum difference from cos(x): ",
        maximum(abs.(y .- y_exact)))