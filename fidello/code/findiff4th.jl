
using LinearAlgebra, DataFrames, SparseArrays

N = 20
h = 2π / N

# Construct sparse coefficient matrix
A = spzeros(N + 1, N + 1)

# Boundary conditions
A[1, 1] = 1
A[end, end] = 1

# First interior row: fourth-order forward stencil
A[2, 1:6] .= [10, -15 + 12h^2, -4, 14, -6, 1]

# Central five-point stencil
for i in 2:N-2
    j = i + 1
    A[j, j-2:j+2] .= [-1, 16, -(30 - 12h^2), 16, -1]
end

# Last interior row: fourth-order backward stencil
A[N, N-4:N+1] .= [1, -6, 14, -4, -15 + 12h^2, 10]

# Right-hand side
b = zeros(N + 1)
b[1] = 1
b[end] = 1

# Solve
x = LinRange(0, 2π, N+1)
y = A \ b

# Exact solution and errors
yex = cos.(x)

df = DataFrame(
    x = x,
    Approximation = y,
    Exact = yex,
    Error = abs.(y .- yex)
)

show(stdout, MIME("text/plain"), df)
println()
