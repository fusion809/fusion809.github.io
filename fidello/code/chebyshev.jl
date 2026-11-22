using LinearAlgebra, SpecialFunctions, DataFrames

# Parameters
n = 20
L = 6.0

# Chebyshev extrema, ordered from -1 to 1
ξ = -cos.(π .* (0:n) ./ n)

# Transformation to the semi-infinite domain
x = [i == n+1 ? Inf : L*(1 + ξ[i])/(1 - ξ[i])
     for i in 1:n+1]

# Derivative transformation
q = (1 .- ξ).^2 ./ (2L)

# Second derivative matrix
D2 = Diagonal(q) * D * Diagonal(q) * D

# Collocation matrix
A = copy(D2)

for i in 2:n
    A[i,i] -= x[i]
end

# Boundary conditions: y(0) = Ai(0), y(∞) = 0
A[1,:] .= 0
A[1,1] = 1

A[end,:] .= 0
A[end,end] = 1

# Right-hand side
b = zeros(n+1)
b[1] = airyai(0)

# Solve
u = A \ b

# Chebyshev coefficients
V = [cos((j-1)*acos(ξ[i])) for i in 1:n+1, j in 1:n+1]
a = V \ u

# Exact solution
yexact = [i == n+1 ? 0.0 : airyai(x[i])
          for i in 1:n+1]

# Results
df = DataFrame(
    ξ = ξ,
    x = x,
    Approximation = u,
    Exact = yexact,
    Error = abs.(u .- yexact)
)

show(stdout, MIME("text/plain"), df)
println()
