using LinearAlgebra, SpecialFunctions, DataFrames

# Parameters
n = 40
L = 10.0

# Chebyshev extrema
ξ = -cos.(π .* (0:n) ./ n)

# Chebyshev differentiation matrix
c = [2.0; ones(n-1); 2.0] .* (-1.0).^(0:n)

X = repeat(ξ, 1, n+1)
dX = X - X'

D = (c * (1 ./ c)') ./ (dX + I)
D -= Diagonal(vec(sum(D, dims=2)))

# Transformation to the semi-infinite domain
x = L*(ξ.+1)./(-ξ.+1)

# Derivative transformation
q = (1 .- ξ).^2 ./ (2L)

# Second derivative matrix
D2 = Diagonal(q) * D * Diagonal(q) * D

# Collocation matrix: y'' - xy = 0
A = copy(D2)-diagm(x);

# Boundary conditions
A[1,:] .= 0
A[1,1] = 1

A[end,:] .= 0
A[end,end] = 1

# Right-hand side: y(0) = Ai(0), y(∞) = 0
b = zeros(n+1)
b[1] = airyai(0)

# Solve
u = A \ b

# Chebyshev coefficients
V = [cos((j-1)*acos(ξ[i])) for i in 1:n+1, j in 1:n+1]
a = V \ u

# Exact solution
uexact = [i == n+1 || x[i] > 100 ? 0.0 :
          airyai(x[i])
          for i in 1:n+1]

# Results
df = DataFrame(
    y = ξ,
    x = x,
    Approximation = u,
    Exact = yexact,
    Error = abs.(u .- uexact)
)

show(stdout, MIME("text/plain"), df)
println()

println("\nMaximum absolute error: ", maximum(abs.(u .- yexact)))
