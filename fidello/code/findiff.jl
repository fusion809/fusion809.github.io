using LinearAlgebra, DataFrames

N = 20
h = 2π / N

d  = vcat(1.0, fill(-2 + h^2, N - 1), 1.0)
dl = vcat(ones(N - 1), 0.0)
du = vcat(0.0, ones(N - 1))

A = Tridiagonal(dl, d, du)

b = zeros(N + 1)
b[1] = 1
b[end] = 1

x = LinRange(0, 2pi, N+1);
y = A \ b
yex = cos.(x);
df = DataFrame(x=x, Approximation=y, Exact=yex, Error=abs.(y.-yex))
show(stdout, MIME("text/plain"), df)
println()
