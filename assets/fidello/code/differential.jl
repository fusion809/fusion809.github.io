# This file was generated, do not modify it. # hide
using DataFrames;
N = 5;
a = 0;
b = 1;
h = (b-a)/N;
x = LinRange(a, b, N+1);
yEul = ones(N+1);
yRK4 = ones(N+1);
function f(x, y)
    return y;
end

for i in 1:N
    yEul[i+1] = yEul[i] + h*f(x[i], yEul[i]);
    k1 = h*f(x[i], yRK4[i]);
    k2 = h*f(x[i]+h/2, yRK4[i]+k1/2);
    k3 = h*f(x[i]+h/2, yRK4[i]+k2/2);
    k4 = h*f(x[i]+h, yRK4[i]+k3);
    yRK4[i+1] = yRK4[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
end
df = DataFrame(x = x, Euler = yEul, RK4 = yRK4, Exact = exp.(x));
show(stdout, MIME("text/plain"), df)
println()