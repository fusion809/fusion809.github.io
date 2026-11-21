# This file was generated, do not modify it. # hide
using DataFrames

a = 0.0;
b = 1.0;
N = 5;
h = (b-a)/N;

atol = 1e-5;
max_factor = 5;

function f(x, y)
    return y;
end

x = [a];
yRKF45 = [1.0];
error = [0.0];

while x[end] < b
    xi = x[end];
    yi = yRKF45[end];

    # Prevent overshooting the endpoint
    hs = min(h, b-xi);

    # Runge-Kutta-Fehlberg coefficients
    k1 = hs*f(xi, yi);

    k2 = hs*f(xi + hs/4,
        yi + k1/4);

    k3 = hs*f(xi + 3hs/8,
        yi + 3k1/32 + 9k2/32);

    k4 = hs*f(xi + 12hs/13,
        yi + 1932k1/2197 - 7200k2/2197 + 7296k3/2197);

    k5 = hs*f(xi + hs,
        yi + 439k1/216 - 8k2 + 3680k3/513 - 845k4/4104);

    k6 = hs*f(xi + hs/2,
        yi - 8k1/27 + 2k2 - 3544k3/2565 +
        1859k4/4104 - 11k5/40);

    # Fourth-order approximation
    y4 = yi + 25k1/216 + 1408k3/2565 +
         2197k4/4104 - k5/5;

    # Fifth-order approximation
    y5 = yi + 16k1/135 + 6656k3/12825 +
         28561k4/56430 - 9k5/50 + 2k6/55;

    # Error estimation
    err = abs(y5-y4);

    # Accept or reject the step
    if err < atol
        push!(x, xi+hs)
        push!(yRKF45, y5)
        push!(error, err)
    end

    # Adjust the step size
    if err == 0
        s = max_factor;
    else
        s = (atol/(2*err))^(1/5);
    end

    global h = h*s
    if h < eps(max(abs(xi), 1.0))
        Base.error("Step size became too small.")
    end
end

df = DataFrame(
    x = x,
    RKF45 = yRKF45,
    Exact = exp.(x),
    Error = error
);

show(stdout, MIME("text/plain"), df)
println()