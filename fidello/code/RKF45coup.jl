using DataFrames

a = 0.0;
b = 2pi;
N = 5;
h = (b-a)/N;

epsilon = 1e-8;

# Initial conditions
x = [a];
y = [0.0];
u = [1.0];

while x[end] < b
    xn = x[end];
    yn = y[end];
    un = u[end];

    # Prevent overshooting the endpoint
    hs = min(h, b-xn);

    # Fehlberg coefficients
    k1 = hs*un;
    l1 = -hs*yn;

    k2 = hs*(un + l1/4);
    l2 = -hs*(yn + k1/4);

    k3 = hs*(un + 3l1/32 + 9l2/32);
    l3 = -hs*(yn + 3k1/32 + 9k2/32);

    k4 = hs*(un + 1932l1/2197 - 7200l2/2197 +
             7296l3/2197);
    l4 = -hs*(yn + 1932k1/2197 - 7200k2/2197 +
              7296k3/2197);

    k5 = hs*(un + 439l1/216 - 8l2 +
             3680l3/513 - 845l4/4104);
    l5 = -hs*(yn + 439k1/216 - 8k2 +
              3680k3/513 - 845k4/4104);

    k6 = hs*(un - 8l1/27 + 2l2 - 3544l3/2565 +
             1859l4/4104 - 11l5/40);
    l6 = -hs*(yn - 8k1/27 + 2k2 - 3544k3/2565 +
              1859k4/4104 - 11k5/40);

    # Fourth-order approximations
    y4 = yn + 25k1/216 + 1408k3/2565 +
         2197k4/4104 - k5/5;

    u4 = un + 25l1/216 + 1408l3/2565 +
         2197l4/4104 - l5/5;

    # Fifth-order approximations
    y5 = yn + 16k1/135 + 6656k3/12825 +
         28561k4/56430 - 9k5/50 + 2k6/55;

    u5 = un + 16l1/135 + 6656l3/12825 +
         28561l4/56430 - 9l5/50 + 2l6/55;

    # Truncation error estimate
    TE = max(abs(y5-y4), abs(u5-u4));

    # Step-size adjustment
    s = TE == 0 ? 5.0 : (epsilon/(2TE))^(1/5);

    # Accept or reject the step
    if TE < epsilon
        xn_new = xn + hs;

        if hs == b-xn
            xn_new = b;
        end

        push!(x, xn_new);
        push!(y, y5);
        push!(u, u5);
    end

    # Update step size
    global h = hs*s;

    if h <= eps(Float64)*max(1.0, abs(xn))
        error("Step size became too small.")
    end
end

# Exact solutions
y_exact = sin.(x);
u_exact = cos.(x);

# Absolute errors
error_y = abs.(y-y_exact);
error_u = abs.(u-u_exact);

# Construct table
df = DataFrame(
    :x => x,
    :y => y,
    Symbol("y'") => u,
    :y_exact => y_exact,
    Symbol("y'_exact") => u_exact,
    :error_y => error_y,
    Symbol("error_y'") => error_u
);

show(stdout, MIME("text/plain"), df)
println()
