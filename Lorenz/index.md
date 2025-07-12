~~~
<head>
<script src="/libs/common/generateTableXYZ.js"></script>
</head>
~~~

@def hassim=true;
@def title = "Lorenz solver"
@def params = (sigma = (val=10,), rho=(val=28,), beta=(val=2.66666666666667,), tf=(val=60,), x0=(val=1,), y0=(val=1,), z0=(val=1,), epsilon=(val=1e-8,))
@def type = "attractor"

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to [Lorenz equations](https://en.wikipedia.org/wiki/Lorenz_system)

\begin{aligned}
    \dfrac{dx}{dt} &= \sigma (y-x) \\
    \dfrac{dy}{dt} &= x (\rho - z) - y \\
    \dfrac{dz}{dt} &= xy - \beta z.
\end{aligned}

~~~
    {{ insert template.html}}
~~~