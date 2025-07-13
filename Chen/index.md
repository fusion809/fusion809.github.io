~~~
<head>
<script src="/libs/common/generateTableXYZ.js"></script>
</head>
~~~

@def hassim=true;
@def title = "Chen solver"
@def params = (a = (val=40, desc="Problem parameter"), b=(val=3, desc="Problem parameter"), c=(val=28, desc="Problem parameter"), tf=(val=120, desc="End time for the simulation in seconds (s)"), x0=(val=-0.1, desc="Initial \\(x\\) coordinate"), y0=(val=0.5, desc="Initial \\(y\\) coordinate"), z0=(val=-0.6, desc="Initial \\(z\\) coordinate"), epsilon=(val=1e-8,))
@def type = "attractor"

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the [Chen equations](https://en.wikipedia.org/wiki/Multiscroll_attractor)

\begin{aligned}
    \dfrac{dx}{dt} &= a (y-x) \\
    \dfrac{dy}{dt} &= x(c-a-z) + cy \\
    \dfrac{dz}{dt} &= xy - b z.
\end{aligned}

~~~
    {{ insert template.html}}
~~~