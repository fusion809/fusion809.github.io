~~~
<head>
<script src="/libs/common/generateTableXYZ.js"></script>
</head>
~~~

@def hassim=true;
@def title = "R&ouml;ssler equations solver"
@def params = (a=(val=0.1,), b=(val=0.1,), c=(val=14,), tf=(val=300,), x0=(val=-0.1,), y0=(val=0.5,), z0=(val=-0.6,), epsilon=(val=1e-8,))
@def type = "attractor"

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to [R&ouml;ssler equations](https://en.wikipedia.org/wiki/Rössler_attractor)

\begin{aligned}
    \dfrac{dx}{dt} &= - y - z \\
    \dfrac{dy}{dt} &= x + ay \\
    \dfrac{dz}{dt} &= b + z(x-c).
\end{aligned}

~~~
    {{ insert template.html}}
~~~