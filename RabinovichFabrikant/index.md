~~~
<head>
<script src="/libs/common/generateTableXYZ.js"></script>
</head>
~~~

@def hassim=true;
@def title = "Rabinovich&ndash;Fabrikant equations solver"
@def params = (alpha = (val=1.1,), gamma=(val=0.87,), tf=(val=60,), x0=(val=-1,), y0=(val=0,), z0=(val=0.5,), epsilon=(val=1e-8,))
@def type = "attractor"

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to [Rabinovich&ndash;Fabrikant equations](https://en.wikipedia.org/wiki/Lorenz_system)

\begin{aligned}
    \dfrac{dx}{dt} &= y(z-1+x^2) + \gamma x \\
    \dfrac{dy}{dt} &= x (3z+1-x^2) + \gamma y \\
    \dfrac{dz}{dt} &= -2z(\alpha+xy).
\end{aligned}

~~~
    {{ insert template.html}}
~~~