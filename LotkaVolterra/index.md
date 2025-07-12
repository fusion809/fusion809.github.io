@def title="Lotka-Volterra equations solver"
@def hassim=true
@def params = (alpha=(val=1.1, desc="Natural growth rate of the population of prey animals."),beta=(val=0.4, desc="The rate at which prey animals are killed by the predators."),gamma=(val=0.4, desc="The rate at which predator animals die in the absence of their prey."),delta=(val=0.1, desc="The rate at which the population of predator animals increases due to the presence of their prey."),t0=(val=0, desc="Starting time for the simulation in seconds (s)."),tf=(val=37, desc="End time for the simulation in seconds."),x0=(val=100, desc="Prey population."),y0=(val=10, desc="Predator population."), epsilon=(val=1e-9,), t1=(val=30,))

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to Lotka-Volterra equations:

\begin{aligned}
\dfrac{dx}{dt} &= \alpha x - \beta xy \\
\dfrac{dy}{dt} &= \delta xy - \gamma y
\end{aligned}

where $x$ is the number of prey animals and $y$ is the number of predator animals and $\alpha, \beta, \gamma$, and $\delta$ describe their interactions with one another.

~~~
{{ insert parameter_form.html }}
<!--Buttons-->
{{ insert 2D_buttons.html }}

<!--Where the table and plot goes-->
{{ insert 2D_output.html }}
{{ insert page_foot_general.html }}
~~~