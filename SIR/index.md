~~~
<head>
<script src="/libs/common/generateTableXYZ.js"></script>
</head>
~~~

@def hassim=true;
@def title = "R&ouml;ssler equations solver"
@def params = (beta=(val=0.33, desc="A parameter that pertains to how many contacts there are per person and how easily the disease spreads from an infected person to an infected person."), gamma=(val=0.25, desc="A parameter that is a measure of how quickly people recover from the disease."), delta=(val=0.5, desc="A parameter with values from 0 to 1 pertaining to how effective quarantine measures are at slowing the disease outbreak. If \\(\\delta = 0\\), the measures are either non-existent or completely ineffective. If \\(\\delta = 1\\), all infected persons are immediately, as soon as they become infected, quarantined."), tf=(val=300,), S0=(val=89,), I0=(val=1,), R0=(val=0,), epsilon=(val=1e-11,))
@def type = "attractor"

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the SIR equations with the $\delta$ parameter to account for quarantine effects:

\begin{aligned}
\dfrac{dS}{dt} &= -\dfrac{\beta I (1-\delta)S}{N} \\
\dfrac{dI}{dt} &= \dfrac{\beta I(1-\delta)S}{N} - \gamma I \\
\dfrac{dR}{dt} &= \gamma I.
\end{aligned}

Where $S$ is the number of susceptible persons, $I$ is the number of infected persons and $R$ is the number of recovered persons. $\beta$ is a parameter that pertains to the average number of contacts per person per time and the rate of transmission for the disease. $\gamma$ is the inverse of the average time a person is infected with the disease. $N$ is the total population.

My original model had $\gamma I$ multiplied by $1-\delta$, but as quarantine should not affect how long it takes for people to recover, it should not affect this term.

~~~
    {{ insert template.html}}
~~~