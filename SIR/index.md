~~~
<head>
<script src="/libs/common/generateTableSIR.js"></script>
</head>
~~~

@def hassim=true;
@def vars = ["S", "I", "R"]
@def title = "SIR equations solver"
@def params = (beta=(val=0.6, desc="A parameter that pertains to how many contacts there are per person and how easily the disease spreads from an infected person to an uninfected person."), gamma=(val=0.1, desc="A parameter that is a measure of how quickly people recover from the disease."), delta=(val=0.1, desc="A parameter with values from 0 to 1 pertaining to how effective quarantine measures are at slowing the disease outbreak. If \\(\\delta = 0\\), the measures are either non-existent or completely ineffective. If \\(\\delta = 1\\), all infected persons are immediately, as soon as they become infected, quarantined."), tf=(val=300, desc="End time of simulation."), S0=(val=89, desc="Initial number of susceptible individuals."), I0=(val=1, desc="Initial number of infected individuals."), R0=(val=0, desc="Initial number of recovered individuals."), epsilon=(val=1e-11,desc="Error tolerance for our numerical solution to the SIR equations."))
@def ids = ["tableOutputs", "phasePlotSIR", "phasePlotSI", "phasePlotSR", "phasePlotIR", "timePlot", "phasePlotSIR", "animation", "animation"]
@def funcs = ["generateTable()", "removeTable()", "generateSIRPhasePlot(solveProblem(RKF45, readInputs()))", "removeSIRPhasePlot()", "generateSIPhasePlot(solveProblem(RKF45, readInputs()))", "removeSIPhasePlot()", "generateSRPhasePlot(solveProblem(RKF45, readInputs()))", "removeSRPhasePlot()", "generateIRPhasePlot(solveProblem(RKF45, readInputs()))", "removeIRPhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))", "removeTimePlot()", "generatePlots(readInputs())", "removePlots()", "generateAnimation()", "removeAnimation()", "generateAllOutputs()", "removeAllOutputs()"]
@def labels = ["Tabulate the solution", "Remove the table", "Generate \\(S\\), \\(I\\) and \\(R\\) 3D phase plot", "Remove \\(S\\), \\(I\\) and \\(R\\) 3D phase plot", "Generate \\(S\\) and \\(I\\) phase plot", "Remove \\(S\\) and \\(I\\) phase plot", "Generate \\(S\\) and \\(R\\) phase plot", "Remove \\(S\\) and \\(R\\) phase plot", "Generate \\(I\\) and \\(R\\) phase plot", "Remove \\(I\\) and \\(R\\) phase plot", "Generate time plot for \\(S\\), \\(I\\) and \\(R\\)", "Remove time plot", "Generate all solution plots", "Remove all plots", "Generate an animation", "Remove animation","Generate all outputs","Remove all outputs"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the SIR equations with the $\delta$ parameter to account for quarantine effects:

\begin{aligned}
\dfrac{dS}{dt} &= -\dfrac{\beta I (1-\delta)S}{N} \\
\dfrac{dI}{dt} &= \dfrac{\beta I(1-\delta)S}{N} - \gamma I \\
\dfrac{dR}{dt} &= \gamma I.
\end{aligned}

Where $S$ is the number of susceptible persons, $I$ is the number of infected persons and $R$ is the number of recovered persons. $\beta$ is a parameter that pertains to the average number of contacts per person per time and the rate of transmission for the disease. $\gamma$ is the inverse of the average time a person is infected with the disease. Consequently, the basic reproduction number, which is typically represented as $R_0$ (not to be confused with the initial population of recovered individuals), is given by $\dfrac{\beta(1-\delta)}{\gamma}$ (Ridenhour et al., 2014). $N$ is the total population.

My original model had $\gamma I$ multiplied by $1-\delta$, but as quarantine should not affect how long it takes for people to recover, it should not affect this term.

If you would like to examine a more realistic model, with birth and death rates, loss of immunity, and incubation periods, then check out the [SEIRS model solver](/SEIRS/) webpage.

~~~
    {{ insert template.html}}
~~~

# Reference list

Ridenhour, B; Kowalik, JM; & Shay, DK (2014). [Unraveling $R_0$: Considerations for Public Health Applications](https://doi.org/10.2105/AJPH.2013.301704). *American Journal of Public Health*, 104(2): e32–e41. 