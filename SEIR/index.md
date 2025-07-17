~~~
<head>
<script src="/libs/common/generateTableSIR.js"></script>
</head>
~~~

@def hassim=true;
@def title = "SEIR solver"
@def params = (a=(val=0.33, desc="Inverse of the average incubation period."), beta=(val=0.33, desc="A parameter that pertains to how many contacts there are per person and how easily the disease spreads from an infected person to an infected person."), gamma=(val=0.25, desc="A parameter that is a measure of how quickly people recover from the disease."), delta=(val=0.5, desc="A parameter with values from 0 to 1 pertaining to how effective quarantine measures are at slowing the disease outbreak. If \\(\\delta = 0\\), the measures are either non-existent or completely ineffective. If \\(\\delta = 1\\), all infected persons are immediately, as soon as they become infected, quarantined."), lambda=(val=1e-4, desc="Birth rate."), mu=(val=1e-5, desc="Death rate."), tf=(val=300,), S0=(val=89,), E0=(val=0,), I0=(val=11,), R0=(val=0,), epsilon=(val=1e-11,))
@def ids = ["tableOutputs", "phasePlotSIR", "phasePlotSEI", "phasePlotSER", "phasePlotEIR", "phasePlotSI", "phasePlotSR", "phasePlotIR", "timePlot", "phasePlotSIR", "animationSIR", "animationSEI", "animationSER", "animationEIR", "animation", "animation"];
@def funcs = ["fillTable(readInputs(), ['Susceptible', 'Exposed', 'Infectious', 'Recovered'])","removeTable()","generateSIRPhasePlot(solveProblem(RKF45, readInputs()))","removeSIRPhasePlot()","generateSEIPhasePlot(solveProblem(RKF45, readInputs()))","removeSEIPhasePlot()","generateSERPhasePlot(solveProblem(RKF45, readInputs()))","removeSERPhasePlot()","generateEIRPhasePlot(solveProblem(RKF45, readInputs()))","removeEIRPhasePlot()","generateSIPhasePlot(solveProblem(RKF45, readInputs()))","removeSIPhasePlot()","generateSRPhasePlot(solveProblem(RKF45, readInputs()))","removeSRPhasePlot()","generateIRPhasePlot(solveProblem(RKF45, readInputs()))","removeIRPhasePlot()","generateTimePlot(solveProblem(RKF45, readInputs()))","removeTimePlot()","generatePlots(readInputs())","removePlots()","generateAnimationSIR()","removeAnimationSIR()","generateAnimationSEI()","removeAnimationSEI()","generateAnimationSER()","removeAnimationSER()","generateAnimationEIR()","removeAnimationEIR()","generateAnimations()","removeAnimations()","generateAllOutputs()", "removeAllOutputs()"] 
@def labels = ["Tabulate the solution","Remove the solution table","Generate a \\(S\\), \\(I\\) and \\(R\\) phase plot","Remove \\(S\\), \\(I\\) and \\(R\\) plot","Generate a \\(S\\), \\(E\\) and \\(I\\) phase plot","Remove \\(S\\), \\(E\\) and \\(I\\) plot","Generate a \\(S\\), \\(E\\) and \\(R\\) phase plot","Remove \\(S\\), \\(E\\) and \\(R\\) plot","Generate a \\(E\\), \\(I\\) and \\(R\\) phase plot","Remove \\(E\\), \\(I\\) and \\(R\\) plot","Generate a \\(S\\) and \\(I\\) phase plot","Remove \\(S\\) and \\(I\\) plot","Generate a \\(S\\) and \\(R\\) phase plot","Remove \\(S\\) and \\(R\\) plot","Generate an \\(I\\) and \\(R\\) phase plot","Remove \\(I\\) and \\(R\\) plot","Generate a \\(S\\), \\(E\\), \\(I\\) and \\(R\\) against time plot","Remove time plot","Generate all solution plots","Remove all plots","Generate a \\(S\\), \\(I\\) and \\(R\\) phase plot animation","Remove \\(S\\), \\(I\\) and \\(R\\) animation","Generate a \\(S\\), \\(E\\) and \\(I\\) phase plot animation","Remove \\(S\\), \\(E\\) and \\(I\\) animation","Generate a \\(S\\), \\(E\\) and \\(R\\) phase plot animation","Remove \\(S\\), \\(E\\) and \\(R\\) animation","Generate \\(E\\), \\(I\\) and \\(R\\) phase plot animation","Remove a \\(E\\), \\(I\\) and \\(R\\) phase plot animation","Generate all animations","Remove all animations","Generate all outputs","Remove all outputs"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking RKF45](/RKF45/) to approximate the solution to the SEIR equations with the $\delta$ parameter to account for quarantine effects:

\begin{aligned}
\frac{dS}{dt} & = \Lambda N - \mu S - \frac{\beta I (1-\delta)S}{N} \\
\frac{dE}{dt} & = \frac{\beta I (1-\delta) S}{N} - (\mu +a ) E \\
\frac{dI}{dt} & = a E - (\gamma +\mu ) I \\
\frac{dR}{dt} & = \gamma I  - \mu R.
\end{aligned}

Where $S$ is the number of susceptible persons, $E$ is the number of exposed persons, $I$ is the number of infectious persons and $R$ is the number of recovered persons. $a$ is the inverse of the average incubation period. $\beta$ is a parameter that pertains to the average number of contacts per person per time and the rate of transmission for the disease. $\gamma$ is the inverse of the average time a person is infected with the disease. $\Lambda$ is the birth rate. $\mu$ is the overall population death rate (not only including the disease death rate). $N$ is the total population.

My original model had $\gamma I$ multiplied by $1-\delta$, but as quarantine should not affect how long it takes for people to recover, it should not affect this term.

~~~
    {{ insert template.html}}
~~~