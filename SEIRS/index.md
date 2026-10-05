~~~
<head>
<script src="/libs/common/generateTableSIR.js"></script>
</head>
~~~

@def hassim=true;
@def title = "SEIRS solver"
@def params = (alpha=(val=0.01, desc="Death by infection rate."), beta=(val=2, desc="Disease transmission rate."), gamma=(val=0.1, desc="Reciprocal of recovery period."), delta=(val=0.1, desc="Quarantine efficacy rate."), lambda=(val=1e-4, desc="Birth rate."), mu=(val=1e-5, desc="Infection-independent death rate."), sigma=(val=0.33, desc="Inverse of incubation period."), omega=(val=1e-4, desc="Rate of immunity loss in recovered individuals."), tf=(val=300, desc="End time of simulation."), S0=(val=89, desc="Initial number of susceptible individuals."), E0=(val=0, desc="Initial number of exposed individuals."), I0=(val=1, desc="Initial number of infected individuals."), R0=(val=0, desc="Initial number of recovered individuals."), epsilon=(val=1e-11,desc="Error tolerance for our numerical solution to the SIR equations."))
@def ids = ["tableOutputs", "phasePlotSIR", "phasePlotSEI", "phasePlotSER", "phasePlotEIR", "phasePlotSI", "phasePlotSR", "phasePlotIR", "timePlot", "phasePlotSIR", "animationSIR", "animationSEI", "animationSER", "animationEIR", "animationSIR", "animationSIR"];
@def vars=["S", "E", "I", "R"]
@def funcs = ["fillTable(readInputs(), ['Susceptible', 'Exposed', 'Infectious', 'Recovered'])","removeTable()","generateSIRPhasePlot(solveProblem(RKF45, readInputs()))","removeSIRPhasePlot()","generateSEIPhasePlot(solveProblem(RKF45, readInputs()))","removeSEIPhasePlot()","generateSERPhasePlot(solveProblem(RKF45, readInputs()))","removeSERPhasePlot()","generateEIRPhasePlot(solveProblem(RKF45, readInputs()))","removeEIRPhasePlot()","generateSIPhasePlot(solveProblem(RKF45, readInputs()))","removeSIPhasePlot()","generateSRPhasePlot(solveProblem(RKF45, readInputs()))","removeSRPhasePlot()","generateIRPhasePlot(solveProblem(RKF45, readInputs()))","removeIRPhasePlot()","generateTimePlot(solveProblem(RKF45, readInputs()))","removeTimePlot()","generatePlots(readInputs())","removePlots()","generateAnimationSIR()","removeAnimationSIR()","generateAnimationSEI()","removeAnimationSEI()","generateAnimationSER()","removeAnimationSER()","generateAnimationEIR()","removeAnimationEIR()","generateAnimations()","removeAnimations()","generateAllOutputs()", "removeAllOutputs()"] 
@def labels = ["Tabulate the solution","Remove the solution table","Generate a \\(S\\), \\(I\\) and \\(R\\) phase plot","Remove \\(S\\), \\(I\\) and \\(R\\) plot","Generate a \\(S\\), \\(E\\) and \\(I\\) phase plot","Remove \\(S\\), \\(E\\) and \\(I\\) plot","Generate a \\(S\\), \\(E\\) and \\(R\\) phase plot","Remove \\(S\\), \\(E\\) and \\(R\\) plot","Generate a \\(E\\), \\(I\\) and \\(R\\) phase plot","Remove \\(E\\), \\(I\\) and \\(R\\) plot","Generate a \\(S\\) and \\(I\\) phase plot","Remove \\(S\\) and \\(I\\) plot","Generate a \\(S\\) and \\(R\\) phase plot","Remove \\(S\\) and \\(R\\) plot","Generate an \\(I\\) and \\(R\\) phase plot","Remove \\(I\\) and \\(R\\) plot","Generate a \\(S\\), \\(E\\), \\(I\\) and \\(R\\) against time plot","Remove time plot","Generate all solution plots","Remove all plots","Generate a \\(S\\), \\(I\\) and \\(R\\) phase plot animation","Remove \\(S\\), \\(I\\) and \\(R\\) animation","Generate a \\(S\\), \\(E\\) and \\(I\\) phase plot animation","Remove \\(S\\), \\(E\\) and \\(I\\) animation","Generate a \\(S\\), \\(E\\) and \\(R\\) phase plot animation","Remove \\(S\\), \\(E\\) and \\(R\\) animation","Generate \\(E\\), \\(I\\) and \\(R\\) phase plot animation","Remove a \\(E\\), \\(I\\) and \\(R\\) phase plot animation","Generate all animations","Remove all animations","Generate all outputs","Remove all outputs"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking RKF45](/RKF45/) to approximate the solution to the SEIRS equations:

\begin{aligned}
\frac{dS}{dt} & = \lambda N + \omega R - \mu S - \frac{\beta I (1-\delta)S}{N} \\
\frac{dE}{dt} & = \frac{\beta I (1-\delta) S}{N} - (\mu + \sigma ) E \\
\frac{dI}{dt} & = \sigma E - (\alpha + \gamma +\mu ) I \\
\frac{dR}{dt} & = \gamma I  - (\mu + \omega) R.
\end{aligned}

Where:
* $t$ is time and is typically taken to be in units of days. 
* $S$ is the number of susceptible persons.
* $E$ is the number of exposed persons, that is those that have been exposed to the disease but have not yet become infectious.
* $I$ is the number of infectious persons.
* $R$ is the number of recovered persons.
* $\alpha$ is the infection mortality rate. This does not perfectly align with how we usually use the term "mortality rate" as it is not just the rate of death per case of infection it is also per time. So more rapidly lethal infections would have a higher $\alpha$ than their slower counterparts. 
* $\beta$ is a parameter that pertains to the average number of contacts per person per time and the rate of transmission for the disease. A higher $\beta$ corresponds to a more contagious disease as it is directly correlated with $R_0$, the basic reproduction number.
* $\gamma$ is the reciprocal of the average time a person is infected with the disease. Diseases with faster recovery therefore have larger $\gamma$ and those with longer recovery periods have lower $\gamma$.
* $\delta$ refers to the efficacy of quarantine. $\delta=0$ means quarantine measures are completely ineffective. $\delta=1$ means they are completely effective.  
* $\lambda$ is the birth rate. It is essentially how many new people join the population per time per existing member of the population. Each new person is assumed to be susceptible, hence this may have limited efficacy in modelling immigration if the immigrants could potentially have some immunity to the disease. For instance, if the disease is a vaccine-preventable disease, immigrants may have been vaccinated and hence would be best place in the "recovered" category for the purpose of the model. 
* $\mu$ is the rate of death in the population from causes other than the infection being modelled. It is essentially how many people die per population per time. 
* $\sigma$ is the latency rate, which is the reciprocal of the incubation period. Smaller incubation periods therefore lead to larger $\sigma$. 
* $\omega$ is the rate of immunity loss in recovered individuals. Like $\lambda$ and $\mu$ it is per population per time, although the population it applies to is the recovered population.  
* $N=S+E+I+R$ is the total population.

This model is heavily based on [Bjørnstad et al. (2020)](#reference_list), with some amendments. Specifically, $\delta$ has been added to account for quarantine effects and $\lambda$ has been accounted for a innate birth rate that may not match the innate death rate of $\mu$. Hence the basic reproduction number $R_0$ (not to be confused with the initial number of recovered individuals) is given by:

\begin{aligned}
R_0 = \dfrac{\beta (1-\delta)\sigma}{(\mu+\sigma)(\alpha + \gamma + \mu)}.
\end{aligned}

If you would like to examine a more simple model, without birth and death rates, immunity loss, and incubation periods, then check out the [SIR model solver](/SIR/) webpage.

~~~
    {{ insert template.html}}
~~~

# Reference list

Bjørnstad ON, Shea K, Krzywinski M, & Altman N (2020). [The SEIRS model for infectious disease dynamics](https://doi.org/10.1038/s41592-020-0856-2). *Nature Methods*, 17(6): 557–558. doi: [10.1038/s41592-020-0856-2](https://doi.org/10.1038/s41592-020-0856-2). PMID [32499633](https://pubmed.ncbi.nlm.nih.gov/32499633/).
