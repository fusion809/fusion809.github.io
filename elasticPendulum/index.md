@def title="Elastic pendulum solver"
@def hassim=true
@def vars = ["z", "dz", "theta", "dtheta"]
@def ids = ["tableOutputs", "phasePlotZTheta","phasePlotDzZ","phasePlotDthetaTheta","timePlot","pendulumPlot", "pendulumTimePlot", "pendulumPlot", "phasePlotZTheta", "animation","animationThetaPhase","animationZPhase","animation","animation"]
@def params = (g=(val=9.81, desc="Acceleration due to gravity in metres (m) per second (s) squared (\\(\\mathrm{m}\\cdot \\mathrm{s}^{-2}\\))."),l0=(val=1, desc="Rest length (m) of pendulum rod."),k=(val=1, desc="Spring coefficient of pendulum."),m=(val=1, desc="Mass (kilograms or kg) of pendulum bob."),tf=(val=60, desc="End time (s) for the simulation."),z0=(val=1, desc="Initial value of \\(z\\) (m)."),dz0=(val=1, desc="Initial value of \\(\\dot{z}\\) (\\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\))."),theta0=(val=0, desc="Initial value of \\(\\theta\\) in radians (r)."),dtheta0=(val=0, desc="Initial value of \\(\\dot{\\theta}\\) (\\(\\mathrm{r}\\cdot \\mathrm{s}^{-1}\\))."),epsilon=(val=1e-9,))
@def funcs = ["fillTable(readInputs(), ['z', 'dz/dt', '&theta;', 'd&theta;/dt'])","removeTable()", "generateZThetaPhasePlot(solveProblem(RKF45, readInputs()))","removeZThetaPhasePlot()", "generateDzZPhasePlot(solveProblem(RKF45, readInputs()))","removeDzZPhasePlot()", "generateDthetaThetaPhasePlot(solveProblem(RKF45, readInputs()))","removeDthetaThetaPhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))","removeTimePlot()", "generatePendulumPlot(readInputs(), solveProblem(RKF45, readInputs()))", "removePendulumPlot()","generatePendulumTimePlot(readInputs(), solveProblem(RKF45, readInputs()))", "removePendulumTimePlot()", "generatePendulumPlots(readInputs(), solveProblem(RKF45, readInputs()))", "removePendulumPlots()", "generatePlots(readInputs())","removePlots()", "generateAnimation()","removeAnimation()", "generateThetaPhaseAnimation()","removeThetaPhaseAnimation()", "generateZPhaseAnimation()","removeZPhaseAnimation()", "generateAnimations()","removeAnimations()", "generateAllOutputs()", "removeAllOutputs()"]
@def labels = ["Tabulate the solution","Remove the solution table","Generate a \\(\\theta\\) against \\(z\\) phase plot","Remove \\(\\theta\\) against \\(z\\) phase plot","Generate a \\(\\dot{z}\\) against \\(z\\) phase plot","Remove \\(\\dot{z}\\) against \\(z\\) phase plot","Generate a \\(\\dot{\\theta}\\) against \\(\\theta\\) phase plot","Remove \\(\\dot{\\theta}\\) against \\(\\theta\\) phase plot","Generate a \\(z\\), \\(\\dot{z}\\), \\(\\theta\\) and \\(\\dot{\\theta}\\) against time plot","Remove \\(z\\), \\(\\dot{z}\\), \\(\\theta\\) and \\(\\dot{\\theta}\\) against time plot", "Generate pendulum position plot", "Remove pendulum position plot", "Generate pendulum position against time plot", "Remove pendulum position against time plot", "Generate pendulum plots", "Remove pendulum plots", "Generate all solution plots","Remove all plots","Generate an animation","Remove animation","Generate a \\(\\theta\\) phase plot animation","Remove \\(\\theta\\) phase plot animation","Generate a \\(z\\) phase plot animation","Remove \\(z\\) phase plot animation","Generate all animations","Remove all animations", "Generate all outputs", "Remove all outputs"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of the [elastic pendulum](https://en.wikipedia.org/wiki/Elastic_pendulum)

\begin{aligned}
    \dfrac{d^{2}z}{dt^2} &= (l_0 + z) \dot{\theta}^2 - \dfrac{kz}{m} + g \sin{\theta} \\
    \dfrac{d^{2} \theta}{dt^2} &= -\dfrac{g}{l_0 + z} \cos{\theta} - \dfrac{2\dot{z}\dot{\theta}}{l_0 + z}.
\end{aligned}

Where:
            
* $z$ is the extension of the pendulum beyond its rest length (in metres).
* $\theta$ is the angle of the pendulum relative to the positive $x$-axis.
* $g$ is the acceleration due to gravity (in metres per second squared).
* $l_0$ is the rest length of the pendulum.

~~~
{{ insert template.html }}
~~~