@def title="Elastic pendulum solver"
@def hassim=true
@def vars = ["r", "dr", "theta", "dtheta"]
@def ids = ["tableOutputs", "phasePlotRTheta","phasePlotDrR","phasePlotDthetaTheta","timePlot","pendulumPlot", "pendulumTimePlot", "pendulumPlot", "phasePlotRTheta", "animation","animationThetaPhase","animationRPhase","animation","animation"]
@def params = (g=(val=9.81, desc="Acceleration due to gravity in metres (m) per second (s) squared (\\(\\mathrm{m}\\cdot \\mathrm{s}^{-2}\\))."),l0=(val=1, desc="Rest length (m) of pendulum rod."),k=(val=1, desc="Spring coefficient of pendulum."),mb=(val=1, desc="Mass (kilograms or kg) of pendulum bob."),mr=(val=1,desc="Mass (kilograms or kg) of pendulum rod."), bb=(val=0, desc="Linear dissipation coefficient for the pendulum bob."), br=(val=0, desc="Linear dissipation coefficient for the pendulum rod."), cb=(val=0, desc="Quadratic dissipation coefficient for the pendulum bob."), cr=(val=0, desc="Quadratic dissipation coefficient for the pendulum rod."), tf=(val=60, desc="End time (s) for the simulation."), r0=(val=1, desc="Initial value of \\(r\\) (m)."), dr0=(val=1, desc="Initial value of \\(\\dot{r}\\) (\\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\))."),theta0=(val=0, desc="Initial value of \\(\\theta\\) in radians (r)."),dtheta0=(val=0, desc="Initial value of \\(\\dot{\\theta}\\) (\\(\\mathrm{r}\\cdot \\mathrm{s}^{-1}\\))."),epsilon=(val=1e-9,))
@def funcs = ["fillTable(readInputs(), ['r', 'dr/dt', '&theta;', 'd&theta;/dt'])","removeTable()", "generateRThetaPhasePlot(solveProblem(RKF45, readInputs()))","removeRThetaPhasePlot()", "generateDrRPhasePlot(solveProblem(RKF45, readInputs()))","removeDrRPhasePlot()", "generateDthetaThetaPhasePlot(solveProblem(RKF45, readInputs()))","removeDthetaThetaPhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))","removeTimePlot()", "generatePendulumPlot(readInputs(), solveProblem(RKF45, readInputs()))", "removePendulumPlot()","generatePendulumTimePlot(readInputs(), solveProblem(RKF45, readInputs()))", "removePendulumTimePlot()", "generatePendulumPlots(readInputs(), solveProblem(RKF45, readInputs()))", "removePendulumPlots()", "generatePlots(readInputs())","removePlots()", "generateAnimation()","removeAnimation()", "generateThetaPhaseAnimation()","removeThetaPhaseAnimation()", "generateRPhaseAnimation()","removeRPhaseAnimation()", "generateAnimations()","removeAnimations()", "generateAllOutputs()", "removeAllOutputs()"]
@def labels = ["Tabulate the solution","Remove the solution table","Generate a \\(\\theta\\) against \\(r\\) phase plot","Remove \\(\\theta\\) against \\(r\\) phase plot","Generate a \\(\\dot{r}\\) against \\(r\\) phase plot","Remove \\(\\dot{r}\\) against \\(r\\) phase plot","Generate a \\(\\dot{\\theta}\\) against \\(\\theta\\) phase plot","Remove \\(\\dot{\\theta}\\) against \\(\\theta\\) phase plot","Generate a \\(r\\), \\(\\dot{r}\\), \\(\\theta\\) and \\(\\dot{\\theta}\\) against time plot","Remove \\(r\\), \\(\\dot{r}\\), \\(\\theta\\) and \\(\\dot{\\theta}\\) against time plot", "Generate pendulum position plot", "Remove pendulum position plot", "Generate pendulum position against time plot", "Remove pendulum position against time plot", "Generate pendulum plots", "Remove pendulum plots", "Generate all solution plots","Remove all plots","Generate an animation","Remove animation","Generate a \\(\\theta\\) phase plot animation","Remove \\(\\theta\\) phase plot animation","Generate a \\(r\\) phase plot animation","Remove \\(r\\) phase plot animation","Generate all animations","Remove all animations", "Generate all outputs", "Remove all outputs"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of the [elastic pendulum](https://en.wikipedia.org/wiki/Elastic_pendulum)

\begin{aligned}
    \dfrac{d^{2}r}{dt^2} &= r \dot{\theta}^2 - \dfrac{\mu_1}{M_1} g\sin{\theta} + \dfrac{Q_r-k(r-l_0)}{M_1} \\
    \dfrac{d^{2} \theta}{dt^2} &= -\dfrac{2\dot{r}\dot{\theta}}{r} - \dfrac{\mu_1 g\cos{\theta}}{M_1 r} + \dfrac{Q_{\theta}}{M_1r^2}.
\end{aligned}

Where
            
* $r$ is the length of the pendulum rod (in metres).
* $\theta$ is the angle of the pendulum relative to the positive $x$-axis.
* $g$ is the acceleration due to gravity (in metres per second squared).
* $l_0$ is the rest length of the pendulum.
* $Q_{r} = -(b_b+c_b v_b)\dot{r} - (b_r + c_r v_r)\dfrac{\dot{r}}{4}$.
* $v_b = \sqrt{\dot{r}^2+r^2\dot{\theta}^2}$ and $v_r = \dfrac{1}{2}v_b$. 
* $Q_{\theta} = -(b_b+c_b v_b)r^2\dot{\theta} - (b_r + c_r v_r)\dfrac{r^2\dot{\theta}}{4}$.

In this calculation, we assumed that the dissipative forces on the rod could be calculated using a centre-of-mass approach.

~~~
{{ insert template.html }}
~~~