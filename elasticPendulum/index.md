@def title="Elastic pendulum solver"
@def hassim=true
@def ids = ["phasePlotZTheta","phasePlotZZDot","phasePlotThetaThetaDot","timePlot","animation","animationThetaPhase","animationZPhase"]
@def params = (g=(val=9.81, desc="Problem parameter."),l0=(val=1, desc="Problem parameter."),k=(val=1, desc="Problem parameter."),m=(val=1, desc="Problem parameter."),tf=(val=60, desc="End time for the simulation in seconds."),z0=(val=1, desc="Initial value of \\(z\\)."),dz0=(val=1, desc="Initial value of \\(\\dot{z}\\)."),theta0=(val=0, desc="Initial value of \\(\\theta\\)."),dtheta0=(val=0, desc="Initial value of \\(\\dot{\\theta}\\)."),epsilon=(val=1e-9,))
@def funcs = ["fillTable(readInputs(), ['z', 'dz/dt', '&theta;', 'd&theta;/dt'])","removeTable()", "generateZThetaPhasePlot(solveProblem(RKF45, readInputs()))","removeZThetaPhasePlot()", "generateZZDotPhasePlot(solveProblem(RKF45, readInputs()))","removeZZDotPhasePlot()", "generateThetaThetaDotPhasePlot(solveProblem(RKF45, readInputs()))","removeThetaThetaDotPhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))","removeTimePlot()", "generatePlots(readInputs())","removePlots()", "generateAnimation()","removeAnimation()", "generateThetaPhaseAnimation()","removeThetaPhaseAnimation()", "generateZPhaseAnimation()","removeZPhaseAnimation()", "generateAnimations()","removeAnimations()"]
@def labels = ["Tabulate the solution","Remove the solution table","Generate a \\(\\theta\\) against \\(z\\) phase plot","Remove \\(\\theta\\) against \\(z\\) phase plot","Generate a \\(\\dot{z}\\) against \\(z\\) phase plot","Remove \\(\\dot{z}\\) against \\(z\\) phase plot","Generate a \\(\\dot{\\theta}\\) against \\(\\theta\\) phase plot","Remove \\(\\dot{\\theta}\\) against \\(\\theta\\) phase plot","Generate a \\(z\\), \\(\\dot{z}\\), \\(\\theta\\) and \\(\\dot{\\theta}\\) against time plot","Remove \\(z\\), \\(\\dot{z}\\), \\(\\theta\\) and \\(\\dot{\\theta}\\) against time plot","Generate all solution plots","Remove all plots","Generate an animation","Remove animation","Generate a \\(\\theta\\) phase plot animation","Remove \\(\\theta\\) phase plot animation","Generate a \\(z\\) phase plot animation","Remove \\(z\\) phase plot animation","Generate all animations","Remove all animations"]

This webpage uses the <a href='/RKF45/' link='_blank'>Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)</a> to approximate the solution to the problem of the <a href='https://en.wikipedia.org/wiki/Elastic_pendulum' link='_blank'>elastic pendulum</a>

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