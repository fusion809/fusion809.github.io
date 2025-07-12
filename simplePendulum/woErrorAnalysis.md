@def title="Simple pendulum Runge-Kutta-Fehlberg solver without error analysis"
@def SP=true;
@def SP.wo=true;
@def hassim=true;
@def params = (g=(val=9.81, desc="Acceleration due to gravity."), l=(val=1, desc="Length of pendulum in metres."), tf=(val=9.4713677903049498, desc="End time for the simulation in seconds. Default is four times the period of the problem."), theta0=(val=0, desc="Initial angle relative to the positive \\(x\\)-axis."), dtheta0=(val=0, desc="Initial rate of change of the aforementioned angle."), N=(val=1e6, desc="Number of quadrature nodes used to approximate period"))
@def ids=["tableOutputs", "phasePlot", "timePlot", "phasePlot", "animation", "animationPhase", "animation"]
@def funcs=["fillTableSP(readInputs())","removeTable()","generatePhasePlot(solveProblemSP(readInputs()))","removePhasePlot()","generateTimePlot(solveProblemSP(readInputs()))","removeTimePlot()","generatePlots(readInputs())","removePlots()","generateAnimation()","removeAnimation()","generatePhaseAnimation()","removePhaseAnimation()","generateAnimations()","removeAnimations()"]
@def labels=["Tabulate the solution","Remove the solution table","Generate a phase plot","Remove phase plot","Generate a time plot","Remove time plot","Generate all solution plots","Remove all solution plots","Generate an animation","Remove animation","Generate a phase plot animation","Remove phase plot animation","Generate all animations","Remove all animations"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of the simple pendulum:

\begin{align}
\dfrac{d^2 \theta}{dt^2} = -\dfrac{g}{l} \cos{\theta}
\end{align}

where $g$ is the acceleration due to gravity in metres per second squared, $l$ is the length of the pendulum in metres, $\theta$ is the angle from the positive $x$ axis (in radians) and $t$ is the time in seconds. Below you can specify the various parameters for the problem we will solve. 

Our $\theta$ approximations are substituted in, and our $\dot{\theta}$ RKF45 approximation is subtracted from this value. From this equation, the period $T$ of the problem is approximated when the conditions for periodicity are satisfied, namely using the equation

\begin{align*}
T &= 2 \left|\int_{\theta_\mathrm{min}}^{\theta_\mathrm{max}} \dfrac{d\theta}{\sqrt{\dot{\theta}_0^2 + \dfrac{2g}{l} (\sin{\theta_0} - \sin{\theta})}}\right|
\end{align*}

where $\theta_\mathrm{min}$ and $\theta_\mathrm{max}$ are the two closest values for which $\dot{\theta} = 0$. $T$ is calculated using [Chebyshev-Gauss quadrature](https://en.wikipedia.org/wiki/Chebyshev-Gauss_quadrature) (Simpson's rule could not be used as there are unremovable singularities at the endpoints which makes Simpson's rule markedly less accurate).

~~~
{{ insert template_SP.html}}
~~~