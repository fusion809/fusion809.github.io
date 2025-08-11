@def title="Coupled duffing oscillator solver"
@def hassim=true
@def vars = ["x", "dx", "y", "dy"]
@def params = (alpha=(val=1, desc="Linear stiffness parameter."),beta=(val=5, desc="Nonlinearity in restoring force parameter."),gamma=(val=8, desc="Amplitude of periodic driving force."),delta=(val=0.02, desc="Damping parameter."),omega=(val=0.5, desc="Angular frequency of periodic driving force."),kappa=(val=0.1, desc="Coupling term coefficient."), tf=(val=20, desc="End time (in seconds or s) for the simulation."),x0=(val=1, desc="Initial \\(x\\) coordinate in metres (m) ."),dx0=(val=0, desc="Initial value of \\(\\dot{x}\\) (\\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\))."),y0=(val=1, desc="Initial \\(y\\) coordinate in metres (m) ."),dy0=(val=0, desc="Initial value of \\(\\dot{y}\\) (\\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\))."))
@def ids = ["tableOutputs", "xPhasePlot", "yPhasePlot", "XYPlot", "timePlot", "XYPlot", "animationX", "animationY", "animationXY", "animationXY", "animationXY"];
@def funcs = ["generateTable()", "removeTable()", "generateXPhasePlot(solveProblem(RKF45, readInputs()))", "removeXPhasePlot()", "generateYPhasePlot(solveProblem(RKF45, readInputs()))", "removeYPhasePlot()", "generateXYPlot(solveProblem(RKF45, readInputs()))", "removeXYPlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))", "removeTimePlot()", "generatePlots(readInputs())", "removePlots()", "generateXAnimation()", "removeXAnimation()","generateYAnimation()", "removeYAnimation()","generateAnimation()", "removeAnimation()","generateAnimations()", "removeAnimations()", "generateAllOutputs()", "removeAllOutputs()"];
@def labels = ["Tabulate the solution", "Remove the solution table", "Generate an X phase plot", "Remove X phase plot", "Generate a Y phase plot", "Remove Y phase plot", "Generate an XY plot", "Remove XY plot", "Generate a time plot", "Remove time plot", "Generate all solution plots", "Remove all solution plots", "Generate an X phase animation", "Remove X phase animation", "Generate a Y phase animation", "Remove Y phase animation", "Generate an XY animation", "Remove XY animation", "Generate all animations", "Remove all animations", "Generate all outputs", "Remove all outputs"];

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of the [coupled Duffing oscillator](https://arxiv.org/pdf/nlin/0304022):

\begin{align*}
\ddot{x} + \delta \dot{x} + \alpha x + \beta x^3 + \kappa (x-y) &= \gamma \cos{(\omega t)} \\
\ddot{y} + \delta \dot{y} + \alpha y + \beta y^3 + \kappa (y-x) &= 0.
\end{align*}

Below you can specify these various parameters, as well as the initial conditions and starting and end times. The default values give chaotic behaviour.

~~~
{{ insert template.html}}
~~~