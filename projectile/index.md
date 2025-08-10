@def title="Projectile motion solver"
@def hassim=true
@def vars = ["x", "dx", "y", "dy"]
@def params = (g=(val=9.81, desc="Gravitational acceleration in metres per second squared."),b=(val=0.1, desc="Linear dissipation coefficient divided by projectile mass"), c=(val=0.0001, desc="Quadratic dissipation coefficient divided by projectile mass"), tf=(val=16, desc="End time for the simulation in seconds."),x0=(val=0, desc="Initial horizontal position of projectile."),dx0=(val=100, desc="Initial horizontal velocity of projectile"), y0=(val=0, desc="Initial vertical position of projectile."),dy0=(val=100, desc="Initial vertical velocity of projectile."), N=(val=100001, desc="Number of steps used for exact solution."), epsilon=(val=1e-9,), tolType=(val=1,), t1=(val=14,))
@def ids = ["tableOutputs", "xPhasePlot", "yPhasePlot", "xYPlot", "timePlot", "xPhasePlot", "animationX", "animationY", "animationXY", "animationXY", "animationXY"];
@def funcs = ["generateTable()", "removeTable()", "generateXPhasePlot(solve(readInputs()))", "removeXPhasePlot()", "generateYPhasePlot(solve(readInputs()))", "removeYPhasePlot()", "generateXYPlot(solve(readInputs()))", "removeXYPlot()", "generateTimePlot(solve(readInputs()))", "removeTimePlot()", "generatePlots(readInputs())", "removePlots()", "generateXAnimation()", "removeXAnimation()","generateYAnimation()", "removeYAnimation()", "generateXYAnimation()", "removeXYAnimation()", "generateAnimations()", "removeAnimations()", "generateAllOutputs()", "removeAllOutputs()"];
@def labels = ["Tabulate the solution", "Remove the solution table", "Generate an X phase plot", "Remove X phase plot", "Generate a Y phase plot", "Remove Y phase plot", "Generate an XY plot", "Remove XY plot", "Generate a time plot", "Remove time plot", "Generate all solution plots", "Remove all solution plots", "Generate an X phase animation", "Remove X phase animation", "Generate a Y phase animation", "Remove Y phase animation", "Generate XY animation", "Remove XY animation", "Generate all animations", "Remove all animations", "Generate all outputs", "Remove all outputs"];

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the equations of motion for a projectile subjected to constant gravity and linear and quadratic drag

\begin{aligned}
\dfrac{d^2 x}{dt^2} &= -\dfrac{dx}{dt}\left(b+c\sqrt{\left(\dfrac{dx}{dt}\right)^2+\left(\dfrac{dy}{dt}\right)^2}\right) \\
\dfrac{d^2 y}{dt^2} &= -g-\dfrac{dy}{dt}\left(b+c\sqrt{\left(\dfrac{dx}{dt}\right)^2+\left(\dfrac{dy}{dt}\right)^2}\right).
\end{aligned}

Where $g$ is the gravitational acceleration constant, $b$ is the linear drag coefficient for the projectile divided by its mass and $c$ is the quadratic drag coefficient for the projectile divided by its mass. $y < 0$ is not allowed in the solver, as it is assumed that $y=0$ is the ground. If $c=0$, exact solutions are used. Namely, if $b\neq0$

\begin{aligned}
x(t) &= x_0 + \dfrac{\dot{x}_0}{b}\left(1-e^{-bt}\right) \\
y(t) &= y_0 - \dfrac{gt^2}{2} + \dfrac{\dot{y}_0}{b}\left(1-e^{-bt}\right),
\end{aligned}

whereas if $b=0$

\begin{aligned}
x(t) &= x_0 + \dot{x}_0 t \\
y(t) &= y_0 + \dot{y}_0 t - \dfrac{gt^2}{2}.
\end{aligned}

~~~
{{ insert template.html }}
~~~