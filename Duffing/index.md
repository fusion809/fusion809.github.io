@def title="Duffing JavaScript integrator"
@def hassim=true
@def params = (alpha=(val=1, desc="Parameter."),beta=(val=5, desc="Parameter."),gamma=(val=8, desc="Parameter."),delta=(val=0.02, desc="Parameter."),omega=(val=0.5, desc="Parameter."),t0=(val=0, desc="Starting time for the simulation in seconds (s)."),tf=(val=20, desc="End time for the simulation in seconds.."),x0=(val=1, desc="x coordinate in metres (m) at time \\(t_0\\)."),dx0=(val=0, desc="First derivative of \\(x\\) with respect to \\(t\\) at \\(t_0\\) (in \\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\))."))

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of the Duffing oscillator:

\begin{align*}
\ddot{x} + \delta \dot{x} + \alpha x + \beta x^3 = \gamma \cos{(\omega t)}
\end{align*}

Below you can specify these various parameters, as well as the initial conditions and starting and end times. The default values give chaotic behaviour.

~~~
{{ insert parameter_form.html }}
<!--Buttons-->
{{ insert 2D_buttons.html }}

<!--Where the table and plot goes-->
{{ insert 2D_output.html }}
~~~