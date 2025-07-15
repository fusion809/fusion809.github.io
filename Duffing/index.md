@def title="Duffing JavaScript integrator"
@def hassim=true
@def params = (alpha=(val=1, desc="Linear stiffness parameter."),beta=(val=5, desc="Nonlinearity in restoring force parameter."),gamma=(val=8, desc="Amplitude of periodic driving force."),delta=(val=0.02, desc="Damping parameter."),omega=(val=0.5, desc="Angular frequency of periodic driving force."),tf=(val=20, desc="End time for the simulation in seconds.."),x0=(val=1, desc="Initial x coordinate in metres (m) ."),dx0=(val=0, desc="Initial first derivative of \\(x\\) with respect to \\(t\\) (in \\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\))."))
@def type="2D"

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of the [Duffing oscillator](https://en.wikipedia.org/wiki/Duffing_oscillator):

\begin{align*}
\ddot{x} + \delta \dot{x} + \alpha x + \beta x^3 = \gamma \cos{(\omega t)}
\end{align*}

Below you can specify these various parameters, as well as the initial conditions and starting and end times. The default values give chaotic behaviour.

~~~
{{ insert template.html}}
~~~