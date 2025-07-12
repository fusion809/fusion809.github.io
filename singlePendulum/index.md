@def title="Single pendulum"
@def hassim=true
@def funcs = ["fillTable(readInputs(), ['θ', 'dθ/dt'])","removeTable()", "generatePhasePlot(solveProblem(RKF45, readInputs()))","removePhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))","removeTimePlot()", "generatePlots(readInputs())","removePlots()", "generateAnimation()","removeAnimation()", "generatePhaseAnimation()","removePhaseAnimation()", "generateAnimations()","removeAnimations()"]
@def params = (g=(val=9.81, desc="Acceleration due to gravity in \\(\\mathrm{m} \\cdot \\mathrm{s}^{-2}\\)."), l=(val=1, desc="Length of the pendulum rod in metres (m)."),mr=(val=1, desc="Mass of rod."), mb=(val=1, desc="Mass of bob."),br=(val=0.05, desc="Linear dissipation coefficient of rod."),bb=(val=0.05, desc="Linear dissipation coefficient of bob."),cr=(val=0.02, desc="Quadratic dissipation coefficient of rod."),cb=(val=0.02, desc="Quadratic dissipation coefficient of bob."),t0=(val=0, desc="Starting time for the simulation in seconds (s)."),tf=(val=100, desc="End time for the simulation in seconds. The default \\(t_f - t_0\\) value is equal to four times the period \\(T\\) of the problem."),theta0=(val=0, desc="Angle from the positive \\(x\\)-axis in radians at time \\(t_0\\)."),dtheta0=(val=0, desc="First derivative of \\(\\theta\\) with respect to \\(t\\) at \\(t_0\\) (in \\(\\mathrm{radians}\\cdot \\mathrm{s}^{-1}\\))."),epsilon=(val=1e-7,))
@def labels = ["Tabulate the solution","Remove the solution table","Generate a phase plot","Remove plot","Generate a time plot","Remove plot","Generate all solution plots","Remove all solution plots","Generate an animation of the system","Remove animation","Generate a phase plot animation","Remove phase plot animation","Generate all animations","Remove all animations"]
@def ids = ["tableOutputs", "phasePlot", "timePlot", "phasePlot", "animation", "animationPhase", "animation"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of a single pendulum. To calculate the equations of motion, we must first determine the Lagrangian which is defined as the kinetic energy minus the potential energy. The kinetic energy of the bob is simply its mass times its velocity squared over two. The kinetic energy of the rod has a term like this, too, (although the velocity we will focus on will be the velocity of its centre of mass) but we must also consider its rotational velocity. As for the potential energy of the bob and rod, it will simply be the mass of them times by the $y$-coordinate of its centre of mass.  

\begin{aligned}
\mathcal{L} &= \dfrac{m_b}{2} |\vec{v}_b|^2 + \dfrac{1}{2} I_{\mathrm{cm}}\omega^2 + \dfrac{m_r}{2} |\vec{v}_{r, \mathrm{cm}}|^2 - m_bg l\sin{\theta} -\dfrac{m_r gl\sin{\theta}}{2}.
\end{aligned}

Where $g$ is the acceleration due to gravity in metres per second squared, $l$ is the length of the pendulum in metres, $\theta$ is the angle from the positive $x$ axis (in radians). $m_b$ is the mass of the bob in kilograms. $\vec{v}_b$ is the velocity vector for the bob in metres per second. $\vec{v}_{r,\mathrm{cm}}$ is the velocity of the rod's centre of mass in metres per second. $m_r$ is the mass of the rod in kilograms. $t$ is the time in seconds. $I_{\mathrm{cm}}$ is the moment of inertia around the centre of mass of the rod. It would be $\dfrac{m_r l^2}{12}$.[^1] Velocity squared, for the bob, in polar coordinates is given by $\dot{l}^2 + l^2\dot{\theta}^2$ but here $l$ is a constant and hence its time derivative is zero. As for the velocity of the rod's centre of mass, we will assume the rod is uniform in mass and hence its centre of mass is in its centre, which we will assume is midway between the bob and the origin. Hence its velocity will be half that of the bob. 

\begin{aligned}
\mathcal{L} &= \dfrac{m_b}{2} l^2\dot{\theta}^2 + \dfrac{m_r}{24}l^2\dot{\theta}^2 + \dfrac{m_r}{8}l^2\dot{\theta}^2 - \left(m_b+\dfrac{m_r}{2}\right)gl\sin{\theta} \\
&= \left[\dfrac{m_b}{2} + \dfrac{m_r}{6}\right]l^2\dot{\theta}^2 - \left(m_b+\dfrac{m_r}{2}\right)gl\sin{\theta}.
\end{aligned}

Hence its generalized momentum, momentum's time derivative and force are:

\begin{aligned}
\dfrac{\partial \mathcal{L}}{\partial \dot{\theta}} &= \left[m_b + \dfrac{m_r}{3}\right]l^2\dot{\theta}\\
\dfrac{d}{dt}\dfrac{\partial \mathcal{L}}{\partial \dot{\theta}} &= \left[m_b + \dfrac{m_r}{3}\right]l^2\ddot{\theta} \\
\dfrac{\partial \mathcal{L}}{\partial \theta} &= - \left(m_b+\dfrac{m_r}{2}\right)gl\cos{\theta}.
\end{aligned}

Next we will define $\dfrac{\delta' f}{\delta' q_i}$ as the negative functional derivative of $f(q_i, \dot{q}_i, t)$ with respect to $q_i$.

\begin{aligned}
\dfrac{\delta' \mathcal{L}}{\delta' \theta} &= \dfrac{d}{dt}\dfrac{\partial \mathcal{L}}{\partial \dot{\theta}} - \dfrac{\partial \mathcal{L}}{\partial \theta} \\
&= \left[m_b + \dfrac{m_r}{3}\right]l^2\ddot{\theta} + \left(m_b+\dfrac{m_r}{2}\right)gl\cos{\theta}.
\end{aligned}

Next we will calculate the generalized dissipation force $Q_i$ which is defined as being the dot product of the dissipation force vector with generalized basis vector for each component of a physical system. Here there are two components of the system: the rod and the bob. The dissipation force vector will be assumed to be $\vec{F}_{D,j} = -(b_j + c_j|\vec{v}_j|)\vec{v}_j$.

 
\begin{aligned}
Q_{\theta} &= -(b_r+c_r|\vec{v}_r|)\vec{v}_r \cdot \dfrac{\partial \vec{r}_r}{\partial \theta} -(b_b+c_b|\vec{v}_b|)\vec{v}_b \cdot \dfrac{\partial \vec{r}_b}{\partial \theta}\\
&= -\left(b_r+c_r\dfrac{l|\dot{\theta}|}{2}\right)\dfrac{l\dot{\theta}}{2}\begin{bmatrix}
-\sin{\theta} \\
\cos{\theta}
\end{bmatrix} \cdot \dfrac{l}{2} \begin{bmatrix}
-\sin{\theta} \\
\cos{\theta}
\end{bmatrix} -(b_{b} +c_{b} l|\dot{\theta}|)l\dot{\theta}\begin{bmatrix}
-\sin{\theta} \\
\cos{\theta}
\end{bmatrix} \cdot l \begin{bmatrix}
-\sin{\theta} \\
\cos{\theta}
\end{bmatrix} \\
&= -\left(b_r+c_r \dfrac{l|\dot{\theta}|}{2}\right)\dfrac{l^2\dot{\theta}}{4} -(b_b+c_b l|\dot{\theta}|)l^2\dot{\theta}.
\end{aligned}


Hence the Euler-Lagrange equation with dissipative forces is

\begin{aligned}
\dfrac{\delta' \mathcal{L}}{\delta' \theta} &= Q_{\theta} \\
\left[m_b + \dfrac{m_r}{3}\right]l^2\ddot{\theta} + \left(m_b+\dfrac{m_r}{2}\right)gl\cos{\theta} &= -\left(b_r+c_r \dfrac{l|\dot{\theta}|}{2}\right)\dfrac{l^2\dot{\theta}}{4} -(b_b+c_b l|\dot{\theta}|)l^2\dot{\theta}.
\end{aligned}


Dividing by the coefficient of $\ddot{\theta}$, namely $\left[m_b + \dfrac{m_r}{3}\right]l^2$ yields

\begin{aligned}
\ddot{\theta} + \dfrac{m_b+\dfrac{m_r}{2}}{\left[m_b + \dfrac{m_r}{3}\right]l}g\cos{\theta} &= -\dfrac{\left(b_r+c_r \dfrac{l|\dot{\theta}|}{2}\right)\dfrac{\dot{\theta}}{4} +(b_b+c_b l|\dot{\theta}|)\dot{\theta}}{m_b + \dfrac{m_r}{3}}.
\end{aligned}

Isolating $\ddot{\theta}$ yields

\begin{aligned}
\ddot{\theta} &= -\dfrac{m_b+\dfrac{m_r}{2}}{\left[m_b + \dfrac{m_r}{3}\right]l}g\cos{\theta} -\dfrac{\left(b_r+c_r \dfrac{l|\dot{\theta}|}{2}\right)\dfrac{\dot{\theta}}{4} +(b_b+c_b l|\dot{\theta}|)\dot{\theta}}{m_b + \dfrac{m_r}{3}} \\
&= -\dfrac{6m_b+3m_r}{\left[6m_b + 2m_r\right]l}g\cos{\theta} -\dfrac{(6b_r+24b_b)\dot{\theta}+(3c_r+24c_b)l|\dot{\theta}|\dot{\theta}}{24m_b + 8m_r}.
\end{aligned}

Below you can specify the various parameters for the problem we will solve. 

# References
[^1] Dourmashkin, P (2022). [16.3 Rotational Kinetic Energy and Moment of Inertia](https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Classical_Mechanics_(Dourmashkin)/16%3A_Two_Dimensional_Rotational_Kinematics/16.03%3A_Rotational_Kinetic_Energy_and_Moment_of_Inertia) in [*Classical Mechanics (Dourmashkin)*](https://commons.libretexts.org/book/phys-24413). Massachusetts Institute of Technology.

~~~
{{ insert parameter_form.html }}
{{ insert button_table.html }}
{{ insert simple_pendulum_output.html }}
~~~