### Without air resistance
To apply Equations \eqref{ELE}, we will treat it as a simple pendulum with no friction, a massless rod and no air resistance. The pendulum will consist of a rod with a bob at one end and freely rotating around a central point. We will define this central point as being the origin in Cartesian coordinates. Say the rod is of fixed length $l$ metres. The pendulum bob will be of mass $m$ kilograms. As you can see in Figure \figref{fig1}, we have defined $\theta$ as the angle between the rod and the positive $x$-axis going counter-clockwise. Hence the $x$ and $y$ coordinates of the pendulum bobs are:
\begin{align*}
    x &= l\cos{\theta} \\
    y &= l\sin{\theta}.
\end{align*}
Hence: 
\begin{align*}
    \dot{x} &= -l\dot{\theta}\sin{\theta} \\
    \dot{y} &= l\dot{\theta}\cos{\theta}.
\end{align*}
As for velocity squared, it is defined as:
\begin{align*}
    v^2 &= \dot{x}^2+\dot{y}^2 \\
    &= l^2\dot{\theta}^2\sin^2{\theta} + l^2\dot{\theta}^2\cos^2{\theta} \\
    &= l^2\dot{\theta}^2.
\end{align*}
Hence $v=l\dot{\theta}$. The kinetic energy is:
\begin{align*}
    T &= \dfrac{mv^2}{2} \\
    &= \dfrac{ml^2\dot{\theta}^2}{2}.
\end{align*}
As for potential energy, given this bob is moving solely under the force of gravity:
\begin{align*}
    V &= mgy \\
    &= mgl\sin{\theta}.
\end{align*}
Hence:
\begin{align*}
    L &= T - V \\
    &= \dfrac{ml^2\dot{\theta}^2}{2} - mgl\sin{\theta}.
\end{align*}
$\dfrac{\partial L}{\partial \dot{\theta}}$ is hence:
\begin{align*}
    \dfrac{\partial L}{\partial \dot{\theta}} &= ml^2 \dot{\theta}.
\end{align*}
And its time derivative is:
\begin{align*}
    \dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{\theta}} = ml^2 \ddot{\theta}.
\end{align*}
As for the $\dfrac{\partial L}{\partial \theta}$ term it is:
\begin{align*}
    \dfrac{\partial L}{\partial \theta} &= -mgl\cos{\theta}.
\end{align*}
So Equations \eqref{ELE} becomes:
\begin{align*}
    ml^2 \ddot{\theta} + mgl\cos{\theta} &= 0.
\end{align*}
Dividing by $ml^2$ and taking the second term to the right-hand side yields:
\begin{align}
\ddot{\theta} &= -\dfrac{g}{l}\cos{\theta}.\label{simppen}
\end{align}
Equation \eqref{simppen} has, as its highest order derivative, a second-order derivative. It also only involves ordinary derivatives (not partial derivatives). Hence it is a second-order ordinary differential equation (ODE). Second-order ODEs are typically obtained in both Newtonian and Lagrangian mechanics. The number of second-order ODEs obtained for a given system depends on the number of degrees of freedom of the system. 

Technically speaking it is possible to analytically solve Equation \eqref{simppen}, but it requires the use of non-standard functions called elliptic integrals. If you would like to check out numerical solutions to this equation, see [this page](/simplePendulum/woErrorAnalysis).