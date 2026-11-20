### Projectile motion
Let us analyse what happens when a projectile is fired into the air. Let $x$ represent the horizontal direction and $y$ represent the vertical direction. Our potential energy here would be $mgy$, hence it is independent of the generalized velocities and we can use the $H=T+V$ definition of the Hamiltonian.
\begin{align*}
H &= \dfrac{m}{2}(\dot{x}^2+\dot{y}^2) + mgy.
\end{align*}
As for our canonical momentum:
\begin{align*}
p_x &= m\dot{x} \hspace{0.2cm}\implies \dot{x} &&= \dfrac{p_x}{m}\\
p_y &= m\dot{y} \hspace{0.2cm}\implies \dot{y} &&= \dfrac{p_y}{m}.
\end{align*}
Hence our Hamiltonian can be rewritten as:
\begin{align*}
H &= \dfrac{p_x^2+p_y^2}{2m} + mgy.
\end{align*}
As for air resistance it could be approximated by:
\begin{align*}
\mathbf{F}_D &= -\dfrac{\rho C_D A|v|\textbf{v}}{2} \\
&= -\dfrac{\rho C_D A\sqrt{\dot{x}^2+\dot{y}^2}}{2} (\dot{x}, \dot{y}).
\end{align*}
As for our covariant basis vectors they are, since we are using Cartesian coordinates:
\begin{align*}
\dfrac{\partial \mathbf{r}}{\partial x} &= (1, 0) \\
\dfrac{\partial \mathbf{r}}{\partial y} &= (0, 1).
\end{align*}
Hence our generalized non-conservative forces are:
\begin{align*}
Q_x &= -\dfrac{\rho C_D A\sqrt{\dot{x}^2+\dot{y}^2}}{2} (\dot{x}, \dot{y}) \cdot (1,0) \\
&= -\dfrac{\rho C_D A\dot{x}\sqrt{\dot{x}^2+\dot{y}^2}}{2}\\
&= -\dfrac{\rho C_D A\dfrac{p_x}{m}\sqrt{\dfrac{p_x^2}{m^2}+\dfrac{p_y^2}{m^2}}}{2}\\
&= -\dfrac{\rho C_D A p_x\sqrt{p_x^2+p_y^2}}{2m^2}\\
Q_y &= -\dfrac{\rho C_D A\sqrt{\dot{x}^2+\dot{y}^2}}{2} (\dot{x}, \dot{y}) \cdot (0, 1)\\
&= -\dfrac{\rho C_D A\dot{y}\sqrt{\dot{x}^2+\dot{y}^2}}{2} \\
&= -\dfrac{\rho C_D A\dfrac{p_y}{m}\sqrt{\dfrac{p_x^2}{m^2}+\dfrac{p_y^2}{m^2}}}{2}\\
&= -\dfrac{\rho C_D A p_y\sqrt{p_x^2+p_y^2}}{2m^2}.
\end{align*}
Our equations of motion are:
\begin{align*}
\dot{x} &= \dfrac{\partial H}{\partial p_x} \\
\dot{y} &= \dfrac{\partial H}{\partial p_y} \\
\dot{p}_x &= -\dfrac{\partial H}{\partial x} + Q_x \\
\dot{p}_y &= -\dfrac{\partial H}{\partial y} + Q_y.
\end{align*}
Hence for projectile motion they become:
\begin{align}
\dot{x} &= \dfrac{p_x}{m} \label{dxpm}\\
\dot{y} &= \dfrac{p_y}{m} \label{dypm}\\
\dot{p}_x &= -\dfrac{\rho C_D A p_x\sqrt{p_x^2+p_y^2}}{2m^2} \label{dpxpm}\\
\dot{p}_y &= -mg -\dfrac{\rho C_D A p_y\sqrt{p_x^2+p_y^2}}{2m^2}.\label{dpypm}
\end{align}

Code Listing \coderef{projmotion} is Julia code for implementing an animation based on Equations \eqref{dxpm}, \eqref{dypm}, \eqref{dpxpm} and \eqref{dpypm}.

\codelisting{projmotion}{Projectile motion simulation code.}{./code/projmotion}{fidello/code/projmotion.jl}
Below is an animation of projectile motion created using the above Julia code.

\figure{projectile_motion}{right}{500px}{/projectile_motion.gif}{Projectile motion animation.}