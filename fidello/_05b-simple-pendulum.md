## Example: simple pendulum
The Hamiltonian for the simple pendulum is:
\begin{align}
H &= T + V \\
&= \dfrac{ml^2 \dot{\theta}^2}{2} + mgl\sin{\theta}. \label{Hsimp}
\end{align}
Problem is, however, that this is written in terms of $q_i$ and $\dot{q}_i$ not $q_i$ and $p_i$. First we need to determine what our canonical momentum is. Equations \eqref{canmom} yields:
\begin{align}
p_{\theta} &= \dfrac{\partial L}{\partial \dot{\theta}} \nonumber \\
p_{\theta} &= ml^2\dot{\theta}. \label{pthetasimp}
\end{align}
Rewriting Equation \eqref{pthetasimp} with $\dot{\theta}$ as its subject:
\begin{align}
\dot{\theta} &= \dfrac{p_{\theta}}{ml^2}. \label{thetadotsimp}
\end{align}
Substituting Equation \eqref{thetadotsimp} into Equation \eqref{Hsimp} yields:
\begin{align*}
H &= \dfrac{ml^2}{2} \left( \dfrac{p_{\theta}}{ml^2} \right)^2 + mgl\sin{\theta} \\
&= \dfrac{ml^2}{2} \dfrac{p_{\theta}^2}{m^2l^4} + mgl\sin{\theta}.
\end{align*}
Cancelling out common factors in the numerator and denominator of the kinetic energy term yields:
\begin{align}
H &= \dfrac{p_{\theta}^2}{2ml^2} + mgl\sin{\theta}. \label{Hsimpf}
\end{align}
Hence the equations of motion are:
\begin{align}
\dot{p}_{\theta} &= -\dfrac{\partial H}{\partial \theta} \nonumber \\
 &= - mgl\cos{\theta}\label{pdotthetasimp}\\
\dot{\theta} &= \dfrac{\partial H}{\partial p_{\theta}} \nonumber\\
&= \dfrac{p_{\theta}}{ml^2}.\label{thetadotsimp}
\end{align}
As we can see, Hamiltonian mechanics yields a series of first-order ordinary differential equations (ODEs) instead of second-order ODEs. Which is important because Lagrangian mechanics and Newtonian mechanics both yield second-order ODE systems. That being said, if you rewrite Equation \eqref{thetadotsimp} to make $p_{\theta}$ its subject and substitute it into Equation \eqref{pdotthetasimp} and rearrange it to make $\ddot{\theta}$ the subject, you will obtain Equation \eqref{simppen}.