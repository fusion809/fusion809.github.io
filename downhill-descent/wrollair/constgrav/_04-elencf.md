### Euler-Lagrange equation with non-conservative forces
The Euler-Lagrange equation with dissipative forces in the $s$ coordinate is:
\begin{align}
	\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{s}} - \dfrac{\partial L}{\partial s} &= Q_s. \label{ELENCF}
\end{align}
The left-hand side of Equation \eqref{ELENCF} can also be denoted as $\dfrac{\delta S}{\delta s}$, the functional derivative of $S$ with respect to $s$. $S$ is the action and is defined by:
\begin{align*}
S=\int_{t_0}^{t_f} L(s, \dot{s}, t) dt
\end{align*}
where $t_0$ is the starting time and $t_f$ is the finishing time. First, we shall calculate $\dfrac{\partial L}{\partial \dot{s}}$, the canonical momentum:
\begin{align*}
	\dfrac{\partial L}{\partial \dot{s}} &= m\dot{s}.
\end{align*}
Hence:
\begin{align}
	\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{s}} &= m\ddot{s}. \label{ma}
\end{align}
As for the second term in Equation \eqref{ELENCF}, it is:
\begin{align*}
	\dfrac{\partial L}{\partial s} &= mg\sin{\theta(s)}.
\end{align*}
So the functional derivative is therefore:
\begin{align*}
\dfrac{\delta S}{\delta s} &= m\ddot{s} - mg\sin{\theta(s)}.
\end{align*}
Hence Equation \eqref{ELENCF} becomes:
\begin{align*}
	m\ddot{s} - mg\sin{\theta(s)} &= -mgC_{rr} \cos{\theta(s)} - \dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
Rearranging and factorizing yields:
\begin{align*}
	m\ddot{s} &= mg(\sin{\theta(s)}+C_{rr}\cos{\theta(s)})-\dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
Dividing both sides by $m$ yields:
\begin{align}
	\ddot{s} &= g(\sin{\theta(s)}+C_{rr}\cos{\theta(s)})-\dfrac{\rho C_D A\dot{s}^2}{2m}\label{eq1}
\end{align}
assuming constant gravitational acceleration.
