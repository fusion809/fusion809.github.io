### Computing Equation \eqref{ELENCF}
The conservative force term of the Euler-Lagrange equation is therefore:
\begin{align}
	\dfrac{\partial L}{\partial s} &= \dfrac{GMm\sin{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2}. \label{Fnc}
\end{align}
Hence the functional derivative of $S$ with respect to $s$ is therefore:
\begin{align*}
	\dfrac{\delta S}{\delta s} &= m\ddot{s} - \dfrac{GMm\sin{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2}.
\end{align*}
The time derivative of generalized momentum is the same as above. As for $Q_s$ it is:
\begin{align}
	Q_s &= -\dfrac{GMmC_{rr}\cos{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} - \dfrac{\rho C_D A \dot{s}^2}{2}. \label{Qs}
\end{align}
Substituting Equations \eqref{Fnc} and \eqref{Qs} into Equation \eqref{ELENCF} yields:
\begin{align*}
	m\ddot{s} - \dfrac{GMm\sin{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} &= -\dfrac{GMmC_{rr}\cos{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} - \dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
Rearranging and factorizing:
\begin{align*}
	m\ddot{s} &= \dfrac{GMm}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right] - \dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
Dividing by $m$:
\begin{align}
	\ddot{s} &= \dfrac{GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right] - \dfrac{\rho C_D A\dot{s}^2}{2m}.\label{eq2}
\end{align}