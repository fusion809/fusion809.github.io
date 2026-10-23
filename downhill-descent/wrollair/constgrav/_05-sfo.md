### Simplification to first order
\begin{align*}
	\ddot{s} &= \dfrac{d\dot{s}}{ds} \dfrac{ds}{dt} \\
	&= \dfrac{ds}{dt} \dfrac{d\dot{s}}{ds} \\
	&= \dot{s} \dfrac{d\dot{s}}{ds}.
\end{align*}
Hence Equation \eqref{eq1} becomes, after dividing by $m$:
\begin{align*}
	\dot{s} \dfrac{d\dot{s}}{ds} &= g(\sin{\theta(s)} - C_{rr} \cos{\theta(s)}) - \dfrac{\rho C_D A\dot{s}^2}{2m}.
\end{align*}
Let $u=\dot{s}^2$, then:
\begin{align*}
	\dfrac{du}{ds} &= \dfrac{d(\dot{s}^2)}{ds}
\end{align*}
Hence by the chain rule:
\begin{align*}
	\dfrac{du}{ds} &= 2\dot{s}\dfrac{d\dot{s}}{ds}.
\end{align*}
Consequently:
\begin{align*}
	\dfrac{1}{2}\dfrac{du}{ds} &= g(\sin{\theta(s)} - C_{rr} \cos{\theta(s)})  - \dfrac{\rho C_D A u}{2m}.    
\end{align*}
Multiplying both sides by 2 yields:
\begin{align*}
	\dfrac{du}{ds} &= 2g(\sin{\theta(s)} - C_{rr} \cos{\theta(s)})  - \dfrac{\rho C_D A u}{m}.
\end{align*}
Adding the $u$ term to the left-hand side:
\begin{align}
    \dfrac{du}{ds} + \dfrac{\rho C_D A u}{m} &= 2g(\sin{\theta(s)} - C_{rr} \cos{\theta(s)}).\label{usc}
\end{align}