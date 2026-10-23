### Simplification to first order
$\ddot{s} = \dot{s} \dfrac{d\dot{s}}{ds}$ hence Equation \eqref{eq2} becomes:
\begin{align*}
	\dot{s} \dfrac{d\dot{s}}{ds} &= \dfrac{GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right] - \dfrac{\rho C_D A\dot{s}^2}{2m}.
\end{align*}
Letting $u=\dot{s}^2$, then $\dfrac{du}{ds} = 2\dot{s}\dfrac{d\dot{s}}{ds}$ and so:
\begin{align*}
	\dfrac{1}{2}\dfrac{du}{ds} &= \dfrac{GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right] - \dfrac{\rho C_D A u}{2m}.
\end{align*}
Multiplying both sides by 2 yields:
\begin{align*}
	\dfrac{du}{ds} &= \dfrac{2GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right] - \dfrac{\rho C_D A u}{m}.
\end{align*}
Adding the $u$ term to the left-hand side:
\begin{align}
	\dfrac{du}{ds} + \dfrac{\rho C_D A}{m}u &= \dfrac{2GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right].\label{dudsnc}
\end{align}