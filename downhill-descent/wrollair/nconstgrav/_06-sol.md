### Solution
We first seek an integration factor, $I(s)$, so that we can re-write Equation \eqref{dudsnc} as:
\begin{align*}
I(s) \dfrac{du}{ds} + I'(s)\dfrac{du}{ds} &= I(s)\dfrac{2GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right].
\end{align*}
Hence $I'(s)$ must equal $I(s)\dfrac{\rho C_D A}{m}$, which therefore implies that 
\begin{align*}
I(s)&=e^{\int_0^s \frac{\rho C_D A}{m}ds} \\
&=e^{\frac{\rho C_D As}{m}}.
\end{align*}
We can replace the left-hand side with $\dfrac{d}{ds}\left[I(s)u(s)\right]$ hence:
\begin{align*}
\dfrac{d}{ds}\left[I(s)u(s)\right] &= I(s)\dfrac{2GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right].\label{ddsIunc}
\end{align*}
Therefore Equation \eqref{ddsIunc} becomes:
\begin{align}
\dfrac{d}{ds}\left[e^{\frac{\rho C_D As}{m}}u(s)\right] &= 2GMe^{\frac{\rho C_D As}{m}}\dfrac{\sin{\theta(s)}-C_{rr}\cos{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2}.\label{ddsIfunc}
\end{align}
Integrating both sides of Equation \eqref{ddsIfunc} from 0 to $s$ with respect to $s$ (calling it a new dummy variable $\sigma$; $\sigma$ in our denominator will be renamed $\xi$) yields:
\begin{align}
e^{\frac{\rho C_D As}{m}}u(s)-e^0 u(0) &= 2GM\int_0^s e^{\frac{\rho C_D A\sigma}{m}}\dfrac{\sin{\theta(\sigma)}-C_{rr}\cos{\theta(\sigma)}}{\left(r_{\mathrm{summit}} - \int_0^{\sigma} \sin{\theta(\xi)}d\xi\right)^2} d\sigma.
\end{align}
The trike's initial velocity is assumed to be zero, hence $u(0)=0$:
\begin{align*}
e^{\frac{\rho C_D As}{m}}u(s) &= 2GM\int_0^s e^{\frac{\rho C_D A\sigma}{m}}\dfrac{\sin{\theta(\sigma)}-C_{rr}\cos{\theta(\sigma)}}{\left(r_{\mathrm{summit}} - \int_0^{\sigma} \sin{\theta(\xi)}d\xi\right)^2} d\sigma.
\end{align*}
Multiplying by $e^{-\frac{\rho C_D As}{m}}$ yields:
\begin{align*}
u(s) &= 2GMe^{-\frac{\rho C_D As}{m}}\int_0^s e^{\frac{\rho C_D A\sigma}{m}}\dfrac{\sin{\theta(\sigma)}-C_{rr}\cos{\theta(\sigma)}}{\left(r_{\mathrm{summit}} - \int_0^{\sigma} \sin{\theta(\xi)}d\xi\right)^2} d\sigma.
\end{align*}
Square rooting $u$ to obtain velocity yields:
\begin{align}
    v(s) &= \sqrt{2GMe^{-\frac{\rho C_D As}{m}}\int_0^s e^{\frac{\rho C_D A\sigma}{m}}\dfrac{\sin{\theta(\sigma)}-C_{rr} \cos{\theta(\sigma)}}{\left(r_{\mathrm{summit}}-\int_0^\sigma \sin{\theta(\xi)} d\xi\right)^2} d \sigma}.
\end{align}
Below are the velocities achieved in the above model for various $C_D A$ and $C_{rr}$ values when descending Castle Hill. The simulation did not end at 0 metres above sea level; it ended at about 36 metres above sea level, which is approximately at the high end of Gregory Street. 

~~~
<figure>
<img src="https://fusion809.github.io/images/SageMath/Castle%20hill%20descent%20speed%20non-constant%20gravity.svg" width="100%">
<figcaption>Velocity while descending Castle Hill, assuming non-constant gravity.</figcaption>
</figure>
~~~