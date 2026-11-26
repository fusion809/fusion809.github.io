### Newton-Kantorovich method
When we have a nonlinear differential equation and we apply a finite difference method, finite element method, spectral element method or spectral method to solve it, we obtain a nonlinear system of equations. We could solve this using Newton's method, or alternatively we could first linearize our differential equation to avoid the need to use Newton's method. Linearizing differential equations can be achieved by essentially using the first two terms (the constant term and linear term) of the differential equation's Taylor series. Specifically, we expand around the our $n$th approximation to the solution to our differential equation. As an example, consider Equation \eqref{simppen}. We will later encounter this differential equation when we study the problem of the simple pendulum. We would define:
\begin{align*}
\theta_{n+1} &\def \theta_{n} + \Delta.
\end{align*}
Where $\Delta$ will be the subject of our linearized differential equation and $\theta_n$ is our $n$th approximation of $\theta$. The Newton-Kantorovich method would yield:
\begin{align*}
\ddot{\theta}_n + \ddot{\Delta} &= -\dfrac{g}{l}\left[\cos{\theta_{n}} + \dfrac{d}{d\theta} \cos{\theta} \Big|_{\theta=\theta_n} \Delta\right] \\
\ddot{\theta}_n + \ddot{\Delta} &= -\dfrac{g}{l}\cos{\theta_{n}} + \dfrac{g}{l}\sin{\theta_{n}}\Delta \\
\ddot{\Delta} - \dfrac{g}{l}\sin{\theta_{n}}\Delta &= -\ddot{\theta}_n - \dfrac{g}{l}\cos{\theta_{n}}.
\end{align*}
Say we had the initial conditions $\theta(0)=\dot{\theta}(0)=0$. Then our initial conditions for $\Delta$ would be:
\begin{align*}
\Delta(0) &= -\theta_n(0) \\
\dot{\Delta}(0) &= -\dot{\theta}_n(0).
\end{align*}
Why? Because then $\theta_{n+1}$ should automatically satisfy the initial conditions. An important question to answer is what should be our first guess as to the solution? In the Newton-Kantorovich method, the first guess matters a lot. A poor first guess can lead to the method failing to converge, which means we will not actually get a solution. We know, from our solution in the [Without air resistance](#without_air_resistance) section, that the simple pendulum under these initial conditions goes to a minimum angle of $-\pi$ and maximum angle of $0$. The curve looks roughly like a cosine wave, too. If the period of the solution is $T$, a decent first guess might be:
\begin{align*}
\theta(t) &= -\dfrac{\pi}{2}\left[1-\cos{\left(\dfrac{2\pi t}{T}\right)}\right].
\end{align*}
How would we determine $T$? We can use Equation \eqref{tsimppen}. Ignoring $\beta$ (which is a product of the indefinite nature of the integral), making our integral definite from $-\pi$ to $0$ and assuming $t$ must be positive gives:
\begin{align*}
t &= \int_{-\pi}^0 \sqrt{\dfrac{ml^2}{2(E-mgl\sin{\theta})}} d\theta
\end{align*}
This should be half our period, as a full period should take us back to $0$, but we cannot integrate from $0$ to $0$ as when the upper and lower bounds are the same, integrals return $0$. $\theta(0)=0$ and hence $V=mgl\sin{\theta}$ should also be zero at $t=0$. $\dot{\theta}(0)=0$ would also imply that $T=\dfrac{ml^2\dot{\theta}^2}{2}$ is zero. Hence the total energy of the system, $E$, is $0$. Hence:
\begin{align*}
t &= \int_{-\pi}^0 \sqrt{\dfrac{ml^2}{2(-mgl\sin{\theta})}} d\theta \\
&= \int_{-\pi}^0 \sqrt{\dfrac{l}{2(-g\sin{\theta})}} d\theta \\
&= \sqrt{\dfrac{l}{2g}} \int_{-\pi}^0 \dfrac{1}{\sqrt{-\sin{\theta}}} d\theta.
\end{align*}
Hence the period should be:
\begin{align*}
T &= 2 \sqrt{\dfrac{l}{2g}} \int_{-\pi}^0 \dfrac{1}{\sqrt{-\sin{\theta}}} d\theta \\
&= \sqrt{\dfrac{2l}{g}}\int_{-\pi}^0 \dfrac{1}{\sqrt{-\sin{\theta}}} d\theta.
\end{align*}

Code Listing \coderef{newtonkantorovich} displays the code we used to implement this solution on $t\in[0, T]$. 
\codelisting{newtonkantorovich}{Script for solving simple pendulum problem with the Newton-Kantorovich method.}{./code/newtonkantorovich}{fidello/code/newton-kantorovich.jl}
\output{./code/newtonkantorovich}
\figure{newtonkantorovichplot}{right}{600px}{/pendulum_theta.svg}{Newton-Kantorovich solution to the simple pendulum equation.}