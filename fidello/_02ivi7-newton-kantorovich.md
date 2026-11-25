### Newton-Kantorovich method
When we have a nonlinear differential equation and we apply a finite difference method, finite element method, spectral element method or spectral method to solve it, we obtain a nonlinear system of equations. We could solve this using Newton's method, or alternatively we could first linearize our differential equation to avoid the need to use Newton's method. Linearizing differential equations can be achieved by essentially using the first two terms (the constant term and linear term) of the differential equation's Taylor series. Specifically, we expand around the our $n$th approximation to the solution to our differential equation. As an example, consider the differential equation:
\begin{align*}
\ddot{\theta} &= -\dfrac{g}{l}\cos{\theta}.
\end{align*}
We will later encounter this differential equation when we study the problem of the simple pendulum. We would define:
\begin{align*}
\theta_{n+1} &= \theta_{n} + \Delta.
\end{align*}
The Newton-Kantorovich method would yield:
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
An important consideration is what should be our first guess as to the solution. Well, we know that simple pendulum under these initial conditions goes to a minimum angle of $-\pi$ and maximum angle of $0$. The curve looks roughly like a cosine wave, too. If the period of the solution is $T$, a decent first guess might be:
\begin{align*}
\theta(t) &= -\dfrac{\pi}{2}\left[1-\cos{\left(\dfrac{2\pi t}{T}\right)}\right].
\end{align*}
How would we determine $T$? We can use Equation \eqref{tsimppen}. Ignoring $\beta$ (which is a product of the indefinite nature of the integral), making our integral definite from $-\pi$ to $0$ (this interval should give us half the period) and assuming $t$ must be positive gives:
\begin{align*}
t &= \int_{-\pi}^0 \sqrt{\dfrac{ml^2}{2(E-mgl\sin{\theta})}} d\theta
\end{align*}
$\theta(0)=0$ and hence $V=mgl\sin{\theta}$ should also be zero at $t=0$. $\dot{\theta}(0)=0$ would also imply that $T=\dfrac{ml^2\dot{\theta}^2}{2}$ is zero. Hence the total energy of the system, $E$, is $0$. Hence:
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
\figure{newtonkantorovichplot}{right}{600px}{/pendulum_theta.svg}{Newton-Kantorovich solution to the simple pendulum equation.}