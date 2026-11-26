#### Without air resistance
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

Technically speaking it is possible to analytically solve Equation \eqref{simppen}, but it requires the use of non-standard functions called elliptic integrals. Another option is we can use the Taylor series around $a=-\dfrac{\pi}{2}$. Why $-\dfrac{\pi}{2}$? Because this corresponds to the pendulum being on the negative $y$ axis. This is an equilibrium point of the simple pendulum, as if the pendulum has no kinetic energy (in other words $\dot{\theta}=0$) when it is at this point, it will stay perfectly still and remain at $\theta=-\dfrac{\pi}{2}$ forever. So what we want to study with the Taylor series is what happens if we only slightly perturb the pendulum from this equilibrium point, as a slight perturbation from this point should keep $\cos{\theta}$ well approximated by the first couple of terms of its Taylor series. Okay, so we already know the $n$th derivative of $\cos{\theta}$ with respect to $\theta$ from the [Maclaurin and Taylor series](#maclaurin_and_taylor_series) section, it is:
\begin{align*}
\dfrac{d^n}{d\theta^n} \cos{\theta} &= \begin{cases}
(-1)^{\lfloor \dfrac{n+1}{2}\rfloor} \cos{\theta}, \hspace{0.2cm}n\hspace{0.1cm}\mathrm{even}\\
(-1)^{\lfloor \dfrac{n+1}{2}\rfloor} \sin{\theta}, \hspace{0.2cm}n\hspace{0.1cm}\mathrm{odd}
\end{cases}.
\end{align*}
If we replace even $n$ values with $2n$ (where $n$ here is simply a non-negative integer) and odd $n$ with $2n+1$:
\begin{align*}
\dfrac{d^{2n}}{d\theta^{2n}} \cos{\theta} &= (-1)^{n} \cos{\theta} \\
\dfrac{d^{2n+1}}{d\theta^{2n+1}} \cos{\theta} &= (-1)^{n+1} \sin{\theta}.
\end{align*}
Substituting in $\theta=a$ which is $-\dfrac{\pi}{2}$:
\begin{align*}
\dfrac{d^{2n}}{d\theta^{2n}} \cos{\theta} \Big|_{\theta=-\dfrac{\pi}{2}} &= (-1)^{n} \cdot 0 \\
&= 0 \\
\dfrac{d^{2n+1}}{d\theta^{2n+1}} \cos{\theta} \Big|_{\theta=-\dfrac{\pi}{2}} &= (-1)^{n+1} \cdot -1 \\
&= (-1)^{n+2} \\
&= (-1)^n.
\end{align*}
Hence our Taylor series is:
\begin{align*}
\cos{\theta} &= \sum_{n=0}^{\infty} \dfrac{(-1)^{n} \left(\theta-(-\dfrac{\pi}{2})\right)^{2n+1}}{(2n+1)!} \\
&= \sum_{n=0}^{\infty} \dfrac{(-1)^{n} \left(\theta+\dfrac{\pi}{2}\right)^{2n+1}}{(2n+1)!}.
\end{align*}
Now we want to make our differential equation as simple as possible, so we will only consider the first term in this Taylor series, the $n=0$ term. In which case $\cos{\theta} \approx \theta + \dfrac{\pi}{2}$. Hence Equation \eqref{simppen} becomes:
\begin{align}
\ddot{\theta} &=-\dfrac{g}{l}\theta - \dfrac{g\pi}{2l}.\label{simppenlin}
\end{align}
What function if double differentiated yields the original function back with $-\dfrac{g}{l}$ coefficient and a constant term? The answer is:
\begin{align*}
\theta(t) &= C_1 \cos{\left(\sqrt{\dfrac{g}{l}}t\right)} + C_2 \sin{\left(\sqrt{\dfrac{g}{l}}t\right)} + C_3.
\end{align*}
Where $C_1$ and $C_2$ are constants that depend on our initial conditions (similarly to the arbitrary constants of integration you get when using indefinite integrals) and $C_3$ is a constant we will work out next, right after differentiating $\theta(t)$ for later use:
\begin{align*}
\dot{\theta}(t) &= -\sqrt{\dfrac{g}{l}}C_1\sin{\left(\sqrt{\dfrac{g}{l}}t\right)} + \sqrt{\dfrac{g}{l}}C_2 \cos{\left(\sqrt{\dfrac{g}{l}}t\right)}.
\end{align*}
$C_3$ is a constant, so double differentiating it will remove it. Hence we know that $C_3$ only enters Equation \eqref{simppenlin} via the $\theta$ term on the right-hand side. So we need to impose the condition that:
\begin{align*}
-\dfrac{g}{l}C_3 - \dfrac{g\pi}{2l} &= 0
\end{align*}
to ensure that the left-hand side and right-hand sides of Equation \eqref{simppenlin} match. Multiplying by $-\dfrac{l}{g}$ yields:
\begin{align*}
C_3 + \dfrac{\pi}{2} &= 0.
\end{align*}
Subtracting $\dfrac{\pi}{2}$ yields:
\begin{align*}
C_3 &= -\dfrac{\pi}{2}.
\end{align*}
At $t=0$, let $\theta = \theta_0$ and $\dot{\theta}=\dot{\theta}_0$. These values simply represent the initial values of $\theta$ and $\dot{\theta}$.
\begin{align*}
\theta_0 &= C_1 \cos{0} + C_2 \sin{0} - \dfrac{\pi}{2} \\
&= C_1 - \dfrac{\pi}{2}.
\end{align*}
Adding $\dfrac{\pi}{2}$ to both sides yields:
\begin{align*}
C_1 &= \theta_0 + \dfrac{\pi}{2}.
\end{align*}
As for $\dot{\theta}_0$:
\begin{align*}
\dot{\theta}_0 &= -\sqrt{\dfrac{g}{l}}\left(\theta_0 + \dfrac{\pi}{2}\right)\sin{0} + \sqrt{\dfrac{g}{l}}C_2 \cos{0} \\
&=  \sqrt{\dfrac{g}{l}}C_2
\end{align*}
Multiplying both sides $\sqrt{\dfrac{l}{g}}$:
\begin{align*}
C_2 &= \sqrt{\dfrac{l}{g}} \dot{\theta}_0.
\end{align*}
Hence:
\begin{align*}
\theta(t) &= \left(\theta_0 + \dfrac{\pi}{2}\right)\cos{\left(\sqrt{\dfrac{g}{l}}t\right)} + \dot{\theta}_0 \sqrt{\dfrac{l}{g}} \sin{\left(\sqrt{\dfrac{g}{l}}t\right)} -\dfrac{\pi}{2}.
\end{align*}
So we can see $\theta(t)$ is approximated, when it deviates only slightly from $\theta=-\dfrac{\pi}{2}$, quite well by sine and cosine functions. More generally, we can use numerical techniques to approximate our solution. Below is some Julia code that numerically approximates the solution for $\theta_0=\dot{\theta}_0=0$, $g=9.81\mathrm{m}\cdot \mathrm{s}^{-2}$, $l=1\mathrm{m}$ and $m=1\mathrm{kg}$. These initial conditions means the pendulum starts stationary on the positive $x$-axis. Then gravity pulls it down and up towards the negative $x$-axis (at $\theta=-\pi$) and then back down and towards the positive $x$-axis ($\theta=0$) and over and over. 
\codelisting{simppen}{Simple pendulum simulation code.}{./code/simppen}{fidello/code/simppen.jl}
Figure \figref{simppen} shows the plot of $\theta$ vs $t$ produced by Code Listing \coderef{simppen}. Figure \figref{simppenan} contains an animation of the simple pendulum simulated by the above Julia code.

\figure{simppen}{left}{900px}{/Simple_pendulum_theta.svg}{Simple pendulum angle plot.}
\figure{simppenan}{right}{900px}{/Simple_pendulum.gif}{Simple pendulum animation.}
~~~
<p style="clear: both;">If you would like to try some of your own initial conditions and calculate the solution numerically for these conditions, see <a href="/simplePendulum/woErrorAnalysis">this page</a>.</p>
~~~

