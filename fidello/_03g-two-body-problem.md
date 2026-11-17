## Example: gravitational two-body problem
~~~
<figure id="two-body" style="float: right;">
<img src="/fidello/Gravitational_two_body_problem.svg" width="400px" />
<figcaption>A system consisting of two bodies in gravitational motion around each other.</figcaption>
</figure>
~~~
In this example, we are modelling the motion of a moving mass around a mass centred at the origin. The relationship between these polar coordinates and our Cartesian coordinates are, for our moving mass:
\begin{align*}
x &= r\cos{\theta} \\
y &= r\sin{\theta}.
\end{align*}
Their time derivatives are:
\begin{align*}
\dot{x} &= \dot{r} \cos{\theta} - r\dot{\theta}\sin{\theta} \\
\dot{y} &= \dot{r}\sin{\theta} + r\dot{\theta}\cos{\theta}.
\end{align*}
Hence velocity squared is:
\begin{align*}
v^2 &= \dot{x}^2 + \dot{y}^2 \\
&= (\dot{r} \cos{\theta} - r\dot{\theta}\sin{\theta})^2 + (\dot{r}\sin{\theta} + r\dot{\theta}\cos{\theta})^2
\end{align*}
Using the binomial theorem, we expand these two brackets to:
\begin{align*}
v^2 &= \dot{r}^2 \cos^2{\theta} + r^2 \dot{\theta}^2 \sin^2{\theta} - 2r\dot{r}\dot{\theta}\cos{\theta}\sin{\theta} + \dot{r}^2 \sin^2{\theta} + r^2 \dot{\theta}^2 \cos^2{\theta} + 2r\dot{r} \dot{\theta}\sin{\theta}\cos{\theta}
\end{align*}
The $2r\dot{r}\dot{\theta}\sin{\theta}\cos{\theta}$ terms have opposite signs and hence cancel each other out:
\begin{align*}
v^2 &= \dot{r}^2 (\cos^2{\theta}+\sin^2{\theta}) + r^2 \dot{\theta}^2 (\sin^2{\theta}+\cos^2{\theta}) \\
&= \dot{r}^2 + r^2 \dot{\theta}^2.
\end{align*}
This gives the kinetic energy:
\begin{align*}
T &= \dfrac{m}{2} \left(\dot{r}^2 + r^2 \dot{\theta}^2\right).
\end{align*}
Where $m$ is the mass of the moving mass. If $M$ is the mass of the central mass and $G$ is the gravitational constant:
\begin{align*}
V &= -\dfrac{GMm}{r}.
\end{align*}
Hence the Lagrangian is:
\begin{align*}
L &= T - V \\
&= \dfrac{m}{2}\left(\dot{r}^2 + r^2 \dot{\theta}^2\right) + \dfrac{GMm}{r}.
\end{align*}
The momentum canonical to $\theta$ is:
\begin{align*}
p_{\theta} &= \dfrac{\partial L}{\partial \dot{\theta}} \\
&= mr^2 \dot{\theta}.
\end{align*}
As for $\dfrac{\partial L}{\partial \theta}$ it is:
\begin{align*}
\dfrac{\partial L}{\partial \theta} &= 0.
\end{align*}
Hence the Euler-Lagrange equation for $\theta$ is:
\begin{align*}
\dfrac{d}{dt}\left(mr^2\dot{\theta}\right) &= 0.
\end{align*}
Which implies that $mr^2 \dot{\theta}$, which is also called the angular momentum, is constant. We will denote this by $p_{\theta}$. Hence:
\begin{align}
\dot{\theta} &= \dfrac{p_{\theta}}{mr^2}. \label{thetadottwobody}
\end{align}
Hence the momentum canonical to $r$ is:
\begin{align*}
p_r &= \dfrac{\partial L}{\partial \dot{r}} \\
&= m\dot{r}.
\end{align*}
As for $\dfrac{\partial L}{\partial r}$ is:
\begin{align*}
\dfrac{\partial L}{\partial r} &= mr\dot{\theta}^2 - \dfrac{GMm}{r^2}.
\end{align*}
Hence the Euler-Lagrange equation for $r$ is:
\begin{align*}
m\ddot{r} - mr\dot{\theta}^2 + \dfrac{GMm}{r^2} &= 0.
\end{align*}
Dividing by $m$ yields:
\begin{align}
\ddot{r} - r\dot{\theta}^2 + \dfrac{GM}{r^2} &= 0.\label{twobodyeq}
\end{align}
Substituting in Equation \eqref{thetadottwobody}:
\begin{align*}
\ddot{r} - r\dfrac{p_{\theta}^2}{m^2r^4} + \dfrac{GM}{r^2} &= 0 \\
\ddot{r} - \dfrac{p_{\theta}^2}{m^2r^3} + \dfrac{GM}{r^2} &= 0.
\end{align*}
Let $u=\dfrac{1}{r}$. Then $r=\dfrac{1}{u}$ and $\dot{r}=-\dfrac{\dot{u}}{u^2}$. We are going to focus on geometrically solving this problem, so we will convert our differential equation from being with respect to $t$ to being with respect to $\theta$. Hence:

\begin{align*}
\dot{u} &= \dfrac{du}{dt} \\
&= \dfrac{du}{d\theta} \dfrac{d\theta}{dt}.
\end{align*}
We will denote $\dfrac{du}{d\theta}$ as $u'$. Using Equation \eqref{thetadottwobody} we get:
\begin{align*}
\dot{u} &= u' \dfrac{p_{\theta}}{mr^2} \\
&= \dfrac{p_{\theta} u'u^2}{m}.
\end{align*}
Hence:
\begin{align*}
\dot{r} &= -\dfrac{p_{\theta} u'u^2}{mu^2} \\
&= -\dfrac{p_{\theta}u'}{m}.
\end{align*}
Differentiating with respect to $t$ again yields:
\begin{align*}
\ddot{r} &= -\dfrac{p_{\theta}\dot{u'}}{m}.
\end{align*}
Hence:
\begin{align*}
\dot{u'} &= \dfrac{d(u')}{dt} \\
&= \dfrac{d(u')}{d\theta} \dfrac{d\theta}{dt} \\
&= u'' \dfrac{p_{\theta}}{mr^2} \\
&= \dfrac{p_{\theta} u'' u^2}{m}.
\end{align*}
Hence:
\begin{align*}
\ddot{r} &= -\dfrac{p_{\theta}^2 u''u^2}{m^2}.
\end{align*}
Substituting this into Equation \eqref{twobodyeq} yields:
\begin{align*}
-\dfrac{p_{\theta}^2 u''u^2}{m^2} - \dfrac{p_{\theta}^2 u^3}{m^2} + GMu^2 &= 0
\end{align*}.
Dividing both sides by $u^2$ yields:
\begin{align*}
-\dfrac{p_{\theta}^2 u''}{m^2} - \dfrac{p_{\theta}^2 u}{m^2} + GM &= 0
\end{align*}
Multiplying by $-\dfrac{m^2}{p_{\theta}^2}$ yields:
\begin{align*}
u'' + u - \dfrac{GMm^2}{p_{\theta}^2} &= 0.
\end{align*}
The solution to this ordinary differential equation is:
\begin{align*}
u(\theta) &= C_1 \cos{\theta} + C_2 \sin{\theta} + \dfrac{GMm^2}{p_{\theta}^2}.
\end{align*}
Let $\mu=\dfrac{GMm^2}{p_{\theta}^2}$. Say that we can rewrite our constants $C_1$ and $C_2$ as:
\begin{align*}
C_1 &= e\mu \cos{\theta_0} \\
C_2 &= e\mu \sin{\theta_0}.
\end{align*}
Where $e$ and $\theta_0$ are new constants we have introduced and will explain once its relevance in $r(\theta)$ becomes clearer. Therefore:
\begin{align*}
u(\theta) &= e\mu\cos{\theta_0}\cos{\theta} + e\mu\sin{\theta_0}\sin{\theta} + \mu \\
&= \mu\left(1+e\left[\cos{\theta_0}\cos{\theta} + \sin{\theta_0}\sin{\theta}\right]\right) \\
&= \mu\left(1+e\cos{(\theta-\theta_0)}\right).
\end{align*}
Hence:
\begin{align*}
r(\theta) &= \dfrac{1}{\mu\left(1+e\cos{(\theta-\theta_0)}\right)}.
\end{align*}
This is the equation for a conic section (parabola, hyperbola or ellipse) with $e$ as its eccentricity and $\theta_0$ is the angle at which the moving mass is closest to the central mass. If $0 \leq e < 1$, it is an ellipse ($e=0$ makes it a perfect circle); if $e=1$, it is a parabola; and if $e>1$, it is a hyperbola. 

So this establishes to us that Lagrangian mechanics predicts that an object in orbit around a larger, central mass, will move in a conic section orbit. For Earth in orbit around the Sun, this conic section has an $e$ value that is almost $0$, which makes it an ellipse. Similarly, for the Moon's orbit around the Earth, $e$ is also almost $0$. For all the planets in the solar system, their $e$ values are between $0$ and $1$. 