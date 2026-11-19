### Gravitational three-body problem
~~~
<figure id="threebody" style="float: right;">
<img src="/fidello/Gravitational_three_body_problem.svg" width="500px"/>
<figcaption>Gravitational three-body problem.</figcaption>
</figure>
~~~
I thought a good example for you would be the gravitational three-body problem. As like the double pendulum problem, it exhibits chaos. Further, it is the premise of a popular scifi Netflix series *3 Body Problem*, which is based on a novel called *Remembrance of Earth's Past*. As you can see in Figure \figref{threebody}, our axes are selected such that the mass $M$ is centred at the origin. Consequently, the central mass does not move on our axes. So only the two other masses contribute to our system's kinetic energy:
\begin{align*}
T &= \dfrac{m_1}{2}(\dot{r}_1^2 + r_1^2 \dot{\theta}_1^2) + \dfrac{m_2}{2}(\dot{r}_2^2+r_2^2\dot{\theta}_2^2).
\end{align*}
I omitted the calculation of the velocity of the two moving bodies as they are identical to the velocity calculation for the [two-body problem](#example_gravitational_two-body_problem).

To calculate the potential energy, we will need to calculate $r_{12}$, the distance between the masses $m_1$ and $m_2$. To calculate this, we will start with the $x$ and $y$ coordinates for the two moving masses:
\begin{align*}
x_1 &= r_1\cos{\theta_1}\\
y_1 &= r_1 \sin{\theta_1} \\
x_2 &= r_2 \cos{\theta_2} \\
y_2 &= r_2 \sin{\theta_2}.
\end{align*}
Hence the difference in $x$ is $r_2\cos{\theta_2} - r_1\cos{\theta_1}$ and the difference in $y$ is $r_2\sin{\theta_2}-r_1\sin{\theta_1}$. Hence:
\begin{align*}
r_{12}^2 &= (r_2\cos{\theta_2} - r_1\cos{\theta_1})^2 + (r_2\sin{\theta_2} - r_1\sin{\theta_1})^2 \\
&= r_2^2 \cos^2{\theta_2} + r_1^2 \cos^2{\theta_1} -2r_1 r_2 \cos{\theta_1}\cos{\theta_2} + r_2^2 \sin^2{\theta_2} + r_1^2\sin^2{\theta_1}-2r_1r_2 \sin{\theta_1}\sin{\theta_2} \\
&= r_1^2 + r_2^2 - 2r_1 r_2 \cos{(\theta_2-\theta_1)}.
\end{align*}
Square-rooting yields:
\begin{align*}
r_{12} &= \sqrt{r_1^2 + r_2^2 - 2r_1 r_2 \cos{(\theta_2-\theta_1)}}.
\end{align*}
We deliberately ignored the negative possibility as $r_{12}$ is defined such that it must be non-negative. Hence the gravitational potential energy of our system:
\begin{align*}
V &= -\dfrac{GMm_1}{r_1} - \dfrac{GMm_2}{r_2} - \dfrac{Gm_1m_2}{r_{12}}.
\end{align*}
The first of these terms is the potential energy for the gravitational attraction between the central mass and $m_1$. The second is the potential energy for the gravitational attraction between the central mass and $m_2$. The third is the potential energy for the gravitational attraction between $m_1$ and $m_2$. 

Hence the Lagrangian is:
\begin{align}
L &= T - V \\
&= \dfrac{m_1}{2}(\dot{r}_1^2 + r_1^2 \dot{\theta}_1^2) + \dfrac{m_2}{2}(\dot{r}_2^2+r_2^2\dot{\theta}_2^2) + \dfrac{GMm_1}{r_1} + \dfrac{GMm_2}{r_2} + \dfrac{Gm_1m_2}{r_{12}}.
\end{align}
We will calculate $\dfrac{\partial r_{12}}{\partial q_i}$ as this will come in handy when finding $\dfrac{\partial L}{\partial q_i}$.
\begin{align*}
\dfrac{\partial r_{12}}{\partial r_1} &= \dfrac{2r_1-2r_2\cos{(\theta_2-\theta_1)}}{2\sqrt{r_1^2 + r_2^2 - 2r_1 r_2 \cos{(\theta_2-\theta_1)}}} \\
&= \dfrac{r_1-r_2\cos{(\theta_2-\theta_1)}}{\sqrt{r_1^2 + r_2^2 - 2r_1 r_2 \cos{(\theta_2-\theta_1)}}} \\
&= \dfrac{r_1-r_2\cos{(\theta_2-\theta_1)}}{r_{12}}\\
\dfrac{\partial r_{12}}{\partial r_2} &= \dfrac{r_2-r_1\cos{(\theta_2-\theta_1)}}{r_{12}}\\
\dfrac{\partial r_{12}}{\partial \theta_1} &= \dfrac{-2r_1r_2\sin{(\theta_2-\theta_1)}}{2r_{12}}\\
&= \dfrac{-r_1r_2\sin{(\theta_2-\theta_1)}}{r_{12}} \\
\dfrac{\partial r_{12}}{\partial \theta_2} &= \dfrac{2r_1r_2\sin{(\theta_2-\theta_1)}}{2r_{12}}\\
&= \dfrac{r_1r_2\sin{(\theta_2-\theta_1)}}{r_{12}}.
\end{align*}
The canonical momentums are:
\begin{align*}
p_{r_1} &= \dfrac{\partial L}{\partial \dot{r}_1} \\
&= m_1 \dot{r}_1 \\
p_{r_2} &= \dfrac{\partial L}{\partial \dot{r}_2} \\
&= m_2 \dot{r}_2 \\
p_{\theta_1} &= \dfrac{\partial L}{\partial \dot{\theta}_1} \\
&= m_1 r_1^2\dot{\theta}_1 \\
p_{\theta_2} &= \dfrac{\partial L}{\partial \dot{\theta}_2} \\
&= m_2 r_2^2 \dot{\theta}_2.
\end{align*}
As for the spatial derivatives of $L$:
\begin{align*}
\dfrac{\partial L}{\partial r_1} &= m_1 r_1 \dot{\theta}_1^2 - \dfrac{GMm_1}{r_1^2}  - \dfrac{Gm_1 m_2}{r_{12}^2} \cdot \dfrac{\partial r_{12}}{\partial r_1} \\
&= m_1 r_1 \dot{\theta}_1^2 - \dfrac{GMm_1}{r_1^2}- \dfrac{Gm_1 m_2( r_1-r_2 \cos{(\theta_2-\theta_1)})}{r_{12}^3}\\
\dfrac{\partial L}{\partial r_2} &= m_2 r_2 \dot{\theta}_2^2 - \dfrac{GMm_2}{r_2^2} - \dfrac{Gm_1 m_2}{r_{12}^2} \cdot \dfrac{\partial r_{12}}{\partial r_2} \\
&= m_2 r_2 \dot{\theta}_2^2 - \dfrac{GMm_2}{r_2^2} - \dfrac{Gm_1 m_2( r_2-r_1 \cos{(\theta_2-\theta_1)})}{r_{12}^3}\\
\dfrac{\partial L}{\partial \theta_1} &= - \dfrac{Gm_1 m_2}{r_{12}^2} \cdot -\dfrac{\partial r_{12}}{\partial \theta_1} \\
&= -\dfrac{Gm_1 m_2 r_1r_2\sin{(\theta_2-\theta_1)}}{r_{12}^3}\\
\dfrac{\partial L}{\partial \theta_2} &= - \dfrac{Gm_1 m_2}{r_{12}^2} \cdot \dfrac{\partial r_{12}}{\partial \theta_2} \\
&= -\dfrac{Gm_1 m_2r_1r_2\sin{(\theta_2-\theta_1)}}{r_{12}^3}.
\end{align*}
Hence the $r_1$ Euler-Lagrange equation is:
\begin{align*}
m_1 \ddot{r}_1 &= m_1 r_1 \dot{\theta}_1^2 - \dfrac{GMm_1}{r_1^2} - \dfrac{Gm_1 m_2( r_1-r_2 \cos{(\theta_2-\theta_1)})}{r_{12}^3} \\
\ddot{r}_1 &= r_1 \dot{\theta}_1^2 - \dfrac{GM}{r_1^2} - \dfrac{Gm_2( r_1-r_2 \cos{(\theta_2-\theta_1)})}{r_{12}^3}.
\end{align*}
Hence the $r_2$ Euler-Lagrange equation is:
\begin{align*}
m_2 \ddot{r}_2 &= m_2 r_2 \dot{\theta}_2^2 - \dfrac{GMm_2}{r_2^2} - \dfrac{Gm_1 m_2( r_2-r_1 \cos{(\theta_2-\theta_1)})}{r_{12}^3} \\
\ddot{r}_2 &= r_2 \dot{\theta}_2^2 - \dfrac{GM}{r_2^2} - \dfrac{Gm_1( r_2-r_1 \cos{(\theta_2-\theta_1)})}{r_{12}^3}.
\end{align*}
As for the $\theta_1$ Euler-Lagrange equation, it is:
\begin{align*}
2m_1 r_1\dot{r}_1 \dot{\theta}_1 + m_1 r_1^2 \ddot{\theta}_1 &= \dfrac{Gm_1 m_2 r_1r_2\sin{(\theta_2-\theta_1)}}{r_{12}^3}.
\end{align*}
Divide both sides by $m_1r_1^2$:
\begin{align*}
\dfrac{2\dot{r}_1 \dot{\theta}_1}{r_1} + \ddot{\theta}_1 &= \dfrac{Gm_2 r_2\sin{(\theta_2-\theta_1)}}{r_1 r_{12}^3} \\
\ddot{\theta}_1 &= -\dfrac{2\dot{r}_1 \dot{\theta}_1}{r_1} + \dfrac{Gm_2 r_2\sin{(\theta_2-\theta_1)}}{r_1 r_{12}^3}.
\end{align*}
As for $\theta_2$:
\begin{align*}
2m_2 r_2\dot{r}_2 \dot{\theta}_2 + m_2 r_2^2 \ddot{\theta}_2 &= \dfrac{Gm_1 m_2 r_1r_2\sin{(\theta_2-\theta_1)}}{r_{12}^3}.
\end{align*}
Dividing both sides by $m_2r_2^2$:
\begin{align*}
\dfrac{2\dot{r}_2 \dot{\theta}_2}{r_2} + \ddot{\theta}_2 &= -\dfrac{Gm_1 r_1\sin{(\theta_2-\theta_1)}}{r_2 r_{12}^3} \\
\ddot{\theta}_2 &= -\dfrac{2\dot{r}_2 \dot{\theta}_2}{r_2} -\dfrac{Gm_1 r_1\sin{(\theta_2-\theta_1)}}{r_2 r_{12}^3}.
\end{align*}

This system of ordinary differential equations cannot be analytically solved, unlike the equations we obtained for the two-body problem. Consequently, we are forced to resort to numerical methods to solve this system. Below is Julia code to numerically solve the system for a made-up example of the three-body problem.

```julia:./code/threebody
\juliacode{fidello/code/3body.jl}
```

The above Julia code produces the animation in Figure \figref{threebodyan}. As you can see, the orbits majorly and seemingly unpredictably change over time. This is pretty typical in the three-body problem: either the orbits are fairly stable and unchanging, or they are unstable and prone to objects being ejected from the system or colliding with one another. 

~~~
<figure id="threebodyan" style="float: right; width: 100%">
<img src="/three_body_system.gif" width="800px"/>
<figcaption>A three-body system simulation.</figcaption>
</figure>
~~~