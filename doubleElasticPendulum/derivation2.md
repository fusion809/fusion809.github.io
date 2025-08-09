@def hassim=false;
@def title="Deriving the equations of motion for the double elastic pendulum"
@def maxtoclevel = 5
@def mintoclevel = 1
~~~
<figure>
    <img src="/doubleElasticPendulum/Double elastic pendulum.svg" width="500px"></img>
    <figcaption style="font-weight: bold; font-size: 18px;">Figure 1: Diagram of the double elastic pendulum.</b></figcaption>
</figure>
~~~

This article focuses on deriving the equations of motion of the coupled elastic pendulums depicted in Figure 1. If you would like to see their solution, they can be found [here](/doubleElasticPendulum/).

We will derive our own equations of motion based on the Euler-Lagrange equation with dissipative forces

\begin{align}
\dfrac{d}{dt} \dfrac{\partial \mathcal{L}}{\partial \dot{q}_i} - \dfrac{\partial \mathcal{L}}{\partial q_i} &= \sum_j \vec{F}_{D,j} \cdot \dfrac{\partial \vec{r}_j}{\partial q_i}.\label{ELD}
\end{align}
        
Where
* $j$ refers to the component of the system we are analysing.
* $(q_i)$ are the generalized coordinates of the system.
* $(\dot{q}_i)$ are the first time derivatives of the generalized coordinates of the system.
* $\vec{r}_j$ is the position vector of component $j$ of the system.
* $\mathcal{L}$ is the Lagrangian &mdash; the difference between the kinetic and potential energy &mdash; of the system.
* $p_i = \dfrac{\partial \mathcal{L}}{\partial \dot{q}_i}$ is the generalized momentum canonical to $q_i$.
* $F_i = \dfrac{\partial \mathcal{L}}{\partial q_i}$ is the generalized force canonical to $q_i$.
* $\vec{F}_{D,j}$ is the dissipative force vector for component $j$.
* $\hat{e}_{j,i} = \dfrac{\partial \vec{r}_j}{\partial q_i}$ is the generalized basis vector canonical to $q_i$ for component $j$ of the system.
* The left-hand side of Equation Equation \eqref{ELD} can also be represented as $-\dfrac{\delta \mathcal{L}}{\delta q_i}$, where $\dfrac{\delta \mathcal{L}}{\delta q_i}$ is the functional derivative of the Lagrangian with respect to $q_i$. To simplify things, we will call $-\dfrac{\delta \mathcal{L}}{\delta q_i} = \dfrac{\delta' \mathcal{L}}{\delta' q_i}$
* The right-hand side of Equation Equation \eqref{ELD} is also called the generalized dissipative force and can be represented as $Q_i$.

\tableofcontents

# Coordinates, velocities and generalized basis vectors
As can be seen, we have four degrees of freedom in this system. The angles the two pendulums make with the positive $x$-axis &mdash; $\theta_1$ and $\theta_2$, respectively &mdash; are among our degrees of freedom. We will also need degrees of freedom corresponding to the lengths of the pendulum rods. These degrees of freedom could either be the extent to which they are extended beyond their rest length or their total length. For the sake of simplicity, we will opt to use their total lengths &mdash; $r_1$ and $r_2$, respectively. Hence

## Bob 1
\begin{align*}
	x_{1b} &= r_1 \cos{\theta_1} & \dot{x}_{1b} &= \dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1}\\
	y_{1b} &= r_1 \sin{\theta_1} & \dot{y}_{1b} &= \dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1}
\end{align*}

This means that the velocity of the first pendulum bob is

\begin{align*}
	\vec{v}_{1b} &= \begin{bmatrix}
		\dot{x}_{1b} \\
		\dot{y}_{1b}
	\end{bmatrix} \\
	&= \begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1}
	\end{bmatrix}.
\end{align*}

Hence 
\begin{align*}
	|\vec{v}_{1b}|^2 &= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2. 
\end{align*}

## Rod 1
We will calculate the rod positions and velocities two different ways. For calculating the kinetic energy, we will parameterize them in terms of $s$ &mdash; the distance along the pendulum &mdash; so that we can integrate over them. For calculating the potential energy and generalized dissipative force, we will calculate them using a centre of mass approach, as we will assume centre of mass approximations will work for these. 

### Parameterization
\begin{align*}
x_{1r} &= s\cos{\theta_1} & \dot{x}_{1r} = \dot{s}\cos{\theta_1} - s\dot{\theta}_1 \sin{\theta_1} \\
y_{1r} &= s\sin{\theta_1} & \dot{y}_{1r} = \dot{s}\sin{\theta_1} + s\dot{\theta}_1 \cos{\theta_1} \\
|\vec{v}_{1r}|^2 &= \dot{s}^2 + s^2\dot{\theta}_1^2 \\
&= \dfrac{s^2\dot{r}_1^2}{r_1^2} + s^2\dot{\theta}_1^2 \\
&= s^2 \left[\dfrac{\dot{r}_1^2}{r_1^2}+\dot{\theta}_1^2\right].
\end{align*}

Here we have assumed that the motion of the pendulum rod is uniform. 

### Centre of mass
\begin{align*}
x_{1r} &= \dfrac{r_1\cos{\theta_1}}{2} & \dot{x}_{1r} = \dfrac{\dot{r_1}\cos{\theta_1}}{2} - \dfrac{r_1\dot{\theta}_1 \sin{\theta_1}}{2} \\
y_{1r} &= \dfrac{r_1\sin{\theta_1}}{2} & \dot{y}_{1r} = \dfrac{\dot{r_1}\sin{\theta_1}}{2} + \dfrac{r_1\dot{\theta}_1 \cos{\theta_1}}{2} \\
\vec{v}_{1r} &= \dfrac{1}{2} \begin{bmatrix}
\dot{r}_1 \cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
\dot{r}_1 \sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}
\end{bmatrix} \\
|\vec{v}_{1r}|^2 &= \dfrac{1}{4} \left(\dot{r}_1^2 + r_1^2 \dot{\theta}_1^2\right)\\
&= \dfrac{|\vec{v}_{1b}|^2}{4}.
\end{align*}

## Bob 2
\begin{align*}
	x_{2b} &= x_{1b} + r_2\cos{\theta_2} & \dot{x}_2 &= \dot{x}_{1b} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
	y_{2b} &= y_{1b} + r_2\sin{\theta_2} & \dot{y}_2 &= \dot{y}_{1b} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}.
\end{align*}

As for the velocity of the second pendulum bob, it is
\begin{align*}
	\vec{v}_{2b} &= \begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix}.
\end{align*}

Let $\Delta = \theta_2-\theta_1$, then square of the velocity is

\begin{align*}
	|\vec{v}_{2b}|^2 &= \dot{r}_1^2 \cos^2{\theta_1} + r_1^2 \dot{\theta}_1^2\sin^2{\theta_1} + \dot{r}_2^2\cos^2{\theta_2} + r_2^2\dot{\theta}_2^2\sin^2{\theta_2} -2r_1\dot{r}_1\dot{\theta}_1 \cos{\theta}_1\sin{\theta_1} + 2\dot{r}_1\dot{r}_2\cos{\theta_1}\cos{\theta_2} - 2\dot{r}_1r_2\dot{\theta}_2\cos{\theta_1}\sin{\theta_2} - 2r_1\dot{r}_2 \dot{\theta}_1 \sin{\theta_1}\cos{\theta_2} + 2r_1r_2 \dot{\theta}_1\dot{\theta}_2\sin{\theta_1}\sin{\theta_2}\\
	&-2\dot{r}_2r_2\dot{\theta}_2\cos{\theta_2}\sin{\theta_2} + \dot{r}_1^2\sin^2{\theta_1} + r_1^2\dot{\theta}_1^2\cos^2{\theta_1} + \dot{r}_2^2\sin^2{\theta_2} + r_2^2\dot{\theta}_2^2\cos^2{\theta_2} + 2r_1\dot{r}_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + 2\dot{r}_1\dot{r}_2\sin{\theta_1}\sin{\theta_2} + 2\dot{r}_1r_2\dot{\theta}_2\sin{\theta_1}\cos{\theta_2} + 2r_1\dot{r}_2\dot{\theta}_1 \cos{\theta_1}\sin{\theta_2}\\
	&+2r_1r_2\dot{\theta}_1\dot{\theta}_2\cos{\theta_1}\cos{\theta_2} + 2r_2\dot{r}_2\dot{\theta}_2\sin{\theta_2}\cos{\theta_2} \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2r_1\dot{r}_1\dot{\theta}_1(-\cos{\theta_1}\sin{\theta_1} + \cos{\theta_1}\sin{\theta_1}) + 2\dot{r}_1\dot{r}_2(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2})+2\dot{r}_1r_2\dot{\theta}_2(-\cos{\theta_1}\sin{\theta_2} + \sin{\theta_1}\cos{\theta_2}) \\
	&+ 2r_1\dot{r}_2\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_2}+\cos{\theta_1}\sin{\theta_2}) + 2r_1r_2\dot{\theta}_1\dot{\theta}_2(\sin{\theta_1}\sin{\theta_2} + \cos{\theta_1}\cos{\theta_2}) + 2r_2\dot{r}_2\dot{\theta}_2 (-\cos{\theta_2}\sin{\theta_2} + \sin{\theta_2}\cos{\theta_2}) \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2\dot{r}_1\dot{r}_2 \cos{(\theta_2-\theta_1)} - 2\dot{r}_1r_2\dot{\theta_2}\sin{(\theta_2-\theta_1)} + 2r_1\dot{r}_2 \dot{\theta}_1 \sin{(\theta_2-\theta_1)} + 2r_1r_2\dot{\theta}_1\dot{\theta}_2\cos{(\theta_2-\theta_1)} \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2\dot{r}_1\dot{r}_2 \cos{\Delta} - 2\dot{r}_1r_2\dot{\theta_2}\sin{\Delta} + 2r_1\dot{r}_2 \dot{\theta}_1 \sin{\Delta} + 2r_1r_2\dot{\theta}_1\dot{\theta}_2\cos{\Delta} \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2\cos{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\sin{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2).
\end{align*}

Let us define $|\vec{v}_{2b}^I|^2 = \dot{r}_2^2+r_2^2\dot{\theta}_2^2$, as this will simplify our Lagrangian later. As for the remaining terms in $|\vec{v}_{2b}|^2$, they are contained with $|\Delta \vec{v}_{2,1}|^2$.

\begin{align*}
	|\Delta \vec{v}_{2,1}|^2 &= 2\cos{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\sin{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2).
\end{align*}

## Rod 2
### Parameterization
\begin{align*}
x_{2r} &= r_1\cos{\theta_1} + s\cos{\theta_2} & \dot{x}_{2r} = \dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1\sin{\theta_1} + \dot{s}\cos{\theta_2} - s\dot{\theta}_2 \sin{\theta_2} \\
y_{2r} &= r_1\sin{\theta_1} + s\sin{\theta_2} & \dot{y}_{2r} = \dot{r}_1 \sin{\theta_1} + r_1 \dot{\theta}_1\cos{\theta_1} + \dot{s}\sin{\theta_2} + s\dot{\theta}_2 \cos{\theta_2}
\end{align*}
\begin{align*}
|\vec{v}_{2r}|^2 &= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{s}^2 + s^2\dot{\theta}_2^2 + 2\cos{\Delta}\left(\dot{r}_1\dot{s} + r_1s\dot{\theta}_1\dot{\theta}_2\right) + 2\sin{\Delta}\left(r_1\dot{s}\dot{\theta}_1-\dot{r}_1s\dot{\theta}_2\right) \\
&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dfrac{\dot{r}_2^2 s^2}{r_2^2} + s^2\dot{\theta}_2^2 + 2\cos{\Delta}\left(\dfrac{\dot{r}_1\dot{r}_2 s}{r_2} + r_1s\dot{\theta}_1\dot{\theta}_2\right) + 2\sin{\Delta}\left(\dfrac{r_1s\dot{r}_2\dot{\theta}_1}{r_2}-\dot{r}_1s\dot{\theta}_2\right).
\end{align*}

### Centre of mass
\begin{align*}
x_{2r} &= r_1\cos{\theta_1} + \dfrac{r_2\cos{\theta_2}}{2} & \dot{x}_{2r} = \dot{r_1}\cos{\theta_1} - r_1\dot{\theta}_1 \sin{\theta_1} + \dfrac{\dot{r_2}\cos{\theta_2}}{2} - \dfrac{r_2\dot{\theta}_2 \sin{\theta_2}}{2} \\
y_{2r} &= r_1\sin{\theta_1} + \dfrac{r_2\sin{\theta_2}}{2} & \dot{y}_{2r} = \dot{r_1}\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dfrac{\dot{r_2}\sin{\theta_2}}{2} + \dfrac{r_2\dot{\theta}_2 \cos{\theta_2}}{2} \\
\end{align*}

\begin{align*}
|\vec{v}_{2r}|^2 &= \dot{r}_1^2\cos^2{\theta_1} - 2\dot{r}_1r_1\dot{\theta}_1 \sin{\theta_1}\cos{\theta_1} + r_1^2\dot{\theta}_1^2\sin^2{\theta_1} + \dfrac{\dot{r}_2^2\cos^2{\theta_2}}{4} + \dfrac{r_2^2\dot{\theta}_2^2\sin^2{\theta_2}}{4} - \dfrac{r_2\dot{r}_2\dot{\theta}_2 \cos{\theta_2}\sin{\theta_2}}{2} + \\
& \dot{r}_1\dot{r}_2\cos{\theta_1}\cos{\theta_2} - \dot{r}_1r_2\dot{\theta}_2 \cos{\theta_1}\sin{\theta_2} + r_1\dot{r}_2 \dot{\theta}_1 \sin{\theta_1}\cos{\theta_2} + r_1r_2 \dot{\theta}_1 \dot{\theta}_2 \sin{\theta_1}\sin{\theta_2} + \dot{r}_1^2\sin^2{\theta_1} + \\
& r_1^2 \dot{\theta}_1^2\cos^2{\theta_1} + \dfrac{\dot{r}_2^2\sin^2{\theta_2}}{4} + \dfrac{r_2^2\dot{\theta}_2^2\cos^2{\theta_2}}{4} + 2\dot{r}_1r_1 \dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + \dot{r}_1\dot{r}_2\sin{\theta_1}\sin{\theta_2} + \dot{r}_1 r_2\dot{\theta}_2 \sin{\theta_1}\cos{\theta_2} + \\
& r_1\dot{r}_2 \dot{\theta}_1 \cos{\theta_1}\sin{\theta_2} + r_1r_2\dot{\theta_1}\dot{\theta_2}\cos{\theta_1}\cos{\theta_2} + \dfrac{r_2\dot{r}_2\dot{\theta}_2\sin{\theta_2}\cos{\theta_2}}{2}\\
&= \dot{r}_1^2 + r_1^2\dot{\theta}_1^2 + \dfrac{\dot{r}_2^2 + r_2^2\dot{\theta}_2^2}{4} + r_1\dot{r}_2\dot{\theta}_1\sin{\Delta} + \dot{r}_1\dot{r}_2 \cos{\Delta} + r_1r_2 \dot{\theta}_1\dot{\theta}_2 \cos{\Delta} - \dot{r}_1r_1 \dot{\theta}_1\sin{\Delta} \\
&= \dot{r}_1^2 + r_1^2\dot{\theta}_1^2 + \dfrac{\dot{r}_2^2 + r_2^2\dot{\theta}_2^2}{4} + \sin{\Delta}(r_1\dot{r}_2\dot{\theta}_1 - \dot{r}_1r_1 \dot{\theta}_1) + \cos{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2 \dot{\theta}_1\dot{\theta}_2) \\
&= |\vec{v}_{1b}|^2 + \dfrac{|\vec{v}_{2b}^I|^2}{4} + \dfrac{|\Delta \vec{v}_{2,1}|^2}{2}
\end{align*}

# Dissipative forces
We will assume that the dissipative forces are proportion to the velocity and velocity squared of the pendulum bobs. Meaning they will have the form

\begin{align*}
\vec{F}_{D,j} &= -(b_j+c_j|\vec{v}_j|)\vec{v}_j.
\end{align*}

Where $j$ is the pendulum bob of interest, $b_j$ and $c_j$ are constants. 

# Kinetic energy
The kinetic energy of the system is given by

\begin{align*}
	T &= \dfrac{m_{1b}}{2}|\vec{v}_{1b}|^2 + \int_0^{r_1}\dfrac{m_{1r}}{2r_1} |\vec{v}_{1r}|^2 ds + \dfrac{m_{2b}}{2}|\vec{v}_{2b}|^2 + \int_0^{r_2} \dfrac{m_{2r}}{2r_2} |\vec{v}_{2r}|^2 ds
\end{align*}

## Rod 1
\begin{align*}
\int_0^{r_1}\dfrac{m_{1r}}{2r_1} |\vec{v}_{1r}|^2 ds &= \int_0^{r_1} \dfrac{m_{1r}}{2r_1} s^2 \left[\dfrac{\dot{r}_1^2}{r_1^2}+\dot{\theta}_1^2\right] ds \\
&= \dfrac{m_{1r}}{2r_1} \left[\dfrac{\dot{r}_1^2}{r_1^2}+\dot{\theta}_1^2\right] \left[\dfrac{s^3}{3}\right]_0^{r_1} \\
&= \dfrac{m_{1r}}{2r_1} \left[\dfrac{\dot{r}_1^2}{r_1^2}+\dot{\theta}_1^2\right] \dfrac{r_1^3}{3} \\
&= \dfrac{m_{1r}}{6}\left[\dot{r}_1^2+r_1^2\dot{\theta}_1^2\right] \\
&= \dfrac{m_{1r}}{6} |\vec{v}_{1b}|^2.
\end{align*}

## Rod 2
\begin{align*}
\int_0^{r_2} \dfrac{m_{2r}}{2r_2} |\vec{v}_{2r}|^2 ds &= \int_0^{r_2} \dfrac{m_{2r}}{2r_2} \left[\dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dfrac{\dot{r}_2^2 s^2}{r_2^2} + s^2\dot{\theta}_2^2 + 2\cos{\Delta}\left(\dfrac{\dot{r}_1\dot{r}_2 s}{r_2} + r_1s\dot{\theta}_1\dot{\theta}_2\right) + 2\sin{\Delta}\left(\dfrac{r_1s\dot{r}_2\dot{\theta}_1}{r_2}-\dot{r}_1s\dot{\theta}_2\right)\right] ds \\
&= \dfrac{m_{2r}}{2} |\vec{v}_{1b}|^2 + \dfrac{m_{2r}}{2r_2}\int_0^{r_2} s^2\left(\dfrac{\dot{r}_2^2}{r_2^2} + \dot{\theta}_2^2\right) + 2s\left[\cos{\Delta}\left(\dfrac{\dot{r}_1\dot{r}_2 }{r_2} + r_1\dot{\theta}_1\dot{\theta}_2\right) + \sin{\Delta}\left(\dfrac{r_1\dot{r}_2\dot{\theta}_1}{r_2}-\dot{r}_1\dot{\theta}_2\right)\right]ds \\
&= \dfrac{m_{2r}}{2} |\vec{v}_{1b}|^2 + \dfrac{m_{2r}}{2r_2}\left[\dfrac{r_2^3}{3}\left(\dfrac{\dot{r}_2^2}{r_2^2} + \dot{\theta}_2^2\right) + r_2^2\left[\cos{\Delta}\left(\dfrac{\dot{r}_1\dot{r}_2 }{r_2} + r_1\dot{\theta}_1\dot{\theta}_2\right) + \sin{\Delta}\left(\dfrac{r_1\dot{r}_2\dot{\theta}_1}{r_2}-\dot{r}_1\dot{\theta}_2\right)\right]\right] \\
&= \dfrac{m_{2r}}{2} |\vec{v}_{1b}|^2 + \dfrac{m_{2r}}{6} |\vec{v}_{2b}^{I}|^2 + \dfrac{m_{2r}}{4}|\Delta \vec{v}_{2,1}^I|^2
\end{align*}

## Final result
\begin{align*}
T &= \dfrac{m_{1b}}{2}|\vec{v}_{1b}|^2 + \dfrac{m_{1r}}{6} |\vec{v}_{1b}|^2 + \dfrac{m_{2b}}{2}(|\vec{v}_{1b}|^2+|\vec{v}_{2b}^I|^2+|\Delta \vec{v}_{2,1}^I|^2) + \dfrac{m_{2r}}{2} |\vec{v}_{1b}|^2 + \dfrac{m_{2r}}{6} |\vec{v}_{2b}^{I}|^2 + \dfrac{m_{2r}}{4}|\Delta \vec{v}_{2,1}^I|^2 \\
&= \dfrac{m_{1b}+\dfrac{m_{1r}}{3}+m_{2b}+m_{2r}}{2} |\vec{v}_{1b}|^2 + \dfrac{m_{2b} + \dfrac{m_{2r}}{3}}{2}|\vec{v}_{2b}^I|^2 + \dfrac{m_{2b}+\dfrac{m_{2r}}{2}}{2}|\Delta \vec{v}_{2,1}^I|^2.
\end{align*}

Let $M_j = \displaystyle \sum_{i=jb}^{2r} m_i \left(1-\dfrac{2\delta_{i,jr}}{3}\right)$ and $\mu_j = \displaystyle \sum_{i=jb}^{2r} m_i \left(1-\dfrac{\delta_{i,jr}}{2}\right)$, where $\delta_{i,j}$ is the Kronecker delta symbol. For instance, $M_1 = m_{1b} + \dfrac{m_{1r}}{3} + m_{2b} + m_{2r}$ and $\mu_2 = m_{2b} + \dfrac{m_{2r}}{2}$. Then

\begin{align*}
T &= \dfrac{M_1}{2} |\vec{v}_{1b}|^2 + \dfrac{M_2}{2}|\vec{v}_{2b}^I|^2 + \dfrac{\mu_2}{2}|\Delta \vec{v}_{2,1}^I|^2.
\end{align*}

# Potential energy
The potential energy of the system is given by

\begin{align*}
	V &= m_{1b} gy_{1b} + m_{1r}gy_{1r} + m_{2b} gy_{2b} + m_{2r}gy_{2r} + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}\\
	&= m_{1b} gr_1\sin{\theta_1} + \dfrac{m_{1r} gr_1\sin{\theta_1}}{2} + m_{2r}g\left(r_1\sin{\theta_1} + \dfrac{r_2\sin{\theta_2}}{2}\right) + m_{2b}g(r_1\sin{\theta_1} + r_2\sin{\theta_2}) + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}\\
	&= \left(m_{1b} + \dfrac{m_{1r}}{2} + m_{2b} + m_{2r}\right)gr_1\sin{\theta_1} + \left(m_{2b}+\dfrac{m_{2r}}{2}\right)gr_2\sin{\theta_2} + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2} \\
	&= \mu_1 gr_1\sin{\theta_1} + \mu_2 gr_2\sin{\theta_2} + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}
\end{align*}

# Lagrangian
Hence the Lagrangian of the system is

\begin{align*}
	\mathcal{L} &= T - V \\
	&= \dfrac{M_1}{2} |\vec{v}_{1b}|^2 + \dfrac{M_2}{2}|\vec{v}_{2b}^I|^2 + \dfrac{\mu_2}{2}|\Delta \vec{v}_{2,1}^I|^2 - \left(\mu_1 gr_1\sin{\theta_1} + \mu_2 gr_2\sin{\theta_2} + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}\right) \\
	&= \dfrac{M_1}{2} |\vec{v}_{1b}|^2 + \dfrac{M_2}{2}|\vec{v}_{2b}^I|^2 + \dfrac{\mu_2}{2}(|\Delta \vec{v}_{2,1}^I|^2 - 2gr_2 \sin{\theta_2}) - \mu_1 gr_1 \sin{\theta_1} - \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}.
\end{align*}

We will not expand this Lagrangian, as doing so just adds to its complexity. Instead, we will calculate the derivatives of each of its components. 

# Derivative of components of the Lagrangian
## Square of the velocity of the first pendulum's bob
The relevant partial and standard derivatives are:
\begin{align*}
	\dfrac{\partial |\vec{v}_{1b}|^2}{\partial r_1} &= 2r_1\dot{\theta}_1^2 & \dfrac{\partial |\vec{v}_{1b}|^2}{\partial r_2} &= 0 & \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \theta_1} &= 0 & \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \theta_2} &= 0\\
	\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{r}_1} &= 2\dot{r}_1 & \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{r}_2} &= 0 & \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{\theta}_1} &= 2r_1^2\dot{\theta}_1 & \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{\theta}_2} &= 0\\
	\dfrac{d}{dt} \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{r}_1} &= 2\ddot{r}_1 & \dfrac{d}{dt}\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{r}_2} &= 0 & \dfrac{d}{dt}\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{\theta}_1} &= 2r_1^2 \ddot{\theta}_1 + 4r_1\dot{r}_1\dot{\theta}_1 & \dfrac{d}{dt}\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{\theta}_2} &= 0
\end{align*}

Hence the negative functional derivatives are

\begin{align*}
	\dfrac{\delta' |\vec{v}_{1b}|^2}{\delta' r_1} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{r}_1} - \dfrac{\partial |\vec{v}_{1b}|^2}{\partial r_1} & \dfrac{\delta' |\vec{v}_{1b}|^2}{\delta' r_2} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{r}_2} - \dfrac{\partial |\vec{v}_{1b}|^2}{\partial r_2} & \dfrac{\delta' |\vec{v}_{1b}|^2}{\delta' \theta_1} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{\theta}_1} - \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \theta_1} & \dfrac{\delta' |\vec{v}_{1b}|^2}{\delta' \theta_2} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_{1b}|^2}{\partial \dot{\theta}_2} - \dfrac{\partial |\vec{v}_{1b}|^2}{\partial \theta_2}\\
	&= 2\ddot{r}_1 - 2r_1\dot{\theta}_1^2 & &= 0 & &=2r_1^2\ddot{\theta}_1 + 4r_1\dot{r}_1\dot{\theta}_1 & &= 0.	
\end{align*}

## Second bob's independent velocity squared
Hence the partial and standard derivatives of the difference in the square of each bob's velocity is
\begin{align*}
	\dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial r_1} &= 0 & \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial r_2} &= 2r_2\dot{\theta}_2^2 \\
	\dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \theta_1} &= 0 & \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \theta_2} &= 0 \\
	\dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{r}_1} &= 0 & \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{r}_2} &= 2\dot{r}_2 \\
	\dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{\theta}_1} &= 0 & \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{\theta}_2} &=2r_2^2\dot{\theta}_2
\end{align*}

\begin{align*}
\dfrac{d}{dt} \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{r}_1} &= 0 & \dfrac{d}{dt} \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{r}_2} &= 2\ddot{r}_2 \\
\dfrac{d}{dt} \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{\theta}_1} &= 0 & \dfrac{d}{dt} \dfrac{\partial |\vec{v}_{2b}^I|^2}{\partial \dot{\theta}_2} &= 2r_2^2\ddot{\theta}_2 + 4r_2\dot{r}_2\dot{\theta}_2.
\end{align*}

Hence

\begin{align*}
\delta_{r_1} |\vec{v}_{2b}^I|^2 &= 0 & \delta_{r_2} |\vec{v}_{2b}^I|^2 &= 2(\ddot{r}_2 - r_2\dot{\theta}_2^2) \\
\delta_{\theta_1} |\vec{v}_{2b}^I|^2 &=0 & \delta_{\theta_2} |\vec{v}_{2b}^I|^2 &= 
2r_2^2 \ddot{\theta}_2 + 4r_2\dot{r}_2\dot{\theta}_2.
\end{align*}

## Difference in independent velocities squared
\begin{align*}
\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial r_1} &= 2r_2\dot{\theta}_1\dot{\theta}_2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_1\sin{\Delta} & \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial r_2} &= 2r_1\dot{\theta}_1\dot{\theta}_2\cos{\Delta}-2\dot{r}_1\dot{\theta}_2\sin{\Delta}\\
\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \theta_1} &= -2\sin{\Delta}\cdot -1(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\cos{\Delta}\cdot -1(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2) & \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \theta_2} &= -2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)  \\
&= 2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) - 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)\\
\dfrac{\partial |\Delta \vec{v}_{2,1}|^2}{\partial \dot{r}_1} &= 2\dot{r}_2\cos{\Delta} - 2r_2\dot{\theta}_2\sin{\Delta} & \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{r}_2} &= 2\dot{r}_1\cos{\Delta} + 2r_1\dot{\theta}_1\sin{\Delta} \\
\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{\theta}_1} &= 2r_1r_2\dot{\theta}_2\cos{\Delta} + 2r_1\dot{r}_2\sin{\Delta} & \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{\theta}_2} &=2r_1r_2\dot{\theta}_1\cos{\Delta} - 2\dot{r}_1r_2\sin{\Delta}
\end{align*}

Let us define $\dot{\Delta}_1 = 2\dot{\theta}_1 - \dot{\theta}_2$ and $\dot{\Delta}_2 = 2\dot{\theta}_2 - \dot{\theta}_1$.

\begin{align*}
	\dfrac{d}{dt} \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{r}_1} &= 2\ddot{r}_2\cos{\Delta} - 2\dot{r}_2\dot{\Delta}\sin{\Delta} - 2\dot{r}_2\dot{\theta}_2\sin{\Delta} - 2r_2\ddot{\theta}_2\sin{\Delta} - 2r_2\dot{\theta}_2\dot{\Delta}\cos{\Delta} \\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2\dot{\Delta} + \dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)\\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2(2\dot{\theta}_2-\dot{\theta}_1)+r_2\ddot{\theta}_2)\\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2\dot{\Delta}_2+r_2\ddot{\theta}_2)\\
	\dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{r}_2} &= 2\ddot{r}_1 \cos{\Delta} -2\dot{r}_1\dot{\Delta}\sin{\Delta} + 2\dot{r}_1\dot{\theta}_1\sin{\Delta} + 2r_1\ddot{\theta}_1\sin{\Delta} + 2r_1\dot{\theta}_1\dot{\Delta}\cos{\Delta} \\
	&= 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1\dot{\theta}_1 - \dot{r}_1\dot{\Delta})\\
	&= 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1(2\dot{\theta}_1 - \ddot{\theta}_2))\\
	&= 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1\dot{\Delta}_1)\\
	\dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{\theta}_1} &= 2\dot{r}_1r_2\dot{\theta}_2\cos{\Delta} + 2r_1\dot{r}_2\dot{\theta}_2\cos{\Delta} + 2r_1r_2\ddot{\theta}_2\cos{\Delta} - 2r_1r_2\dot{\theta}_2\dot{\Delta}\sin{\Delta} + 2\dot{r}_1\dot{r}_2\sin{\Delta} + 2r_1\ddot{r}_2\sin{\Delta} + 2r_1\dot{r}_2\dot{\Delta}\cos{\Delta} \\
	&= 2\cos{\Delta}(\dot{r}_1r_2\dot{\theta}_2 + r_1\dot{r}_2 (\dot{\theta}_2+\dot{\Delta})+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2\dot{\Delta})\\
	&= 2\cos{\Delta}(\dot{r}_1r_2\dot{\theta}_2 + r_1\dot{r}_2 \dot{\Delta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2\dot{\Delta})\\
	\dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{\theta}_2} &= 2\dot{r}_1r_2\dot{\theta}_1\cos{\Delta} + 2r_1\dot{r}_2\dot{\theta}_1\cos{\Delta} + 2r_1r_2\ddot{\theta}_1\cos{\Delta} - 2r_1r_2\dot{\theta}_1\dot{\Delta}\sin{\Delta} - 2\ddot{r}_1r_2\sin{\Delta} - 2\dot{r}_1\dot{r}_2\sin{\Delta} - 2\dot{r}_1r_2\dot{\Delta}\cos{\Delta} \\
	&=2\cos{\Delta}(\dot{r}_1r_2(\dot{\theta}_1-\dot{\Delta})+r_1\dot{r}_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)-2\sin{\Delta}(r_1r_2\dot{\theta}_1\dot{\Delta} + \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2) \\
	&=2\cos{\Delta}(\dot{r}_1r_2\dot{\Delta}_1+r_1\dot{r}_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)-2\sin{\Delta}(r_1r_2\dot{\theta}_1\dot{\Delta} + \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2).
\end{align*}

Hence the negative functional derivative for $r_1$ is
\begin{align*}
	\delta'_{r_1} |\Delta \vec{v}_{2,1}^I|^2 &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{r}_1} - \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial r_1} \\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2\dot{\Delta}_2+r_2\ddot{\theta}_2) - \left(2r_2\dot{\theta}_1\dot{\theta}_2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_1\sin{\Delta}\right) \\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2(\dot{\Delta}+\dot{\theta}_1)) - 2\sin{\Delta} (\dot{r}_2(\dot{\Delta}_2+\dot{\theta}_1)+r_2\ddot{\theta}_2).
\end{align*}

Where $\dot{\Delta} + \dot{\theta}_1 = \dot{\theta}_2 - \dot{\theta}_1 + \dot{\theta}_1 = \dot{\theta}_2$ and $\dot{\Delta}_2 + \dot{\theta}_1 = 2\dot{\theta}_2 - \dot{\theta}_1 + \dot{\theta}_1 = 2\dot{\theta}_2$. (Confirmed with SymPy)

\begin{align*}
	\delta'_{r_1} |\Delta \vec{v}_{2,1}^I|^2 &= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - 2\sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2).
\end{align*}

As for $r_2$ (checked)

\begin{align*}
	\delta'_{r_2} |\Delta \vec{v}_{2,1}^I|^2 &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{r}_2} - \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial r_2} \\
	&= 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1\dot{\Delta}_1) - \left[2r_1\dot{\theta}_1\dot{\theta}_2\cos{\Delta}-2\dot{r}_1\dot{\theta}_2\sin{\Delta}\right]\\
	&= 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1( \dot{\Delta}-\dot{\theta}_2))+2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1(\dot{\Delta}_1+\dot{\theta}_2)).
\end{align*}

Hence $\dot{\Delta}-\dot{\theta}_2 = \dot{\theta}_2-\dot{\theta}_1-\dot{\theta}_2 = -\dot{\theta}_1$ and $\dot{\Delta}_1 + \dot{\theta}_2 = 2\dot{\theta}_1 - \dot{\theta}_2 + \dot{\theta}_2 = 2\dot{\theta}_1$. 

\begin{align*}
	\delta'_{r_2} |\Delta \vec{v}_{2,1}^I|^2 &= 2\cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+2\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1).
\end{align*}

As for $\theta_1$ (confirmed by SymPy)

\begin{align*}
	\delta'_{\theta_1} |\Delta \vec{v}_{2,1}^I|^2 &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{\theta}_1} - \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \theta_1} \\
	&= 2\cos{\Delta}(\dot{r}_1r_2\dot{\theta}_2 + r_1\dot{r}_2 \dot{\Delta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2\dot{\Delta}) - \left[2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) - 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)\right] \\
	&= 2\cos{\Delta}(\dot{r}_1r_2(\dot{\theta}_2 - \dot{\theta}_2)+ r_1\dot{r}_2 (\dot{\Delta}_2+\dot{\theta}_1)+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2-\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2(\dot{\Delta}+\dot{\theta}_1)) \\
	&= 2\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2).
\end{align*}

As for $\theta_2$ (correct, checked with SymPy)

\begin{align*}
	\delta'_{\theta_2} |\Delta \vec{v}_{2,1}^I|^2 &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \dot{\theta}_2} - \dfrac{\partial |\Delta \vec{v}_{2,1}^I|^2}{\partial \theta_2} \\
	&= 2\cos{\Delta}(\dot{r}_1r_2\dot{\Delta}_1+r_1\dot{r}_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)-2\sin{\Delta}(r_1r_2\dot{\theta}_1\dot{\Delta} + \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2) - \left[-2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)\right] \\
	&= 2\cos{\Delta}(\dot{r}_1r_2(\dot{\Delta}_1+\dot{\theta}_2)+r_1\dot{r}_2(\dot{\theta}_1-\dot{\theta}_1)+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1(\dot{\theta}_2-\dot{\Delta}) - \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2-\dot{r}_1\dot{r}_2) \\
	&= 2\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1(\dot{\theta}_2-(\dot{\theta}_2-\dot{\theta}_1))-\ddot{r}_1r_2) \\
	&= 2\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2).
\end{align*}

# Euler-Lagrange equations with dissipation
## $r_1$
It is important to note that $\dfrac{\delta' f(q_i)}{\delta' q_i} = -\dfrac{\partial f}{\partial q_i}$ and of course if a term does not depend on $q_i$ or $\dot{q}_i$ its functional derivative with respect to $q_i$ is zero. Hence

\begin{align*}
	\delta'_{r_1} \mathcal{L} &= \dfrac{M_1}{2} \delta'_{r_1} |\vec{v}_{1b}|^2 + \dfrac{M_2}{2} \delta'_{r_1}|\vec{v}_{2b}^I|^2 + \dfrac{\mu_2}{2}\delta'_{r_1} |\Delta \vec{v}_{2,1}^I|^2 + \mu_1 g\sin{\theta_1} + k_1(r_1-l_1)
\end{align*}
We have deliberately ignored the $m_2gr_2\sin{\theta_2}$ and $-\dfrac{k_2(r_2-l_2)^2}{2}$ as they are independent of $r_1$.
\begin{align*}
	\delta'_{r_1} \mathcal{L} &= \dfrac{M_1}{2}(2\ddot{r}_1-2r_1\dot{\theta}_1^2) + \dfrac{M_2}{2} \cdot 0 + \dfrac{\mu_2}{2} \left(2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2(\dot{\Delta}+\dot{\theta}_1)) - 2\sin{\Delta} (\dot{r}_2(\dot{\Delta}_2+\dot{\theta}_1)+r_2\ddot{\theta}_2)\right) + \mu_1g\sin{\theta_1} + k_1(r_1-l_1) \\
	&= M_1(\ddot{r}_1-r_1\dot{\theta}_1^2) + \mu_2 \left(\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2(\dot{\Delta}+\dot{\theta}_1)) - \sin{\Delta} (\dot{r}_2(\dot{\Delta}_2+\dot{\theta}_1)+r_2\ddot{\theta}_2)\right) + \mu_1g\sin{\theta_1} + k_1(r_1-l_1) \\
	&= M_1(\ddot{r}_1-r_1\dot{\theta}_1^2) + \mu_2 \left(\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - \sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)\right) + \mu_1g\sin{\theta_1} + k_1(r_1-l_1).
\end{align*}
	
The generalized dissipation force canonical to $r_1$ is hence
\begin{align*}
	Q_{r_1} &= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\vec{v}_{1b}\cdot \dfrac{\partial \vec{r}_{1b}}{\partial r_{1}} -(b_{1r}+c_{1r}|\vec{v}_{1r}|)\vec{v}_{1r}\cdot \dfrac{\partial \vec{r}_{1r}}{\partial r_{1}} - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\vec{v}_{2b}\cdot \dfrac{\partial \vec{r}_{2b}}{\partial r_1} - (b_{2r}+c_{2r}|\vec{v}_{2r}|)\vec{v}_{2r}\cdot \dfrac{\partial \vec{r}_{2r}}{\partial r_1} \\
	&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}
	\end{bmatrix} \cdot \begin{bmatrix}
		\cos{\theta_1} \\
		\sin{\theta_1}
	\end{bmatrix} -(b_{1r}+c_{1r}|\vec{v}_{1r}|)\dfrac{1}{2}\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}
	\end{bmatrix} \cdot \dfrac{1}{2} \begin{bmatrix}
		\cos{\theta_1} \\
		\sin{\theta_1}
	\end{bmatrix} - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot \begin{bmatrix}
		\cos{\theta_1} \\
		\sin{\theta_1}
	\end{bmatrix}- (b_{2r}+c_{2r}|\vec{v}_{2r}|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dfrac{\dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2}}{2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dfrac{\dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}}{2}
	\end{bmatrix} \cdot \begin{bmatrix}
		\cos{\theta_1} \\
		\sin{\theta_1}
	\end{bmatrix} \\
	&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\left[\dot{r}_1\cos^2{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + \dot{r}_1\sin^2{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_1}\right] \\
	&- \dfrac{1}{4}(b_{1r}+c_{1r}|\vec{v}_{1r}|)\left[\dot{r}_1\cos^2{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + \dot{r}_1\sin^2{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_1}\right] \\
	&- (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1\cos^2{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + \dot{r}_2\cos{\theta_1}\cos{\theta_2}-r_2\dot{\theta}_2\cos{\theta_1}\sin{\theta_2} \right.\\
	&\left.+ \dot{r}_1\sin^2{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_1} + \dot{r}_2\sin{\theta_1}\sin{\theta_2} + r_2\dot{\theta}_2\sin{\theta_1}\cos{\theta_2} \right] \\
	&- (b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1\cos^2{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + \dfrac{1}{2}\left[\dot{r}_2\cos{\theta_1}\cos{\theta_2}-r_2\dot{\theta}_2\cos{\theta_1}\sin{\theta_2}\right]\right.\\
	&\left.+ \dot{r}_1\sin^2{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_1} + \dfrac{\dot{r}_2\sin{\theta_1}\sin{\theta_2} + r_2\dot{\theta}_2\sin{\theta_1}\cos{\theta_2}}{2} \right] \\
	&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\left[\dot{r}_1(\cos^2{\theta_1}+\sin^2{\theta_1}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_1} + \cos{\theta_1}\sin{\theta_1})\right] \\
	&- \dfrac{1}{4}(b_{1r}+c_{1r}|\vec{v}_{1r}|)\left[\dot{r}_1(\cos^2{\theta_1}+\sin^2{\theta_1}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_1} + \cos{\theta_1}\sin{\theta_1})\right] \\
	&- (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1(\cos^2{\theta_1}+\sin^2{\theta_1}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_1} + \cos{\theta_1}\sin{\theta_1}) \right.\\
	&\left.+ \dot{r}_2(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2})+r_2\dot{\theta}_2(-\cos{\theta_1}\sin{\theta_2} + \sin{\theta_1}\cos{\theta_2}) \right] \\
	&- (b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1(\cos^2{\theta_1}+\sin^2{\theta_1}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_1} + \cos{\theta_1}\sin{\theta_1})\right.\\
	&\left.+ \dfrac{1}{2}\left[\dot{r}_2(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2})+r_2\dot{\theta}_2(-\cos{\theta_1}\sin{\theta_2}+\sin{\theta_1}\cos{\theta_2})\right]\right]\\
	&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\dot{r}_1 - \dfrac{\dot{r}_1}{4}(b_{1r}+c_{1r}|\vec{v}_{1r}|) - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1 + \dot{r}_2\cos{\Delta}-r_2\dot{\theta}_2\sin{\Delta} \right] \\
	&- (b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1+ \dfrac{\dot{r}_2\cos{\Delta}-r_2\dot{\theta}_2\sin{\Delta}}{2}\right].
\end{align*}

Hence the Euler-Lagrange equation for $r_1$ with dissipative forces is

\begin{align*}
	M_1(\ddot{r}_1-r_1\dot{\theta}_1^2) + \mu_2 \left(\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - \sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)\right) + \mu_1g\sin{\theta_1} + k_1(r_1-l_1) &= Q_{r_1}.
\end{align*}

So, this line of our matrix equation for $\mathbf{\ddot{q}}$ will be

\begin{align*}
\begin{bmatrix}
M_1 & \mu_2 \cos{\Delta} & 0 & -\mu_2r_2 \sin{\Delta}
\end{bmatrix} \begin{bmatrix}
\ddot{r}_1 \\
\ddot{r}_2 \\
\ddot{\theta}_1 \\
\ddot{\theta}_2
\end{bmatrix} &= \begin{bmatrix}
M_1r_1^2 \dot{\theta}_1^2 + \mu_2 r_2 \dot{\theta}_2^2\cos{\Delta} + 2\mu_2 \dot{r}_2 \dot{\theta}_2\sin{\Delta} - \mu_1g\sin{\theta_1} - k_1(r_1-l_1) + Q_{r_1}
\end{bmatrix}.
\end{align*}

## $r_2$
As for $r_2$

\begin{align*}
	\delta'_{r_2} \mathcal{L} &= \dfrac{M_1}{2} \delta'_{r_2} |\vec{v}_{1b}|^2 + \dfrac{M_2}{2} \delta'_{r_2} |\vec{v}_{2b}^I|^2 + \dfrac{\mu_2}{2}(\delta'_{r_2} |\Delta \vec{v}_{2,1}^I|^2+2g\sin{\theta_2}) + k_2(r_2-l_2) \\
	&=  \dfrac{M_2}{2} (2(\ddot{r}_2 - r_2\dot{\theta}_2^2)) + \dfrac{\mu_2}{2}(2\cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+2\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1)+2g\sin{\theta_2}) + k_2(r_2-l_2) \\
	&= M_2 (\ddot{r}_2 - r_2\dot{\theta}_2^2) + \mu_2 (\cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1)+g\sin{\theta_2}) + k_2(r_2-l_2) \\
	&= M_2 (\ddot{r}_2 - r_2\dot{\theta}_2^2) + \mu_2 (\cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1)+g\sin{\theta_2}) + k_2(r_2-l_2).
\end{align*}

Generalized dissipative forces

\begin{align*}
Q_{r_2} &= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\vec{v}_{1b} \cdot \dfrac{\partial \vec{r}_{1b}}{\partial r_2}-(b_{1r}+c_{1r}|\vec{v}_{1r}|)\vec{v}_{1r} \cdot \dfrac{\partial \vec{r}_{1r}}{\partial r_2} - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\vec{v}_{2b} \cdot \dfrac{\partial \vec{r}_{2b}}{\partial r_{2}} - (b_{2r}+c_{2r}|\vec{v}_{2r}|)\vec{v}_{2r} \cdot \dfrac{\partial \vec{r}_{2r}}{\partial r_{2}} \\
	&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}
	\end{bmatrix} \cdot \vec{0} -(b_{1r}+c_{1r}|\vec{v}_{1r}|)\dfrac{1}{2}\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}
	\end{bmatrix} \cdot \vec{0} - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot \begin{bmatrix}
		\cos{\theta_2} \\
		\sin{\theta_2}
	\end{bmatrix} - (b_{2r}+c_{2r}|\vec{v}_{2r}|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dfrac{\dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2}}{2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dfrac{\dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}}{2}
	\end{bmatrix} \cdot \dfrac{1}{2}\begin{bmatrix}
		\cos{\theta_2} \\
		\sin{\theta_2}
	\end{bmatrix} \\
	&=  - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1\cos{\theta_1}\cos{\theta_2} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_2} + \dot{r}_2\cos^2{\theta_2}-r_2\dot{\theta}_2\sin{\theta_2}\cos{\theta_2} + \dot{r}_1\sin{\theta_1}\sin{\theta_2} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_2} + \dot{r}_2\sin^2{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}\sin{\theta_2}\right] - \dfrac{1}{2}(b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1\cos{\theta_1}\cos{\theta_2} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_2} + \dfrac{\dot{r}_2\cos^2{\theta_2}-r_2\dot{\theta}_2\sin{\theta_2}\cos{\theta_2}}{2} + \dot{r}_1\sin{\theta_1}\sin{\theta_2} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_2} + \dfrac{\dot{r}_2\sin^2{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}\sin{\theta_2}}{2}\right] \\
	&= - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_2} + \cos{\theta_1}\sin{\theta_2}) + \dot{r}_2(\cos^2{\theta_2} + \sin^2{\theta_2})+r_2\dot{\theta}_2(-\sin{\theta_2}\cos{\theta_2} +  \cos{\theta_2}\sin{\theta_2})\right] - \dfrac{1}{2}(b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_2} + \cos{\theta_1}\sin{\theta_2}) + \dfrac{\dot{r}_2(\cos^2{\theta_2} + \sin^2{\theta_2})+r_2\dot{\theta}_2(-\sin{\theta_2}\cos{\theta_2} +  \cos{\theta_2}\sin{\theta_2})}{2}\right] \\
	&= - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1\cos{(\theta_2-\theta_1)} + r_1\dot{\theta}_1\sin{(\theta_2-\theta_1)} + \dot{r}_2\right] - \dfrac{1}{2}(b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1\cos{\Delta} + r_1\dot{\theta}_1\sin{\Delta} + \dfrac{\dot{r}_2}{2}\right]\\
	&= - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1\cos{\Delta} + r_1\dot{\theta}_1\sin{\Delta} + \dot{r}_2\right]- \dfrac{1}{2}(b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1\cos{\Delta} + r_1\dot{\theta}_1\sin{\Delta} + \dfrac{\dot{r}_2}{2}\right].
\end{align*}

Hence the Euler-Lagrange equation for $r_2$ with dissipative forces is

\begin{align*}
	M_2 (\ddot{r}_2 - r_2\dot{\theta}_2^2) + \mu_2 (\cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1)+g\sin{\theta_2}) + k_2(r_2-l_2) &= Q_{r_2}.
\end{align*}

Or, in matrix form

\begin{align*}
\begin{bmatrix}
\mu_2 \cos{\Delta} & M_2 & \mu_2 r_1\sin{\Delta} & 0
\end{bmatrix} \begin{bmatrix}
\ddot{r}_1 \\
\ddot{r}_2 \\
\ddot{\theta}_1 \\
\ddot{\theta}_2
\end{bmatrix} &= \begin{bmatrix}
M_2r_2\dot{\theta}_2^2 + \mu_2 (r_1\dot{\theta}_1^2\cos{\Delta}-2\dot{r}_1\dot{\theta}_1\sin{\Delta} - g\sin{\theta_2}) - k_2(r_2-l_2) + Q_{r_2}
\end{bmatrix}.
\end{align*}

## $\theta_1$
As for $\theta_1$

\begin{align*}
\delta'_{\theta_1} \mathcal{L} &= \dfrac{M_1}{2} \delta'_{\theta_1} |\vec{v}_{1b}|^2 + \dfrac{M_2}{2} \delta'_{\theta_1} |\vec{v}_{2b}^I|^2 + \dfrac{\mu_2}{2}(\delta'_{\theta_1}|\Delta \vec{v}_{2,1}^I|^2) + \mu_1 gr_1 \cos{\theta_1} \\
&= \dfrac{M_1}{2} (2r_1^2 \ddot{\theta}_1 + 4r_1\dot{r}_1\dot{\theta}_1) + \dfrac{\mu_2}{2}(2\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2)) + \mu_1 gr_1 \cos{\theta_1} \\
&= M_1 (r_1^2 \ddot{\theta}_1 + 2r_1\dot{r}_1\dot{\theta}_1) + \mu_2(\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2)) + \mu_1 gr_1 \cos{\theta_1}.
\end{align*}

Generalized dissipative force

\begin{align*}
	Q_{\theta_1} &= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\vec{v}_{1b} \cdot \dfrac{\partial \vec{r}_{1b}}{\partial \theta_1} -(b_{1r}+c_{1r}|\vec{v}_{1r}|)\vec{v}_{1r} \cdot \dfrac{\partial \vec{r}_{1r}}{\partial \theta_1} - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\vec{v}_{2b} \cdot \dfrac{\partial \vec{r}_{2b}}{\partial \theta_1}- (b_{2r}+c_{2r}|\vec{v}_{2r}|)\vec{v}_{2r} \cdot \dfrac{\partial \vec{r}_{2r}}{\partial \theta_1} \\
	&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta}_1
	\end{bmatrix} \cdot r_1\begin{bmatrix}
		-\sin{\theta_1} \\
		\cos{\theta_1}
	\end{bmatrix} -\dfrac{1}{4}(b_{1r}+c_{1r}|\vec{v}_{1r}|)\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta}_1
	\end{bmatrix} \cdot r_1\begin{bmatrix}
		-\sin{\theta_1} \\
		\cos{\theta_1}
	\end{bmatrix} \\
	&- (b_{2b}+c_{2b}|\vec{v}_{2b}|) \begin{bmatrix}
	\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
	\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot r_1\begin{bmatrix}
	-\sin{\theta_1} \\
	\cos{\theta_1}
	\end{bmatrix} \\
	&- (b_{2r}+c_{2r}|\vec{v}_{2r}|) \begin{bmatrix}
	\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dfrac{\dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2}}{2} \\
	\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dfrac{\dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}}{2}
	\end{bmatrix} \cdot r_1\begin{bmatrix}
	-\sin{\theta_1} \\
	\cos{\theta_1}
	\end{bmatrix}
&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\left[-r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1}+r_1^2\dot{\theta}_1\sin^2{\theta_1} + r_1\dot{r}_1\sin{\theta_1}\cos{\theta_1}+r_1^2\dot{\theta}_1\cos^2{\theta_1}\right] \\
&-\dfrac{1}{4}(b_{1r}+c_{1r}|\vec{v}_{1r}|)\left[-r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1}+r_1^2\dot{\theta}_1\sin^2{\theta_1} + r_1\dot{r}_1\sin{\theta_1}\cos{\theta_1}+r_1^2\dot{\theta}_1\cos^2{\theta_1}\right] \\
&- (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[-r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1} + r_1^2\dot{\theta_1}\sin^2{\theta_1} -r_1\dot{r}_2\sin{\theta_1}\cos{\theta_2}+r_1r_2\dot{\theta}_2\sin{\theta_1}\sin{\theta_2} \right.\\
&\left.+r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1}+r_1^2\dot{\theta}_1\cos^2{\theta_1} + r_1\dot{r}_2\cos{\theta_1}\sin{\theta_2}+r_1r_2\dot{\theta}_2\cos{\theta_1}\cos{\theta_2}\right]\\
&- (b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[-r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1} + r_1^2\dot{\theta_1}\sin^2{\theta_1} -\dfrac{r_1\dot{r}_2\sin{\theta_1}\cos{\theta_2}+r_1r_2\dot{\theta}_2\sin{\theta_1}\sin{\theta_2}}{2} \right.\\
&\left.+r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1}+r_1^2\dot{\theta}_1\cos^2{\theta_1} + \dfrac{r_1\dot{r}_2\cos{\theta_1}\sin{\theta_2}+r_1r_2\dot{\theta}_2\cos{\theta_1}\cos{\theta_2}}{2}\right]
&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\left[r_1\dot{r}_1(-\cos{\theta_1}\sin{\theta_1}+\sin{\theta_1}\cos{\theta_1}) + r_1^2\dot{\theta}_1(\sin^2{\theta_1} +\cos^2{\theta_1})\right] - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[r_1\dot{r}_1(-\cos{\theta_1}\sin{\theta_1} + \cos{\theta_1}\sin{\theta_1}) +r_1^2\dot{\theta_1}(\sin^2{\theta_1}+\cos^2{\theta_1})\right.\\
&\left.+r_1\dot{r}_2(-\sin{\theta_1}\cos{\theta_2}+\cos{\theta_1}\sin{\theta_2})+r_1r_2\dot{\theta}_2(\sin{\theta_1}\sin{\theta_2}+\cos{\theta_1}\cos{\theta_2})\right] \\
&= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)r_1^2\dot{\theta}_1 -\dfrac{1}{4}(b_{1r}+c_{1r}|\vec{v}_{1r}|)r_1^2\dot{\theta}_1 - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[ r_1^2\dot{\theta_1}+r_1\dot{r}_2\sin{\Delta}+r_1r_2\dot{\theta}_2\cos{\Delta}\right] \\
&- (b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[ r_1^2\dot{\theta_1}+\dfrac{r_1\dot{r}_2\sin{\Delta}+r_1r_2\dot{\theta}_2\cos{\Delta}}{2}\right].
\end{align*}

Hence Equation \eqref{ELD} is

\begin{align*}
	M_1 (r_1^2 \ddot{\theta}_1 + 2r_1\dot{r}_1\dot{\theta}_1) + \mu_2(\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2)) + \mu_1 gr_1 \cos{\theta_1} &= Q_{\theta_1}.
\end{align*}

Hence the corresponding line of our matrix equation is

\begin{align*}
\begin{bmatrix}
0 & \mu_2 r_1\sin{\Delta} & M_1r_1^2 & \mu_2 r_1r_2\cos{\Delta}
\end{bmatrix} \begin{bmatrix}
\ddot{r}_1 \\
\ddot{r}_2 \\
\ddot{\theta}_1 \\
\ddot{\theta}_2
\end{bmatrix} &= \begin{bmatrix}
-2M_1r_1\dot{r}_1 \dot{\theta}_1 + \mu_2(r_1r_2\dot{\theta}_2^2 \sin{\Delta} -2r_1\dot{r}_2 \dot{\theta}_2\cos{\Delta}) - \mu_1 gr_1\cos{\theta_1} + Q_{\theta_1}
\end{bmatrix}.
\end{align*}

## $\theta_2$
As for $\theta_2$

\begin{align*}
	\delta'_{\theta_2} \mathcal{L} &= \dfrac{M_1}{2} \delta'_{\theta_2} |\vec{v}_{1b}|^2 + \dfrac{M_2}{2} \delta'_{\theta_2} |\vec{v}_{2b}^I|^2 + \dfrac{\mu_2}{2}(\delta'_{\theta_2} |\Delta \vec{v}_{2,1}^I|^2+2gr_2\cos{\theta_2}) \\
	&= \dfrac{M_2}{2}(2r_2^2 \ddot{\theta}_2 + 4r_2 \dot{r}_2 \dot{\theta}_2) + \dfrac{\mu_2}{2}\left(2\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2)+2gr_2\cos{\theta_2}\right) \\
	&= M_2(r_2^2 \ddot{\theta}_2 + 2r_2 \dot{r}_2 \dot{\theta}_2) + \mu_2 (\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2)+gr_2\cos{\theta_2}).
\end{align*}

\begin{align*}
	Q_{\theta_2} &= -(b_{1b}+c_{1b}|\vec{v}_{1b}|)\vec{v}_{1b} \cdot \dfrac{\partial \vec{r}_{1b}}{\partial \theta_2}-(b_{1r}+c_{1r}|\vec{v}_{1r}|)\vec{v}_{1r} \cdot \dfrac{\partial \vec{r}_{1r}}{\partial \theta_2} - (b_{2b}+c_{2b}|\vec{v}_{2b}|)\vec{v}_{2b}\cdot \dfrac{\partial \vec{r}_{2b}}{\partial \theta_2} -(b_{2r}+c_{2r}|\vec{v}_{2r}|)\vec{v}_{2r}\cdot \dfrac{\partial \vec{r}_{2r}}{\partial \theta_2} \\
	&= -(b_{2b}+c_{2b}|\vec{v}_{2b}|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot r_2\begin{bmatrix}
		-\sin{\theta_2}\\
		\cos{\theta_2}
	\end{bmatrix} -(b_{2r}+c_{2r}|\vec{v}_{2r}|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dfrac{\dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2}}{2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dfrac{\dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}}{2}
	\end{bmatrix} \cdot \dfrac{r_2}{2}\begin{bmatrix}
		-\sin{\theta_2}\\
		\cos{\theta_2}
	\end{bmatrix} \\
	&= -(b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[-\dot{r}_1r_2\cos{\theta}_1\sin{\theta_2} + r_1r_2\dot{\theta}_1\sin{\theta_1}\sin{\theta_2} - r_2\dot{r}_2\cos{\theta_2}\sin{\theta_2} + r_2^2\dot{\theta}_2\sin^2{\theta_2} + \dot{r}_1r_2\sin{\theta_1}\cos{\theta_2} + r_1r_2\dot{\theta}_1\cos{\theta_1}\cos{\theta_2} + r_2\dot{r}_2\sin{\theta_2}\cos{\theta_2} + r_2^2\dot{\theta}_2\cos^2{\theta_2}\right]\\
	&-\dfrac{1}{2}(b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[-\dot{r}_1r_2\cos{\theta}_1\sin{\theta_2} + r_1r_2\dot{\theta}_1\sin{\theta_1}\sin{\theta_2} - \dfrac{r_2\dot{r}_2\cos{\theta_2}\sin{\theta_2} + r_2^2\dot{\theta}_2\sin^2{\theta_2}}{2} + \dot{r}_1r_2\sin{\theta_1}\cos{\theta_2} + r_1r_2\dot{\theta}_1\cos{\theta_1}\cos{\theta_2} + \dfrac{r_2\dot{r}_2\sin{\theta_2}\cos{\theta_2} + r_2^2\dot{\theta}_2\cos^2{\theta_2}}{2}\right]\\
	&= -(b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[\dot{r}_1r_2(-\cos{\theta}_1\sin{\theta_2} + \sin{\theta_1}\cos{\theta_2}) + r_1r_2\dot{\theta}_1(\sin{\theta_1}\sin{\theta_2}+\cos{\theta_1}\cos{\theta_2}) + r_2\dot{r}_2(-\cos{\theta_2}\sin{\theta_2} + \sin{\theta_2}\cos{\theta_2}) + r_2^2\dot{\theta}_2(\sin^2{\theta_2} +\cos^2{\theta_2})\right]\\
	&-\dfrac{1}{2}(b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dot{r}_1r_2(-\cos{\theta}_1\sin{\theta_2} + \sin{\theta_1}\cos{\theta_2}) + r_1r_2\dot{\theta}_1(\sin{\theta_1}\sin{\theta_2}+\cos{\theta_1}\cos{\theta_2}) + \dfrac{r_2\dot{r}_2(-\cos{\theta_2}\sin{\theta_2} + \sin{\theta_2}\cos{\theta_2}) + r_2^2\dot{\theta}_2(\sin^2{\theta_2} +\cos^2{\theta_2})}{2}\right]\\
	&= -(b_{2b}+c_{2b}|\vec{v}_{2b}|)\left[ r_2^2\dot{\theta}_2-\dot{r}_1r_2\sin{\Delta} + r_1r_2\dot{\theta}_1\cos{\Delta}\right]-\dfrac{1}{2}(b_{2r}+c_{2r}|\vec{v}_{2r}|)\left[\dfrac{r_2^2\dot{\theta}_2}{2}-\dot{r}_1r_2\sin{\Delta} + r_1r_2\dot{\theta}_1\cos{\Delta}\right].\\
\end{align*}

Hence Equation \eqref{ELD} for $\theta_2$ is

\begin{align*}
	M_2(r_2^2 \ddot{\theta}_2 + 2r_2 \dot{r}_2 \dot{\theta}_2) + \mu_2 (\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2)+gr_2\cos{\theta_2}) &= Q_{\theta_2}.
\end{align*}

Or, in matrix form

\begin{align*}
\begin{bmatrix}
-\mu_2r_2 \sin{\Delta} & 0 & \mu_2r_1 r_2\cos{\Delta} & M_2r_2^2
\end{bmatrix} \begin{bmatrix}
\ddot{r}_1\\
\ddot{r}_2\\
\ddot{\theta}_1\\
\ddot{\theta}_2\\
\end{bmatrix} &= \begin{bmatrix}
-2M_2 r_2\dot{r}_2 \dot{\theta}_2 - \mu_2(2\dot{r}_1r_2\dot{\theta}_1\cos{\Delta} + r_1r_2\dot{\theta}_1^2\sin{\Delta} + gr_2\cos{\theta_2}) + Q_{\theta_2}
\end{bmatrix}.
\end{align*}

## Matrix form
\begin{align*}
\begin{bmatrix}
M_1 & \mu_2 \cos{\Delta} & 0 & -\mu_2r_2 \sin{\Delta}\\
\mu_2 \cos{\Delta} & M_2 & \mu_2 r_1\sin{\Delta} & 0 \\
0 & \mu_2 r_1\sin{\Delta} & M_1r_1^2 & \mu_2 r_1r_2\cos{\Delta}\\
-\mu_2r_2 \sin{\Delta} & 0 & \mu_2r_1 r_2\cos{\Delta} & M_2r_2^2
\end{bmatrix} \begin{bmatrix}
\ddot{r}_1\\
\ddot{r}_2\\
\ddot{\theta}_1\\
\ddot{\theta}_2\\
\end{bmatrix} &= \begin{bmatrix}
M_1r_1\dot{\theta}_1^2 + \mu_2 (r_2 \dot{\theta}_2^2\cos{\Delta} + 2 \dot{r}_2 \dot{\theta}_2\sin{\Delta}) - \mu_1g\sin{\theta_1} - k_1(r_1-l_1) + Q_{r_1} \\
M_2r_2\dot{\theta}_2^2 + \mu_2 (r_1\dot{\theta}_1^2\cos{\Delta}-2\dot{r}_1\dot{\theta}_1\sin{\Delta} - g\sin{\theta_2}) - k_2(r_2-l_2) + Q_{r_2}\\
-2M_1r_1\dot{r}_1 \dot{\theta}_1 + \mu_2(r_1r_2\dot{\theta}_2^2 \sin{\Delta} -2r_1\dot{r}_2 \dot{\theta}_2\cos{\Delta}) - \mu_1 gr_1\cos{\theta_1} + Q_{\theta_1} \\
-2M_2 r_2\dot{r}_2 \dot{\theta}_2 - \mu_2(2\dot{r}_1r_2\dot{\theta}_1\cos{\Delta} + r_1r_2\dot{\theta}_1^2\sin{\Delta} + gr_2\cos{\theta_2}) + Q_{\theta_2}
\end{bmatrix}
\end{align*}