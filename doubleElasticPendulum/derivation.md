@def title="Deriving the equations of motion for the double elastic pendulum"

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
* The left-hand side of Equation Equation \eqref{ELD} can also be represented as $\dfrac{\delta \mathcal{L}}{\delta q_i}$, the functional derivative of the Lagrangian with respect to $q_i$. 
* The right-hand side of Equation Equation \eqref{ELD} is also called the generalized dissipative force and can be represented as $Q_i$.

The masses of the pendulum rods (or springs) are ignored as including them into the calculation for [rigid double pendulums](/doublePendulum/) does not make things more interesting. 
 
\tableofcontents

As can be seen, we have four degrees of freedom in this system. The angles the two pendulums make with the positive $x$-axis &mdash; $\theta_1$ and $\theta_2$, respectively &mdash; are among our degrees of freedom. We will also need degrees of freedom corresponding to the lengths of the pendulum rods. These degrees of freedom could either be the extent to which they are extended beyond their rest length or their total length. For the sake of simplicity, we will opt to use their total lengths &mdash; $r_1$ and $r_2$, respectively. Hence

\begin{align*}
	x_1 &= r_1 \cos{\theta_1} & \dot{x}_1 &= \dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1}\\
	y_1 &= r_1 \sin{\theta_1} & \dot{y}_1 &= \dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} \\
	x_2 &= x_1 + r_2\cos{\theta_2} & \dot{x}_2 &= \dot{x}_1 + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
	y_2 &= y_1 + r_2\sin{\theta_2} & \dot{y}_2 &= \dot{y}_1 + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}.
\end{align*}

This means that the velocity of the first pendulum bob is

\begin{align*}
	\vec{v}_1 &= \begin{bmatrix}
		\dot{x}_1 \\
		\dot{y}_1
	\end{bmatrix} \\
	&= \begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1}
	\end{bmatrix}.
\end{align*}

Hence 
\begin{align*}
	|\vec{v}_1|^2 &= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2. 
\end{align*}

As for the velocity of the second pendulum bob, it is
\begin{align*}
	\vec{v}_2 &= \begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix}.
\end{align*}

Let $\Delta = \theta_2-\theta_1$, then square of the velocity is

\begin{align*}
	|\vec{v}_2|^2 &= \dot{r}_1^2 \cos^2{\theta_1} + r_1^2 \dot{\theta}_1^2\sin^2{\theta_1} + \dot{r}_2^2\cos^2{\theta_2} + r_2^2\dot{\theta}_2^2\sin^2{\theta_2} -2r_1\dot{r}_1\dot{\theta}_1 \cos{\theta}_1\sin{\theta_1} + 2\dot{r}_1\dot{r}_2\cos{\theta_1}\cos{\theta_2} - 2\dot{r}_1r_2\dot{\theta}_2\cos{\theta_1}\sin{\theta_2} - 2r_1\dot{r}_2 \dot{\theta}_1 \sin{\theta_1}\cos{\theta_2} + 2r_1r_2 \dot{\theta}_1\dot{\theta}_2\sin{\theta_1}\sin{\theta_2}\\
	&-2\dot{r}_2r_2\dot{\theta}_2\cos{\theta_2}\sin{\theta_2} + \dot{r}_1^2\sin^2{\theta_1} + r_1^2\dot{\theta}_1^2\cos^2{\theta_1} + \dot{r}_2^2\sin^2{\theta_2} + r_2^2\dot{\theta}_2^2\cos^2{\theta_2} + 2r_1\dot{r}_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + 2\dot{r}_1\dot{r}_2\sin{\theta_1}\sin{\theta_2} + 2\dot{r}_1r_2\dot{\theta}_2\sin{\theta_1}\cos{\theta_2} + 2r_1\dot{r}_2\dot{\theta}_1 \cos{\theta_1}\sin{\theta_2}\\
	&+2r_1r_2\dot{\theta}_1\dot{\theta}_2\cos{\theta_1}\cos{\theta_2} + 2r_2\dot{r}_2\dot{\theta}_2\sin{\theta_2}\cos{\theta_2} \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2r_1\dot{r}_1\dot{\theta}_1(-\cos{\theta_1}\sin{\theta_1} + \cos{\theta_1}\sin{\theta_1}) + 2\dot{r}_1\dot{r}_2(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2})+2\dot{r}_1r_2\dot{\theta}_2(-\cos{\theta_1}\sin{\theta_2} + \sin{\theta_1}\cos{\theta_2}) \\
	&+ 2r_1\dot{r}_2\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_2}+\cos{\theta_1}\sin{\theta_2}) + 2r_1r_2\dot{\theta}_1\dot{\theta}_2(\sin{\theta_1}\sin{\theta_2} + \cos{\theta_1}\cos{\theta_2}) + 2r_2\dot{r}_2\dot{\theta}_2 (-\cos{\theta_2}\sin{\theta_2} + \sin{\theta_2}\cos{\theta_2}) \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2\dot{r}_1\dot{r}_2 \cos{(\theta_2-\theta_1)} - 2\dot{r}_1r_2\dot{\theta_2}\sin{(\theta_2-\theta_1)} + 2r_1\dot{r}_2 \dot{\theta}_1 \sin{(\theta_2-\theta_1)} + 2r_1r_2\dot{\theta}_1\dot{\theta}_2\cos{(\theta_2-\theta_1)} \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2\dot{r}_1\dot{r}_2 \cos{\Delta} - 2\dot{r}_1r_2\dot{\theta_2}\sin{\Delta} + 2r_1\dot{r}_2 \dot{\theta}_1 \sin{\Delta} + 2r_1r_2\dot{\theta}_1\dot{\theta}_2\cos{\Delta} \\
	&= \dot{r}_1^2 + r_1^2 \dot{\theta}_1^2 + \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2\cos{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\sin{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2).
\end{align*}

Let us define $|\Delta \vec{v}_{21}|^2 = |\vec{v}_2|^2 - |\vec{v}_1|^2$, as this will simplify our Lagrangian later.

\begin{align*}
	|\Delta \vec{v}_{21}|^2 &= \dot{r}_2^2 + r_2^2\dot{\theta}_2^2 + 2\cos{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\sin{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2).
\end{align*}

## Dissipative forces
We will assume that the dissipative forces are proportion to the velocity and velocity squared of the pendulum bobs. Meaning they will have the form

\begin{align*}
\vec{F}_{D,j} &= -(b_j+c_j|\vec{v}_j|)\vec{v}_j.
\end{align*}

Where $j$ is the pendulum bob of interest, $b_j$ and $c_j$ are constants. 

## Kinetic energy
The kinetic energy of the system is given by

\begin{align*}
	T &= \dfrac{m_1}{2}|\vec{v}_1|^2 + \dfrac{m_2}{2}|\vec{v}_2|^2 \\
	&= \dfrac{m_1+m_2}{2}|\vec{v}_1|^2 + \dfrac{m_2}{2}|\Delta \vec{v}_{21}|^2.
\end{align*}

## Potential energy
The potential energy of the system is given by

\begin{align*}
	V &= m_1 gy_1 + m_2gy_2 + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}\\
	&= m_1 gr_1\sin{\theta_1} + m_2g(r_1\sin{\theta_1} + r_2\sin{\theta_2}) + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}\\
	&= (m_1+m_2)gr_1\sin{\theta_1} + m_2gr_2\sin{\theta_2} + \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}.
\end{align*}

## Lagrangian
Hence the Lagrangian of the system is

\begin{align*}
	\mathcal{L} &= T - V \\
	&= \dfrac{m_1+m_2}{2}|\vec{v}_1|^2 + \dfrac{m_2}{2}|\Delta \vec{v}_{21}|^2 - \left[(m_1+m_2)gr_1\sin{\theta_1} + m_2gr_2\sin{\theta_2}\right] - \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2} \\
	&= (m_1+m_2)\left[\dfrac{|\vec{v}_1|^2}{2} - gr_1\sin{\theta_1}\right] + m_2\left[\dfrac{|\Delta \vec{v}_{21}|^2}{2} - gr_2\sin{\theta_2}\right]  - \dfrac{k_1(r_1-l_1)^2+k_2(r_2-l_2)^2}{2}.
\end{align*}

We will not expand this Lagrangian, as doing so just adds to its complexity. Instead, we will calculate the derivatives of each of its components. 

## Derivative of components of the Lagrangian
### Square of the velocity of the first pendulum's bob
The relevant partial and standard derivatives are:
\begin{align*}
	\dfrac{\partial |\vec{v}_1|^2}{\partial r_1} &= 2r_1\dot{\theta}_1^2 & \dfrac{\partial |\vec{v}_1|^2}{\partial r_2} &= 0 & \dfrac{\partial |\vec{v}_1|^2}{\partial \theta_1} &= 0 & \dfrac{\partial |\vec{v}_1|^2}{\partial \theta_2} &= 0\\
	\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{r}_1} &= 2\dot{r}_1 & \dfrac{\partial |\vec{v}_1|^2}{\partial \dot{r}_2} &= 0 & \dfrac{\partial |\vec{v}_1|^2}{\partial \dot{\theta}_1} &= 2r_1^2\dot{\theta}_1 & \dfrac{\partial |\vec{v}_1|^2}{\partial \dot{\theta}_2} &= 0\\
	\dfrac{d}{dt} \dfrac{\partial |\vec{v}_1|^2}{\partial \dot{r}_1} &= 2\ddot{r}_1 & \dfrac{d}{dt}\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{r}_2} &= 0 & \dfrac{d}{dt}\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{\theta}_1} &= 2r_1^2 \ddot{\theta}_1 + 4r_1\dot{r}_1\dot{\theta}_1 & \dfrac{d}{dt}\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{\theta}_2} &= 0
\end{align*}

Hence the functional derivatives are

\begin{align*}
	\dfrac{\delta |\vec{v}_1|^2}{\delta r_1} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{r}_1} - \dfrac{\partial |\vec{v}_1|^2}{\partial r_1} & \dfrac{\delta |\vec{v}_1|^2}{\delta r_2} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{r}_2} - \dfrac{\partial |\vec{v}_1|^2}{\partial r_2} & \dfrac{\delta |\vec{v}_1|^2}{\delta \theta_1} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{\theta}_1} - \dfrac{\partial |\vec{v}_1|^2}{\partial \theta_1} & \dfrac{\delta |\vec{v}_1|^2}{\delta \theta_2} &= \dfrac{d}{dt}\dfrac{\partial |\vec{v}_1|^2}{\partial \dot{\theta}_2} - \dfrac{\partial |\vec{v}_1|^2}{\partial \theta_2}\\
	&= 2\ddot{r}_1 - 2r_1\dot{\theta}_1^2 & &= 0 & &=2r_1^2\ddot{\theta}_1 + 4r_1\dot{r}_1\dot{\theta}_1 & &= 0.	
\end{align*}

### Difference in the square of each bob's velocity
Hence the partial and standard derivatives of the difference in the square of each bob's velocity is
\begin{align*}
	\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial r_1} &= 2r_2\dot{\theta}_1\dot{\theta}_2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_1\sin{\Delta} & \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial r_2} &= 2r_2\dot{\theta}_2^2+2r_1\dot{\theta}_1\dot{\theta}_2\cos{\Delta}-2\dot{r}_1\dot{\theta}_2\sin{\Delta} \\
	\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \theta_1} &= -2\sin{\Delta}\cdot -1(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\cos{\Delta}\cdot -1(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2) & \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \theta_2} &= -2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)  \\
	&= 2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) - 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)\\
	\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{r}_1} &= 2\dot{r}_2\cos{\Delta} - 2r_2\dot{\theta}_2\sin{\Delta} & \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{r}_2} &= 2\dot{r}_2 + 2\dot{r}_1\cos{\Delta} + 2r_1\dot{\theta}_1\sin{\Delta} \\
	\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{\theta}_1} &= 2r_1r_2\dot{\theta}_2\cos{\Delta} + 2r_1\dot{r}_2\sin{\Delta} & \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{\theta}_2} &=2r_2^2\dot{\theta}_2 + 2r_1r_2\dot{\theta}_1\cos{\Delta} - 2\dot{r}_1r_2\sin{\Delta}
\end{align*}

Let us define $\dot{\Delta}_1 = 2\dot{\theta}_1 - \dot{\theta}_2$ and $\dot{\Delta}_2 = 2\dot{\theta}_2 - \dot{\theta}_1$.

\begin{align*}
	\dfrac{d}{dt} \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{r}_1} &= 2\ddot{r}_2\cos{\Delta} - 2\dot{r}_2\dot{\Delta}\sin{\Delta} - 2\dot{r}_2\dot{\theta}_2\sin{\Delta} - 2r_2\ddot{\theta}_2\sin{\Delta} - 2r_2\dot{\theta}_2\dot{\Delta}\cos{\Delta} \\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2\dot{\Delta} + \dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)\\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2(2\dot{\theta}_2-\dot{\theta}_1)+r_2\ddot{\theta}_2)\\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2\dot{\Delta}_2+r_2\ddot{\theta}_2)\\
	\dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{r}_2} &= 2\ddot{r}_2 + 2\ddot{r}_1 \cos{\Delta} -2\dot{r}_1\dot{\Delta}\sin{\Delta} + 2\dot{r}_1\dot{\theta}_1\sin{\Delta} + 2r_1\ddot{\theta}_1\sin{\Delta} + 2r_1\dot{\theta}_1\dot{\Delta}\cos{\Delta} \\
	&= 2\ddot{r}_2 + 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1\dot{\theta}_1 - \dot{r}_1\dot{\Delta})\\
	&= 2\ddot{r}_2 + 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1(2\dot{\theta}_1 - \ddot{\theta}_2))\\
	&= 2\ddot{r}_2 + 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1\dot{\Delta}_1)\\
	\dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{\theta}_1} &= 2\dot{r}_1r_2\dot{\theta}_2\cos{\Delta} + 2r_1\dot{r}_2\dot{\theta}_2\cos{\Delta} + 2r_1r_2\ddot{\theta}_2\cos{\Delta} - 2r_1r_2\dot{\theta}_2\dot{\Delta}\sin{\Delta} + 2\dot{r}_1\dot{r}_2\sin{\Delta} + 2r_1\ddot{r}_2\sin{\Delta} + 2r_1\dot{r}_2\dot{\Delta}\cos{\Delta} \\
	&= 2\cos{\Delta}(\dot{r}_1r_2\dot{\theta}_2 + r_1\dot{r}_2 (\dot{\theta}_2+\dot{\Delta})+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2\dot{\Delta})\\
	&= 2\cos{\Delta}(\dot{r}_1r_2\dot{\theta}_2 + r_1\dot{r}_2 \dot{\Delta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2\dot{\Delta})\\
	\dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{\theta}_2} &= 4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 + 2\dot{r}_1r_2\dot{\theta}_1\cos{\Delta} + 2r_1\dot{r}_2\dot{\theta}_1\cos{\Delta} + 2r_1r_2\ddot{\theta}_1\cos{\Delta} - 2r_1r_2\dot{\theta}_1\dot{\Delta}\sin{\Delta} - 2\ddot{r}_1r_2\sin{\Delta} - 2\dot{r}_1\dot{r}_2\sin{\Delta} - 2\dot{r}_1r_2\dot{\Delta}\cos{\Delta} \\
	&=4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 +2\cos{\Delta}(\dot{r}_1r_2(\dot{\theta}_1-\dot{\Delta})+r_1\dot{r}_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)-2\sin{\Delta}(r_1r_2\dot{\theta}_1\dot{\Delta} + \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2) \\
	&=4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 +2\cos{\Delta}(\dot{r}_1r_2\dot{\Delta}_1+r_1\dot{r}_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)-2\sin{\Delta}(r_1r_2\dot{\theta}_1\dot{\Delta} + \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2).
\end{align*}

Hence functional derivative for $r_1$ is
\begin{align*}
	\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta r_1} &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{r}_1} - \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial r_1} \\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2\dot{\Delta}) - 2\sin{\Delta}(\dot{r}_2\dot{\Delta}_2+r_2\ddot{\theta}_2) - \left(2r_2\dot{\theta}_1\dot{\theta}_2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_1\sin{\Delta}\right) \\
	&= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2(\dot{\Delta}+\dot{\theta}_1)) - 2\sin{\Delta} (\dot{r}_2(\dot{\Delta}_2+\dot{\theta}_1)+r_2\ddot{\theta}_2).
\end{align*}

Where $\dot{\Delta} + \dot{\theta}_1 = \dot{\theta}_2 - \dot{\theta}_1 + \dot{\theta}_1 = \dot{\theta}_2$ and $\dot{\Delta}_2 + \dot{\theta}_1 = 2\dot{\theta}_2 - \dot{\theta}_1 + \dot{\theta}_1 = 2\dot{\theta}_2$.

\begin{align*}
	\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta r_1} &= 2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - 2\sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2).
\end{align*}

As for $r_2$

\begin{align*}
	\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta r_2} &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{r}_2} - \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial r_2} \\
	&= 2\ddot{r}_2 + 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1 \dot{\Delta}) + 2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1\dot{\Delta}_1) - \left[2r_2\dot{\theta}_2^2+2r_1\dot{\theta}_1\dot{\theta}_2\cos{\Delta}-2\dot{r}_1\dot{\theta}_2\sin{\Delta}\right]\\
	&= 2\ddot{r}_2 -2r_2\dot{\theta}_2^2 + 2\cos{\Delta}(\ddot{r}_1 + r_1\dot{\theta}_1( \dot{\Delta}-\dot{\theta}_2))+2\sin{\Delta}(r_1\ddot{\theta}_1 + \dot{r}_1(\dot{\Delta}_1+\dot{\theta}_2)).
\end{align*}

Hence $\dot{\Delta}-\dot{\theta}_2 = \dot{\theta}_2-\dot{\theta}_1-\dot{\theta}_2 = -\dot{\theta}_1$ and $\dot{\Delta}_1 + \dot{\theta}_2 = 2\dot{\theta}_1 - \dot{\theta}_2 + \dot{\theta}_2 = 2\dot{\theta}_1$. 

\begin{align*}
	\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta r_2} &= 2\ddot{r}_2 -2r_2\dot{\theta}_2^2 + 2\cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+2\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1).
\end{align*}

As for $\theta_1$

\begin{align*}
	\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta \theta_1} &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{\theta}_1} - \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \theta_1} \\
	&= 2\cos{\Delta}(\dot{r}_1r_2\dot{\theta}_2 + r_1\dot{r}_2 \dot{\Delta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2\dot{\Delta}) - \left[2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) - 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)\right] \\
	&= 2\cos{\Delta}(\dot{r}_1r_2(\dot{\theta}_2 - \dot{\theta}_2)+ r_1\dot{r}_2 (\dot{\Delta}_2+\dot{\theta}_1)+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(\dot{r}_1\dot{r}_2-\dot{r}_1\dot{r}_2+r_1\ddot{r}_2-r_1r_2\dot{\theta}_2(\dot{\Delta}+\dot{\theta}_1)) \\
	&= 2\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2).
\end{align*}

As for $\theta_2$

\begin{align*}
	\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta \theta_2} &= \dfrac{d}{dt}\dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \dot{\theta}_2} - \dfrac{\partial |\Delta \vec{v}_{21}|^2}{\partial \theta_2} \\
	&= 4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 +2\cos{\Delta}(\dot{r}_1r_2\dot{\Delta}_1+r_1\dot{r}_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)-2\sin{\Delta}(r_1r_2\dot{\theta}_1\dot{\Delta} + \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2) - \left[-2\sin{\Delta}(\dot{r}_1\dot{r}_2 + r_1r_2\dot{\theta}_1\dot{\theta}_2) + 2\cos{\Delta}(r_1\dot{r}_2\dot{\theta}_1-\dot{r}_1r_2\dot{\theta}_2)\right] \\
	&= 4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 +2\cos{\Delta}(\dot{r}_1r_2(\dot{\Delta}_1+\dot{\theta}_2)+r_1\dot{r}_2(\dot{\theta}_1-\dot{\theta}_1)+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1(\dot{\theta}_2-\dot{\Delta}) - \ddot{r}_1r_2 + \dot{r}_1\dot{r}_2-\dot{r}_1\dot{r}_2) \\
	&= 4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 +2\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1(\dot{\theta}_2-(\dot{\theta}_2-\dot{\theta}_1))-\ddot{r}_1r_2) \\
	&= 4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 +2\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2).
\end{align*}

## Euler-Lagrange equations with dissipation
### $r_1$
It is important to note that $\dfrac{\delta f(q_i)}{\delta q_i} = -\dfrac{\partial f}{\partial q_i}$ and of course if a term does not depend on $q_i$ or $\dot{q}_i$ its functional derivative with respect to $q_i$ is zero. Hence

\begin{align*}
	\dfrac{\delta \mathcal{L}}{\delta r_1} &= (m_1+m_2)\left[\dfrac{\dfrac{\delta |\vec{v}_1|^2}{\delta r_1}}{2} -g\sin{\theta_1}\dfrac{\delta r_1}{\delta r_1}\right] + m_2\left[\dfrac{\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta r_1}}{2}\right] - \dfrac{k_1}{2} \dfrac{\delta (r_1-l_1)^2}{\delta r_1}.
\end{align*}
We have deliberately ignored the $m_2gr_2\sin{\theta_2}$ and $-\dfrac{k_2(r_2-l_2)^2}{2}$ as they are independent of $r_1$.
\begin{align*}
	\dfrac{\delta \mathcal{L}}{\delta r_1} &= (m_1+m_2)\left[\dfrac{2\ddot{r}_1-2r_1\dot{\theta}_1^2}{2} + g\sin{\theta_1}\right] + m_2\left[\dfrac{2\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - 2\sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)}{2}\right] + k_1(r_1-l_1)\\
	&= (m_1+m_2)\left[\ddot{r}_1-r_1\dot{\theta}_1^2 + g\sin{\theta_1}\right] + m_2\left[\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - \sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)\right] + k_1(r_1-l_1).
\end{align*}
	
The generalized dissipation force canonical to $r_1$ is hence
\begin{align*}
	Q_{r_1} &= -(b_1+c_1|\vec{v}_1|)\vec{v}_1\cdot \dfrac{\partial \vec{r}_1}{\partial r_1} - (b_2+c_2|\vec{v}_2|)\vec{v}_2\cdot \dfrac{\partial \vec{r}_2}{\partial r_1} \\
	&= -(b_1+c_1|\vec{v}_1|)\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}
	\end{bmatrix} \cdot \begin{bmatrix}
		\cos{\theta_1} \\
		\sin{\theta_1}
	\end{bmatrix} - (b_2+c_2|\vec{v}_2|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot \begin{bmatrix}
		\cos{\theta_1} \\
		\sin{\theta_1}
	\end{bmatrix} \\
	&= -(b_1+c_1|\vec{v}_1|)\left[\dot{r}_1\cos^2{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + \dot{r}_1\sin^2{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_1}\right] - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1\cos^2{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_1} + \dot{r}_2\cos{\theta_1}\cos{\theta_2}-r_2\dot{\theta}_2\cos{\theta_1}\sin{\theta_2} \right.\\
	&\left.+ \dot{r}_1\sin^2{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_1} + \dot{r}_2\sin{\theta_1}\sin{\theta_2} + r_2\dot{\theta}_2\sin{\theta_1}\cos{\theta_2} \right] \\
	&= -(b_1+c_1|\vec{v}_1|)\left[\dot{r}_1(\cos^2{\theta_1}+\sin^2{\theta_1})+r_1\dot{\theta}_1(-\sin{\theta}_1\cos{\theta_1}+\sin{\theta_1}\cos{\theta_1})\right] - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1(\cos^2{\theta_1}+\sin^2{\theta_1}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_1}+\sin{\theta_1}\cos{\theta_1}) \right.\\
	&\left. +  \dot{r}_2(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2})+r_2\dot{\theta}_2(-\cos{\theta_1}\sin{\theta_2} + \sin{\theta_1}\cos{\theta_2}) \right] \\
	&= -(b_1+c_1|\vec{v}_1|)\dot{r}_1 - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1+\dot{r}_2\cos{(\theta_2-\theta_1)}-r_2\dot{\theta}_2\sin{(\theta_2-\theta_1)} \right]. \\
	&= -(b_1+c_1|\vec{v}_1|)\dot{r}_1 - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1+\dot{r}_2\cos{\Delta}-r_2\dot{\theta}_2\sin{\Delta} \right]. \\
\end{align*}

Hence the Euler-Lagrange equation for $r_1$ with dissipative forces is

\begin{align*}
	(m_1+m_2)\left[\ddot{r}_1-r_1\dot{\theta}_1^2 + g\sin{\theta_1}\right] + m_2\left[\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - \sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)\right] + k_1(r_1-l_1) &= Q_{r_1}.
\end{align*}

Dividing by $m_1+m_2$

\begin{align*}
	\ddot{r}_1-r_1\dot{\theta}_1^2 + g\sin{\theta_1} + \dfrac{m_2}{m_1+m_2}\left[\cos{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2) - \sin{\Delta} (2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2)\right] + \dfrac{k_1}{m_1+m_2}(r_1-l_1) &= \dfrac{Q_{r_1}}{m_1+m_2}.
\end{align*}

Expanding out second derivative terms and placing them first on the left-hand side yields
\begin{align*}
	\ddot{r}_1 + \dfrac{m_2\cos{\Delta}}{m_1+m_2}\ddot{r}_2 + 0\ddot{\theta}_1 - \dfrac{m_2r_2\sin{\Delta}}{m_1+m_2}\ddot{\theta}_2 - r_1\dot{\theta}_1^2 + g\sin{\theta_1} - \dfrac{m_2}{m_1+m_2}\left[r_2\dot{\theta}_2^2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_2\sin{\Delta}\right] + \dfrac{k_1}{m_1+m_2}(r_1-l_1) &= \dfrac{Q_{r_1}}{m_1+m_2}.
\end{align*}

Moving all terms that do not involve second derivatives to the right-hand side yields

\begin{align*}
	\ddot{r}_1 + \dfrac{m_2\cos{\Delta}}{m_1+m_2}\ddot{r}_2 + 0\ddot{\theta}_1 - \dfrac{m_2r_2\sin{\Delta}}{m_1+m_2}\ddot{\theta}_2 &= r_1\dot{\theta}_1^2 - g\sin{\theta_1} + \dfrac{m_2}{m_1+m_2}\left[r_2\dot{\theta}_2^2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_2\sin{\Delta}\right] + \dfrac{Q_{r_1}-k_1(r_1-l_1)}{m_1+m_2}.
\end{align*}

### $r_2$
As for $r_2$

\begin{align*}
	\dfrac{\delta \mathcal{L}}{\delta r_2} &= (m_1+m_2)\left[\dfrac{\dfrac{\delta |\vec{v}_1|^2}{\delta r_2}}{2}\right] + m_2\left[\dfrac{\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta r_2}}{2} - g\sin{\theta}_2\dfrac{\delta r_2}{\delta r_2}\right] - \dfrac{k_2}{2} \dfrac{\delta (r_2-l_2)^2}{\delta r_2} \\
	&= (m_1+m_2)\left[\dfrac{0}{2}\right] + m_2\left[\dfrac{2\ddot{r}_2 -2r_2\dot{\theta}_2^2 + 2\cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+2\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1)}{2} + g\sin{\theta_2}\right] + k_2(r_2-l_2) \\
	&= m_2\left[\ddot{r}_2 - r_2\dot{\theta}_2^2 + \cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1) + g\sin{\theta_2}\right] + k_2(r_2-l_2)\\
	Q_{r_2} &= -(b_1+c_1|\vec{v}_1|)\vec{v}_1 \cdot \dfrac{\partial \vec{r}_1}{\partial r_2} - (b_2+c_2|\vec{v}_2|)\vec{v}_2 \cdot \dfrac{\partial \vec{r}_2}{\partial r_2} \\
	&= -(b_1+c_1|\vec{v}_1|)\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta_1}
	\end{bmatrix} \cdot \vec{0}  - (b_2+c_2|\vec{v}_2|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot \begin{bmatrix}
		\cos{\theta_2} \\
		\sin{\theta_2}
	\end{bmatrix} \\
	&=  - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1\cos{\theta_1}\cos{\theta_2} - r_1\dot{\theta}_1\sin{\theta_1}\cos{\theta_2} + \dot{r}_2\cos^2{\theta_2}-r_2\dot{\theta}_2\sin{\theta_2}\cos{\theta_2} + \dot{r}_1\sin{\theta_1}\sin{\theta_2} + r_1\dot{\theta}_1\cos{\theta_1}\sin{\theta_2} + \dot{r}_2\sin^2{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}\sin{\theta_2}\right] \\
	&= - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1(\cos{\theta_1}\cos{\theta_2}+\sin{\theta_1}\sin{\theta_2}) + r_1\dot{\theta}_1(-\sin{\theta_1}\cos{\theta_2} + \cos{\theta_1}\sin{\theta_2}) + \dot{r}_2(\cos^2{\theta_2} + \sin^2{\theta_2})+r_2\dot{\theta}_2(-\sin{\theta_2}\cos{\theta_2} +  \cos{\theta_2}\sin{\theta_2})\right] \\
	&= - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1\cos{(\theta_2-\theta_1)} + r_1\dot{\theta}_1\sin{(\theta_2-\theta_1)} + \dot{r}_2\right] \\
	&= - (b_2+c_2|\vec{v}_2|)\left[\dot{r}_1\cos{\Delta} + r_1\dot{\theta}_1\sin{\Delta} + \dot{r}_2\right].
\end{align*}

Hence the Euler-Lagrange equation for $r_2$ with dissipative forces is

\begin{align*}
	m_2\left[\ddot{r}_2 - r_2\dot{\theta}_2^2 + \cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1) + g\sin{\theta_2}\right] + k_2(r_2-l_2) &= Q_{r_2}
\end{align*}

Dividing by $m_2$ yields

\begin{align*}
	\ddot{r}_2 - r_2\dot{\theta}_2^2 + \cos{\Delta}(\ddot{r}_1 - r_1\dot{\theta}_1^2)+\sin{\Delta}(r_1\ddot{\theta}_1 + 2\dot{r}_1\dot{\theta}_1) + g\sin{\theta_2} + \dfrac{k_2(r_2-l_2)}{m_2} &= \dfrac{Q_{r_2}}{m_2}.
\end{align*}

Next we will expand out second time derivatives and moving everything else to the right-hand side
\begin{align*}
	\cos{\Delta}\ddot{r}_1 + \ddot{r}_2 + r_1\sin{\Delta}\ddot{\theta}_1 + 0\ddot{\theta}_2 &= r_2\dot{\theta}_2^2  - g\sin{\theta_2} + r_1\dot{\theta}_1^2\cos{\Delta} - 2\dot{r}_1\dot{\theta}_1\sin{\Delta} + \dfrac{Q_{r_2}-k_2(r_2-l_2)}{m_2}.
\end{align*}

### $\theta_1$
As for $\theta_1$

\begin{align*}
	\dfrac{\delta \mathcal{L}}{\delta \theta_1} &= (m_1+m_2)\left[\dfrac{\dfrac{\delta |\vec{v}_1|^2}{\delta \theta_1}}{2} - gr_1\dfrac{\delta \sin{\theta_1}}{\delta \theta_1}\right] + m_2\left[\dfrac{\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta \theta_1}}{2}\right] \\
	&= (m_1+m_2)\left[\dfrac{2r_1^2\ddot{\theta}_1 + 4r_1\dot{r}_1\dot{\theta}_1}{2} + gr_1\cos{\theta}_1\right] + m_2\left[\dfrac{2\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +2\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2)}{2}\right]\\
	&= (m_1+m_2)\left[r_1^2\ddot{\theta}_1 + 2r_1\dot{r}_1\dot{\theta}_1 + gr_1\cos{\theta}_1\right] + m_2\left[\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2)\right]\\
	Q_{\theta_1} &= -(b_1+c_1|\vec{v}_1|)\vec{v}_1 \cdot \dfrac{\partial \vec{r}_1}{\partial \theta_1} - (b_2+c_2|\vec{v}_2|)\vec{v}_2 \cdot \dfrac{\partial \vec{r}_2}{\partial \theta_1} \\
	&=  -(b_1+c_1|\vec{v}_1|)\begin{bmatrix}
		\dot{r}_1\cos{\theta_1} - r_1\dot{\theta}_1\sin{\theta_1} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1\cos{\theta}_1
	\end{bmatrix} \cdot r_1\begin{bmatrix}
		-\sin{\theta_1} \\
		\cos{\theta_1}
	\end{bmatrix} - (b_2+c_2|\vec{v}_2|) \begin{bmatrix}
	\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
	\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot r_1\begin{bmatrix}
	-\sin{\theta_1} \\
	\cos{\theta_1}
	\end{bmatrix} \\
	&= -(b_1+c_1|\vec{v}_1|)\left[-r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1}+r_1^2\dot{\theta}_1\sin^2{\theta_1} + r_1\dot{r}_1\sin{\theta_1}\cos{\theta_1}+r_1^2\dot{\theta}_1\cos^2{\theta_1}\right] - (b_2+c_2|\vec{v}_2|)\left[-r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1} + r_1^2\dot{\theta_1}\sin^2{\theta_1} -r_1\dot{r}_2\sin{\theta_1}\cos{\theta_2} \right. \\
	&\left.+r_1r_2\dot{\theta}_2\sin{\theta_1}\sin{\theta_2}+r_1\dot{r}_1\cos{\theta_1}\sin{\theta_1}+r_1^2\dot{\theta}_1\cos^2{\theta_1} + r_1\dot{r}_2\cos{\theta_1}\sin{\theta_2}+r_1r_2\dot{\theta}_2\cos{\theta_1}\cos{\theta_2}\right] \\
	&=  -(b_1+c_1|\vec{v}_1|)\left[r_1\dot{r}_1(-\cos{\theta_1}\sin{\theta_1}+\sin{\theta_1}\cos{\theta_1}) + r_1^2\dot{\theta}_1(\sin^2{\theta_1} +\cos^2{\theta_1})\right] - (b_2+c_2|\vec{v}_2|)\left[r_1\dot{r}_1(-\cos{\theta_1}\sin{\theta_1} + \cos{\theta_1}\sin{\theta_1}) +r_1^2\dot{\theta_1}(\sin^2{\theta_1}+\cos^2{\theta_1})\right.\\
	&\left.+r_1\dot{r}_2(-\sin{\theta_1}\cos{\theta_2}+\cos{\theta_1}\sin{\theta_2})+r_1r_2\dot{\theta}_2(\sin{\theta_1}\sin{\theta_2}+\cos{\theta_1}\cos{\theta_2})\right] \\
	&= -(b_1+c_1|\vec{v}_1|)r_1^2\dot{\theta}_1 - (b_2+c_2|\vec{v}_2|)\left[ r_1^2\dot{\theta_1}+r_1\dot{r}_2\sin{(\theta_2-\theta_1)}+r_1r_2\dot{\theta}_2\cos{(\theta_2-\theta_1)}\right]
\end{align*}

Hence Equation \eqref{ELD} is

\begin{align*}
	(m_1+m_2)\left[r_1^2\ddot{\theta}_1 + 2r_1\dot{r}_1\dot{\theta}_1 + gr_1\cos{\theta}_1\right] + m_2\left[\cos{\Delta}(2r_1\dot{r}_2\dot{\theta}_2+r_1r_2\ddot{\theta}_2) +\sin{\Delta}(r_1\ddot{r}_2-r_1r_2\dot{\theta}_2^2)\right] &= Q_{\theta_1}.
\end{align*}

Dividing by $(m_1+m_2)r_1^2$ yields

\begin{align*}
	\ddot{\theta}_1 + \dfrac{2\dot{r}_1\dot{\theta}_1}{r_1} + \dfrac{g\cos{\theta_1}}{r_1} + \dfrac{m_2}{(m_1+m_2)r_1}\left[\cos{\Delta}(2\dot{r}_2\dot{\theta}_2+r_2\ddot{\theta}_2) +\sin{\Delta}(\ddot{r}_2-r_2\dot{\theta}_2^2)\right] &= \dfrac{Q_{\theta_1}}{(m_1+m_2)r_1^2}.
\end{align*}

Expanding out all second time derivatives and moving all other terms to the right-hand side yields

\begin{align*}
	0\ddot{r}_1 + \dfrac{m_2\sin{\Delta}}{(m_1+m_2)r_1}\ddot{r}_2 + \ddot{\theta}_1 + \dfrac{m_2r_2\cos{\Delta}}{(m_1+m_2)r_1}\ddot{\theta}_2 &= -\dfrac{2\dot{r}_1\dot{\theta}_1}{r_1} - \dfrac{g\cos{\theta_1}}{r_1} - \dfrac{m_2}{(m_1+m_2)r_1}\left[2\dot{r}_2\dot{\theta}_2\cos{\Delta} -r_2\dot{\theta}_2^2\sin{\Delta}\right] + \dfrac{Q_{\theta_1}}{(m_1+m_2)r_1^2}.
\end{align*}

### $\theta_2$
As for $\theta_2$

\begin{align*}
	\dfrac{\delta \mathcal{L}}{\delta \theta_2} &= (m_1+m_2)\left[\dfrac{\dfrac{\delta |\vec{v}_1|^2}{\delta \theta_2}}{2}\right] + m_2\left[\dfrac{\dfrac{\delta |\Delta \vec{v}_{21}|^2}{\delta \theta_2}}{2} - gr_2\dfrac{\delta \sin{\theta_2}}{\delta \theta_2}\right] \\
	&= (m_1+m_2)\left[\dfrac{0}{2}\right]+m_2\left[\dfrac{4r_2\dot{r}_2\dot{\theta}_2 + 2r_2^2\ddot{\theta}_2 +2\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+2\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2)}{2} + gr_2\cos{\theta_2}\right]\\
	&= m_2\left[2r_2\dot{r}_2\dot{\theta}_2 + r_2^2\ddot{\theta}_2 +\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2) + gr_2\cos{\theta_2}\right]\\
	Q_{\theta_2} &= -(b_1+c_1|\vec{v}_1|)\vec{v}_1 \cdot \dfrac{\partial \vec{r}_1}{\partial \theta_2} - (b_2+c_2|\vec{v}_2|)\vec{v}_2\cdot \dfrac{\partial \vec{r}_2}{\partial \theta_2} \\
	&= -(b_2+c_2|\vec{v}_2|)\begin{bmatrix}
		\dot{r}_1 \cos{\theta_1} - r_1 \dot{\theta}_1 \sin{\theta_1} + \dot{r}_2\cos{\theta_2} - r_2\dot{\theta}_2 \sin{\theta_2} \\
		\dot{r}_1\sin{\theta_1} + r_1\dot{\theta}_1 \cos{\theta_1} + \dot{r}_2\sin{\theta_2} + r_2\dot{\theta}_2 \cos{\theta_2}
	\end{bmatrix} \cdot r_2\begin{bmatrix}
		-\sin{\theta_2}\\
		\cos{\theta_2}
	\end{bmatrix} \\
	&= -(b_2+c_2|\vec{v}_2|)\left[-\dot{r}_1r_2\cos{\theta}_1\sin{\theta_2} + r_1r_2\dot{\theta}_1\sin{\theta_1}\sin{\theta_2} - r_2\dot{r}_2\cos{\theta_2}\sin{\theta_2} + r_2^2\dot{\theta}_2\sin^2{\theta_2} + \dot{r}_1r_2\sin{\theta_1}\cos{\theta_2} + r_1r_2\dot{\theta}_1\cos{\theta_1}\cos{\theta_2} + r_2\dot{r}_2\sin{\theta_2}\cos{\theta_2} + r_2^2\dot{\theta}_2\cos^2{\theta_2}\right]\\
	&= -(b_2+c_2|\vec{v}_2|)\left[\dot{r}_1r_2(-\cos{\theta}_1\sin{\theta_2} + \sin{\theta_1}\cos{\theta_2}) + r_1r_2\dot{\theta}_1(\sin{\theta_1}\sin{\theta_2}+\cos{\theta_1}\cos{\theta_2}) + r_2\dot{r}_2(-\cos{\theta_2}\sin{\theta_2} + \sin{\theta_2}\cos{\theta_2}) + r_2^2\dot{\theta}_2(\sin^2{\theta_2} +\cos^2{\theta_2})\right]\\
	&= -(b_2+c_2|\vec{v}_2|)\left[ r_2^2\dot{\theta}_2-\dot{r}_1r_2\sin{(\theta_2-\theta_1)} + r_1r_2\dot{\theta}_1\cos{(\theta_2-\theta_1)}\right].
\end{align*}

Hence Equation \eqref{ELD} for $\theta_2$ is

\begin{align*}
	m_2\left[2r_2\dot{r}_2\dot{\theta}_2 + r_2^2\ddot{\theta}_2 +\cos{\Delta}(2\dot{r}_1r_2\dot{\theta}_1+r_1r_2\ddot{\theta}_1)+\sin{\Delta}(r_1r_2\dot{\theta}_1^2-\ddot{r}_1r_2) + gr_2\cos{\theta_2}\right] &= Q_{\theta_2}.
\end{align*}

Dividing by $m_2r_2^2$ yields

\begin{align*}
	\ddot{\theta}_2 + \dfrac{2\dot{r}_2\dot{\theta}_2}{r_2} + \dfrac{\cos{\Delta}}{r_2} (2\dot{r}_1\dot{\theta}_1+r_1\ddot{\theta}_1)+\dfrac{\sin{\Delta}}{r_2}(r_1\dot{\theta}_1^2-\ddot{r}_1) + \dfrac{g\cos{\theta_2}}{r_2} &= \dfrac{Q_{\theta_2}}{m_2r_2^2}.
\end{align*}

Next we will expand out second time derivatives and moving everything else to the right-hand side

\begin{align*}
	-\dfrac{\sin{\Delta}}{r_2}\ddot{r}_1 + 0\ddot{r}_2 + \dfrac{r_1\cos{\Delta}}{r_2}\ddot{\theta}_1 + \ddot{\theta}_2 &= -\dfrac{2\dot{r}_2\dot{\theta}_2}{r_2} - \dfrac{g\cos{\theta_2}}{r_2} - \dfrac{2\dot{r}_1\dot{\theta}_1\cos{\Delta}}{r_2}-\dfrac{r_1\dot{\theta}_1^2\sin{\Delta}}{r_2} + \dfrac{Q_{\theta_2}}{m_2r_2^2}.
\end{align*}

## Analysis
There are three ways we could solve this problem; each of which involves numerical integration to obtain the final solution. Firstly, we could algebraically manipulate this ordinary differential equation (ODE) system until only one second time derivative appears in each equation. The final ODE system after this manipulation could, in turn, be numerically integrated using any standard scheme (e.g. the Runge-Kutta-Fehlberg method). Secondly, we could numerically integrate it as is using a differential-algebraic equation (DAE) solver, but these tend to be more prone to give convergence errors in my experience. Finally, we could convert the system to matrix form and invert it to obtain numerical approximations of each of our second time derivatives and use these to numerically integrate the system with a standard ODE solver. We will opt for this last method, as the first approach is almost guaranteed to introduce errors and the second gives convergence errors, at least in Julia. 

Essentially, we will write our differential equation system as
\begin{align*}
	\mathbf{A}\mathbf{\ddot{q}} &= \mathbf{b}.
\end{align*}
where $\mathbf{\ddot{q}}$ is a vector containing the second time derivatives of our generalized coordinates. Hence $\mathbf{\ddot{q}}=\mathbf{A}^{-1} \mathbf{b}$. Or, in full form, this above equation is

\begin{align*}
	\begin{bmatrix}
		1 & \dfrac{m_2\cos{\Delta}}{m_1+m_2} & 0 & -\dfrac{m_2r_2\sin{\Delta}}{m_1+m_2} \\
		\cos{\Delta} & 1 & r_1\sin{\Delta} & 0 \\
		0 & \dfrac{m_2\sin{\Delta}}{(m_1+m_2)r_1} & 1 & \dfrac{m_2r_2\cos{\Delta}}{(m_1+m_2)r_1} \\
		-\dfrac{\sin{\Delta}}{r_2} & 0 & \dfrac{r_1\cos{\Delta}}{r_2} & 1
	\end{bmatrix} \begin{bmatrix}
	\ddot{r}_1 \\
	\ddot{r}_2 \\
	\ddot{\theta}_1 \\
	\ddot{\theta}_2
	\end{bmatrix} &= \begin{bmatrix}
	r_1\dot{\theta}_1^2-g\sin{\theta_1} + \dfrac{m_2}{m_1+m_2}\left(r_2\dot{\theta}_2^2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_2\sin{\Delta}\right)  + \dfrac{Q_{r_1}-k_1(r_1-l_1)}{m_1+m_2}\\
	r_2\dot{\theta}_2^2-g\sin{\theta_2} + r_1\dot{\theta}_1^2\cos{\Delta} - 2\dot{r}_1\dot{\theta}_1\sin{\Delta} + \dfrac{Q_{r_2}-k_2(r_2-l_2)}{m_2} \\
	-\dfrac{2\dot{r}_1\dot{\theta}_1}{r_1} - \dfrac{g\cos{\theta_1}}{r_1} - \dfrac{m_2}{(m_1+m_2)r_1}\left[2\dot{r}_2\dot{\theta}_2\cos{\Delta} -r_2\dot{\theta}_2^2\sin{\Delta}\right] + \dfrac{Q_{\theta_1}}{(m_1+m_2)r_1^2} \\
	-\dfrac{2\dot{r}_2\dot{\theta}_2}{r_2}- \dfrac{g\cos{\theta_2}}{r_2} - \dfrac{2\dot{r}_1\dot{\theta}_1\cos{\Delta}}{r_2} - \dfrac{r_1\dot{\theta}_1^2\sin{\Delta}}{r_2} + \dfrac{Q_{\theta_2}}{m_2r_2^2}
	\end{bmatrix}.
\end{align*}


So the solution is
    \begin{align*}
        \begin{bmatrix}
\ddot{r}_1 \\
\ddot{r}_2 \\
\ddot{\theta}_1 \\
\ddot{\theta}_2
        \end{bmatrix} &= \begin{bmatrix}
        1 & \dfrac{m_2\cos{\Delta}}{m_1+m_2} & 0 & -\dfrac{m_2r_2\sin{\Delta}}{m_1+m_2} \\
        \cos{\Delta} & 1 & r_1\sin{\Delta} & 0 \\
        0 & \dfrac{m_2\sin{\Delta}}{(m_1+m_2)r_1} & 1 & \dfrac{m_2r_2\cos{\Delta}}{(m_1+m_2)r_1} \\
        -\dfrac{\sin{\Delta}}{r_2} & 0 & \dfrac{r_1\cos{\Delta}}{r_2} & 1
        \end{bmatrix}^{-1} \begin{bmatrix}
        r_1\dot{\theta}_1^2-g\sin{\theta_1} + \dfrac{m_2}{m_1+m_2}\left(r_2\dot{\theta}_2^2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_2\sin{\Delta}\right)  + \dfrac{Q_{r_1}-k_1(r_1-l_1)}{m_1+m_2}\\
        r_2\dot{\theta}_2^2-g\sin{\theta_2} + r_1\dot{\theta}_1^2\cos{\Delta} - 2\dot{r}_1\dot{\theta}_1\sin{\Delta} + \dfrac{Q_{r_2}-k_2(r_2-l_2)}{m_2} \\
        -\dfrac{2\dot{r}_1\dot{\theta}_1}{r_1} - \dfrac{g\cos{\theta_1}}{r_1} - \dfrac{m_2}{(m_1+m_2)r_1}\left[2\dot{r}_2\dot{\theta}_2\cos{\Delta} -r_2\dot{\theta}_2^2\sin{\Delta}\right] + \dfrac{Q_{\theta_1}}{(m_1+m_2)r_1^2} \\
        -\dfrac{2\dot{r}_2\dot{\theta}_2}{r_2}- \dfrac{g\cos{\theta_2}}{r_2} - \dfrac{2\dot{r}_1\dot{\theta}_1\cos{\Delta}}{r_2} - \dfrac{r_1\dot{\theta}_1^2\sin{\Delta}}{r_2} + \dfrac{Q_{\theta_2}}{m_2r_2^2}
        \end{bmatrix}.
    \end{align*}

Ordinarily, I would use a symplectic method to integrate a classical mechanical system, as it minimizes cumulative errors in the Hamiltonian over long-term integration. That being said, symplectic methods typically are either implicit (which significantly increases their computational complexity) or require a separable Hamiltonian. The Hamiltonian of this system would not separable as the kinetic energy cannot be written entirely in terms of momenta. Additionally, rewriting the kinetic energy in terms of the momenta would be tedious and prone to error.
