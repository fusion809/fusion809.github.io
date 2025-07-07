@def title = "Derivation of the equations of motion of the double pendulum"
@def mintoclevel=1
@def maxtoclevel=2

This webpage derives the equations of motion for the double pendulum. To see their solution, obtained via the  [Runge-Kutta-Fehlberg](https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method) fourth-order method with fifth-order error checking (RKF45), go to [this page](/doublePendulum/).

The problem of the double pendulum has been widely analysed in the mechanics literature. Its most simplified version involves two pendulums, with the second attached to the first at the bob, the rods being massless and there being no friction and air resistance. Even this simplified version is subject to chaos and has equations of motion far more complicated than that of the simple pendulum. In this webpage, we will be deriving the equations of motion with dissipative forces and rods with mass. 

~~~
<figure>
    <img src="/doublePendulum/Double pendulum.svg" width="500px"></img>
    <figcaption style="font-weight: bold; font-size: 18px;">Figure 1: Diagram of the double pendulum.</b></figcaption>
</figure>
~~~

Figure 1 depicts the double pendulum system being analysed. Most physicists treat $\theta_1$ and $\theta_2$ as being relative to the negative y-axis instead of the positive x-axis. I, as primarily a mathematician, prefer the positive x-axis as the point of reference for angles as I find it more intuitive. 

\tableofcontents
	
# Euler-Lagrange equations with generalized dissipation force

To obtain our equations of motion, we must apply the Euler-Lagrange equations with generalized dissipation force ($Q_i$ which is on the right-hand side of the following equation)
	
\begin{align}
	\dfrac{d}{dt}\dfrac{\partial \mathcal{L}}{\partial \dot{q_i}} - \dfrac{\partial \mathcal{L}}{\partial q_i} &=  Q_i. \label{ELD}
\end{align}
	
Where
* $\mathcal{L}$ is the Lagrangian &mdash; the difference between the kinetic energy and potential energy of the system. 
* $(q_i)=(\theta_1,\theta_2)$ are the generalized coordinates of the system.
* $Q_i=\displaystyle \sum_j \vec{F}_{D,j} \cdot \hat{e}_{j,i}$ is the generalized dissipation force.
* $j=1b,1r,2b,2r$ refers to the moving components of our double pendulum system, namely the two bobs and two rods.
* $\vec{F}_{D,j}$ is the dissipation forces acting on the $j$th component of the system.
* $\vec{r}_j$ is the position vector for the $j$th component of the system.
* $\hat{e}_{j,i}=\dfrac{\partial \vec{r}_j}{\partial q_i}$ is the generalized basis vector.
	
# Position, velocity and generalized basis vector
## Bob 1
The first pendulum bob would have the coordinates, velocity and generalized basis vector
\begin{align*}
		x_{1b} &= l_1 \cos{\theta_1} &\implies \dot{x}_{1b} &= -l_1\dot{\theta}_1 \sin{\theta_1}\\
		y_{1b} &= l_1 \sin{\theta_1} &\implies \dot{y}_{1b} &= l_1 \dot{\theta}_1 \cos{\theta_1} \\
		\vec{r}_{1b} &= \begin{bmatrix}
			x_{1b}\\
			y_{1b}
		\end{bmatrix}  &\therefore \vec{v}_{1b} &= \dfrac{d \vec{r}_{1b}}{dt}\\
		&= l_1\begin{bmatrix}
				-\sin{\theta_1} \\
				\cos{\theta_1}
			\end{bmatrix} & &=l_1\dot{\theta}_1 \begin{bmatrix}
			-\sin{\theta_1} \\
			\cos{\theta_1}
		\end{bmatrix} \\
		\vec{e}_{1b,\theta_1} &= \dfrac{\partial \vec{r}_{1b}}{\partial \theta_1}& v_{1b}^2 &= \dot{x}_{1b}^2 + \dot{y}_{1b}^2 \\
		&= l_1 \begin{bmatrix}
			-\sin{\theta_1}\\
			\cos{\theta_1}
		\end{bmatrix} & &= \left(-l_1\dot{\theta}_1 \sin{\theta_1}\right)^2 + \left(l_1 \dot{\theta}_1 \cos{\theta_1}\right)^2 \\
		\vec{e}_{1b,\theta_2} &= \dfrac{\partial \vec{r}_{1b}}{\partial \theta_2} &&= l_1^2 \dot{\theta}_1^2\\
		&= \begin{bmatrix}
			0\\
			0
		\end{bmatrix}.
\end{align*}
	
## Rod 1
The first pendulum rod would have the coordinates, velocity and generalized basis vector below. It is important to note that we are analysing its effect on the motion of the pendulum based on the centre of mass approach. 
	\begin{align*}
		x_{1r} &= \dfrac{l_1\cos{\theta_1}}{2} &\implies \dot{x}_{1r} &= -\dfrac{l_1\dot{\theta}_1\sin{\theta_1}}{2}\\
		y_{1r} &= \dfrac{l_1\sin{\theta_1}}{2} &\implies \dot{y}_{1r} &= \dfrac{l_1\dot{\theta}_1\cos{\theta_1}}{2}\\
		\vec{r}_{1r} &= \begin{bmatrix}
			x_{1r} \\
			y_{1r}
		\end{bmatrix} &\therefore \vec{v}_{1r} &= \dfrac{d\vec{r}_{1r}}{dt} \\
		&= \dfrac{l_1}{2}\begin{bmatrix}
			\cos{\theta_1}\\
			\sin{\theta_1}
		\end{bmatrix} & &= \dfrac{l_1\dot{\theta}_1}{2} \begin{bmatrix}
			-\sin{\theta_1}\\
			\cos{\theta_1}
		\end{bmatrix} \\
		\vec{e}_{1r,\theta_1} &= \dfrac{\partial \vec{r}_{1r}}{\partial \theta_1} & v_{1r}^2 &= \dot{x}_{1r}^2 + \dot{y}_{1r}^2 \\
		&= \dfrac{l_1}{2} \begin{bmatrix}
			-\sin{\theta_1}\\
			\cos{\theta_1}
		\end{bmatrix} & &= \left(-\dfrac{l_1\dot{\theta}_1 \sin{\theta_1}}{2}\right)^2 + \left(\dfrac{l_1 \dot{\theta}_1 \cos{\theta_1}}{2}\right)^2 \\
	\vec{e}_{1r,\theta_2} &= \dfrac{\partial \vec{r}_{1r}}{\partial \theta_2} &&= \dfrac{l_1^2 \dot{\theta}_1^2}{4} \\
	&= \begin{bmatrix}
		0 \\
		0
	\end{bmatrix}.
\end{align*}
	
## Bob 2
	\begin{align*}
		x_{2b} &= x_{1b} + l_2 \cos{\theta_2} &\dot{x}_{2b} &= -l_1\dot{\theta}_1 \sin{\theta_1}-l_2\dot{\theta}_2 \sin{\theta_2}\\
		&= l_1 \cos{\theta_1} + l_2 \cos{\theta_2} & \dot{y}_{2b} &= l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2} \\ 
		y_{2b} &= y_{1b} + l_2 \sin{\theta_2} & \therefore \vec{v}_{2b} &= \begin{bmatrix}
			-l_1\dot{\theta}_1 \sin{\theta_1}-l_2\dot{\theta}_2 \sin{\theta_2}\\
			l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2}
		\end{bmatrix}\\
		&= l_1 \sin{\theta_1} + l_2 \sin{\theta_2} & \therefore v_{2b}^2 &= \dot{x}_{2b}^2 + \dot{y}_{2b}^2 \\
		\vec{r}_{2b}&= \begin{bmatrix}
			x_{2b} \\
			y_{2b}
		\end{bmatrix}& &= \left(-l_1\dot{\theta}_1 \sin{\theta_1}-l_2\dot{\theta}_2 \sin{\theta_2}\right)^2 + \left(l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2}\right)^2 \\
		&= \begin{bmatrix}
			l_1 \cos{\theta_1} + l_2 \cos{\theta_2}\\
			l_1 \sin{\theta_1} + l_2 \sin{\theta_2}
		\end{bmatrix} & &= l_1^2 \dot{\theta}_1^2 \sin^2{\theta_1} + l_2^2 \dot{\theta}_2^2 \sin^2{\theta_2} + 2l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{\theta_1}\sin{\theta_2} + l_1^2 \dot{\theta}_1^2 \cos^2{\theta_1} + l_2^2 \dot{\theta}_2^2 \cos^2{\theta_2} + \\
		\vec{e}_{2b,\theta_1} &= \dfrac{\partial \vec{r}_{2b}}{\partial \theta_1} && 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\theta_1}\cos{\theta_2} \\
		&= l_1\begin{bmatrix}
			-\sin{\theta_1} \\
			\cos{\theta_1}
		\end{bmatrix}& &= l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\left(\theta_1-\theta_2\right)}\\
		\vec{e}_{2b,\theta_2} &= \dfrac{\partial \vec{r}_{2b}}{\partial \theta_2} \\
		&= l_2 \begin{bmatrix}
			-\sin{\theta_2} \\
			\cos{\theta_2}
		\end{bmatrix}.
	\end{align*}
	
## Rod 2
	\begin{align*}
		x_{2r} &= x_{1b} + \dfrac{l_2 \cos{\theta_2}}{2} &\dot{x}_{2r} &= -l_1\dot{\theta}_1 \sin{\theta_1}-\dfrac{l_2\dot{\theta}_2 \sin{\theta_2}}{2}\\
		&= l_1 \cos{\theta_1} + \dfrac{l_2 \cos{\theta_2}}{2} & \dot{y}_{2r} &= l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2} \\ 
		y_{2r} &= y_{1b} + \dfrac{l_2 \sin{\theta_2}}{2} & \therefore \vec{v}_{2r} &= \begin{bmatrix}
			-l_1\dot{\theta}_1 \sin{\theta_1}-\dfrac{l_2\dot{\theta}_2 \sin{\theta_2}}{2}\\
			l_1\dot{\theta}_1 \cos{\theta_1}+\dfrac{l_2\dot{\theta}_2 \cos{\theta_2}}{2}
		\end{bmatrix}\\
		&= l_1 \sin{\theta_1} + \dfrac{l_2 \sin{\theta_2}}{2} & \therefore v_{2r}^2 &= \dot{x}_{2r}^2 + \dot{y}_{2r}^2 \\
		\vec{r}_{2r}&= \begin{bmatrix}
			x_{2r} \\
			y_{2r}
		\end{bmatrix}& &= \left(-l_1\dot{\theta}_1 \sin{\theta_1}-\dfrac{l_2\dot{\theta}_2 \sin{\theta_2}}{2}\right)^2 + \left(l_1\dot{\theta}_1 \cos{\theta_1}+\dfrac{l_2\dot{\theta}_2 \cos{\theta_2}}{2}\right)^2 \\
		&= \begin{bmatrix}
			l_1 \cos{\theta_1} + \dfrac{l_2 \cos{\theta_2}}{2}\\
			l_1 \sin{\theta_1} + \dfrac{l_2 \sin{\theta_2}}{2}
		\end{bmatrix} & &= l_1^2 \dot{\theta}_1^2 \sin^2{\theta_1} + \dfrac{l_2^2 \dot{\theta}_2^2 \sin^2{\theta_2}}{4} + l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{\theta_1}\sin{\theta_2} + l_1^2 \dot{\theta}_1^2 \cos^2{\theta_1} + \dfrac{l_2^2 \dot{\theta}_2^2 \cos^2{\theta_2}}{4} + \\
		\vec{e}_{2b,\theta_1} &= \dfrac{\partial \vec{r}_{2b}}{\partial \theta_1} && l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\theta_1}\cos{\theta_2} \\
		&= l_1\begin{bmatrix}
			-\sin{\theta_1} \\
			\cos{\theta_1}
		\end{bmatrix}& &= l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\left(\theta_1-\theta_2\right)}\\
		\vec{e}_{2b,\theta_2} &= \dfrac{\partial \vec{r}_{2b}}{\partial \theta_2} \\
		&= \dfrac{l_2}{2} \begin{bmatrix}
			-\sin{\theta_2} \\
			\cos{\theta_2}
		\end{bmatrix}.
	\end{align*}
	 
# Kinetic energy
	
## Bob 1
First, we will calculate the kinetic energy of the first pendulum bob. 
	\begin{align*}
		T_{1b} &= \dfrac{m_{1b}v_{1b}^2}{2}.
	\end{align*}
	
	Substituting in what we previously found $v_{1b}^2$ to be, we obtain
	
	\begin{align*}
		T_{1b} &= \dfrac{m_{1b} l_1^2 \dot{\theta}_1^2}{2}.
	\end{align*}
## Rod 1
	First we must determine the kinetic energy of the first rod, it is given by
	\begin{align*}
		T_{1r} &= \dfrac{1}{2} I_{\mathrm{cm}} \omega_{\mathrm{cm}}^2.
	\end{align*}
	Where $I_{\mathrm{cm}}$ is the moment of inertia around the centre of mass and $\omega_{\mathrm{cm}}$ is the angular velocity around the centre of mass. As we have a rod, $I_{\mathrm{cm}} = \dfrac{m_{1r}l_1^2}{12}$ and $\omega_{\mathrm{cm}} = \dot{\theta}_1$. This yields the following kinetic energy
	\begin{align*}
		T_{1r} &= \dfrac{1}{2}\dfrac{m_{1r} l_1^2}{12}\dot{\theta}_1^2 \\
		&= \dfrac{m_{1r} l_1^2\dot{\theta}_1^2}{24}.
	\end{align*}
	
## Bob 2
	The kinetic energy of the second bob is
	
	\begin{align*}
		T_{2b} &= \dfrac{m_{2b} v_{2b}^2}{2} \\
		&= \dfrac{m_{2b}}{2} \left(l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\left(\theta_1-\theta_2\right)}\right).
	\end{align*}
	
## Rod 2
	The second rod should be the same except as rod 1 but with its own parameters
	\begin{align*}
		T_{2r} &= \dfrac{m_{2r} l_2^2\dot{\theta}_2^2}{24}.
	\end{align*}
	
## Total
	Hence the total kinetic energy is
	
	\begin{align*}
		T &= \sum_j T_j \\
		&= T_{1r} + T_{1b} + T_{2r} + T_{2b} \\
		&= \dfrac{m_{1r}l_1^2\dot{\theta}_1^2}{24} + \dfrac{m_{1b} l_1^2 \dot{\theta}_1^2}{2} + \dfrac{m_{2r} l_2^2\dot{\theta}_2^2}{24} + \dfrac{m_{2b}}{2} \left(l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\left(\theta_1-\theta_2\right)}\right) \\
		&= \dfrac{m_{1r} l_1^2 \dot{\theta}_1^2 + m_{2r}l_2^2 \dot{\theta}_2^2}{24} + \dfrac{m_{1b}+m_{2b}}{2}l_1^2 \dot{\theta}_1^2 + \dfrac{m_{2b}}{2} \left(l_2\dot{\theta}_2^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}\right)
	\end{align*}
	\begin{align*}
		T &= \dfrac{l_1^2 \dot{\theta}_1^2}{2} \left(\dfrac{m_{1r}}{12} + m_{1b} + m_{2b}\right) + \dfrac{l_2^2 \dot{\theta}_2^2}{2} \left(\dfrac{m_{2r}}{12} + m_{2b}\right) + m_{2b}l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}.
	\end{align*}
	
# Potential energy
Next, we must calculate the potential energy. In this problem, there is only one source of potential energy---gravity. This means that the potential energy of each component of the pendulum will be $V_j = m_jgy_j$, where $y_j$ is the $y$-coordinate of the centre of mass of the component and $m_j$ is the mass of that component. 
	
## Bob 1
	First, we have the gravitational potential energy of the first pendulum bob. 
	\begin{align*}
		V_{1b} &= m_{1b} gy_{1b} \\
		&= m_{1b}gl_1 \sin{\theta_1}.
	\end{align*}
	
## Rod 1
	Next, we will calculate the gravitational potential energy of the first rod. For it, we will use the centre of mass of the pendulum for the $y$-coordinate we use to calculate the gravitational potential energy. We will assume it has uniform mass, so its midpoint will be where its centre of mass is. Therefore $y_{1r} = \dfrac{1}{2}y_{1b}$. 
	\begin{align*}
		V_{1r} &= m_{1r} gy_{1r} \\
		&= \dfrac{m_{1r}gl_1 \sin{\theta_1}}{2}.
	\end{align*}
	
## Bob 2
	The potential energy of the second pendulum bob is
	\begin{align*}
		V_{2b} &= m_{2b} gy_{2b} \\
		&= m_{2b} g \left(l_1 \sin{\theta_1} + l_2 \sin{\theta_2}\right).
	\end{align*}
	
## Rod 2
	The potential energy of the second pendulum rod is
	\begin{align*}
		V_{2r} &= m_{2r} gy_{2r} \\
		 &= m_{2r}g \left(l_1 \sin{\theta_1} + \dfrac{l_2\sin{\theta_2}}{2}\right).
	\end{align*}
	
## Total
	Therefore, the total potential energy is
	\begin{align*}
		V &= \sum_j V_j \\
		&= V_{1r} + V_{2r} + V_{1b} + V_{2b} \\
		&= \dfrac{m_{1r}gl_1 \sin{\theta_1}}{2} + m_{2r}g \left(l_1 \sin{\theta_1} + \dfrac{l_2\sin{\theta_2}}{2}\right) + m_{1b}gl_1 \sin{\theta_1} + m_{2b} g \left(l_1 \sin{\theta_1} + l_2 \sin{\theta_2}\right) \\
		&= gl_1 \sin{\theta_1}\left(\dfrac{m_{1r}}{2} + m_{2r} + m_{1b} + m_{2b}\right) + gl_2 \sin{\theta_2}\left(\dfrac{m_{2r}}{2} + m_{2b}\right).
	\end{align*}
	
# Lagrangian
Hence the Lagrangian is
	
	\begin{align*}
		\mathcal{L} &= T - V \\
		&= \dfrac{l_1^2 \dot{\theta}_1^2}{2} \left(\dfrac{m_{1r}}{12} + m_{1b} + m_{2b}\right) + \dfrac{l_2^2 \dot{\theta}_2^2}{2} \left(\dfrac{m_{2r}}{12} + m_{2b}\right) + m_{2b}l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)} - gl_1 \sin{\theta_1}\left(\dfrac{m_{1r}}{2} \right.\\
		&\left.+ m_{2r} + m_{1b} + m_{2b}\right) - gl_2 \sin{\theta_2}\left(\dfrac{m_{2r}}{2} + m_{2b}\right).
	\end{align*}
	
# Applying the Euler-Lagrange equations with generalized dissipative force
## $\theta_1$
### Left-hand side
	Calculating the left-hand side terms of Equation \eqref{ELD} with regards to $\theta_1$. First we calculate the generalized momenta canonical to $\theta_1$
	
	\begin{align*}
		p_{\theta_1} &= \dfrac{\partial \mathcal{L}}{\partial \dot{\theta}_1} \\
		&= \dfrac{m_{1r}l_1^2 \dot{\theta}_1}{12} + (m_{1b}+m_{2b}) l_1^2 \dot{\theta}_1 + m_{2b}l_1 l_2 \dot{\theta}_2 \cos{\left(\theta_1-\theta_2\right)}.
	\end{align*}
	
	Hence its time derivative is
	\begin{align*}
		\dot{p}_{\theta_1} &= \dfrac{m_{1r} l_1^2 \ddot{\theta}_1}{12} + (m_{1b}+m_{2b})l_1^2 \ddot{\theta}_1 + m_{2b}l_1 l_2 \ddot{\theta}_2\cos{(\theta_1-\theta_2)} - m_{2b}l_1 l_2 \dot{\theta}_2\left(\dot{\theta}_1 - \dot{\theta}_2\right)\sin{(\theta_1-\theta_2)}.
	\end{align*}
	
	As for the other term in Equation \eqref{ELD}
	\begin{align*}		
		F_{\theta_1} &= \dfrac{\partial \mathcal{L}}{\partial \theta_1} \\
		&= -m_{2b}l_1l_2\dot{\theta}_1\dot{\theta}_2 \sin{(\theta_1-\theta_2)} - \dfrac{m_{1r}gl_1 \cos{\theta_1}}{2} -m_{2r}gl_1 \cos{\theta_1} -m_{1b}gl_1 \cos{\theta_1} -m_{2b}gl_1 \cos{\theta_1} \\
		&= -m_{2b}l_1l_2\dot{\theta}_1\dot{\theta}_2 \sin{(\theta_1-\theta_2)} -gl_1\cos{\theta_1}\left(\dfrac{m_{1r}}{2}+m_{2r}+m_{1b} + m_{2b}\right).
	\end{align*}
	
	Hence the LHS of Equation \eqref{ELD} for $\theta_1$ is
	\begin{align*}
		\dot{p}_{\theta_1} - F_{\theta_1} &= \dfrac{m_{1r} l_1^2 \ddot{\theta}_1}{12} + (m_{1b}+m_{2b})l_1^2 \ddot{\theta}_1 + m_{2b}l_1 l_2 \ddot{\theta}_2\cos{(\theta_1-\theta_2)} - m_{2b}l_1 l_2 \dot{\theta}_2\left(\dot{\theta}_1 - \dot{\theta}_2\right)\sin{(\theta_1-\theta_2)} + m_{2b}l_1l_2\dot{\theta}_1\dot{\theta}_2\sin{(\theta_1-\theta_2)} +gl_1\cos{\theta_1}\left(\dfrac{m_{1r}}{2}+m_{2r}+m_{1b} + m_{2b}\right) \\
		&= \dfrac{m_{1r} l_1^2 \ddot{\theta}_1}{12} + (m_{1b}+m_{2b})l_1^2 \ddot{\theta}_1 + m_{2b}l_1 l_2 \ddot{\theta}_2\cos{(\theta_1-\theta_2)} + m_{2b}l_1 l_2 \dot{\theta}_2^2\sin{(\theta_1-\theta_2)} +gl_1\cos{\theta_1}\left(\dfrac{m_{1r}}{2}+m_{2r}+m_{1b} + m_{2b}\right) \\
		&= l_1^2 \ddot{\theta}_1 \left(\dfrac{m_{1r}}{12} + m_{1b} + m_{2b}\right) + m_{2b}l_1 l_2 \ddot{\theta}_2\cos{(\theta_1-\theta_2)} + m_{2b}l_1 l_2 \dot{\theta}_2^2\sin{(\theta_1-\theta_2)} +gl_1\cos{\theta_1}\left(\dfrac{m_{1r}}{2}+m_{2r}+m_{1b} + m_{2b}\right)\\
		&= l_1^2 \ddot{\theta}_1 \left(\dfrac{m_{1r}}{12} + m_{1b} + m_{2b}\right) + m_{2b}l_1 l_2\left( \ddot{\theta}_2\cos{(\theta_1-\theta_2)} + \dot{\theta}_2^2\sin{(\theta_1-\theta_2)}\right) +gl_1\cos{\theta_1}\left(\dfrac{m_{1r}}{2}+m_{2r}+m_{1b} + m_{2b}\right).
	\end{align*}
	
	Hence $\ddot{\theta}_1$ is
	\begin{align*}
		\ddot{\theta}_1 &= \dfrac{1}{l_1^2\left(\dfrac{m_{1r}}{12} + m_{1b} + m_{2b}\right)}\left[Q_{\theta_1} - m_{2b}l_1 l_2\left( \ddot{\theta}_2\cos{(\theta_1-\theta_2)} + \dot{\theta}_2^2\sin{(\theta_1-\theta_2)}\right) -gl_1\cos{\theta_1}\left(\dfrac{m_{1r}}{2}+m_{2r}+m_{1b} + m_{2b}\right)\right].
	\end{align*}
	
### Right-hand side
	As for the generalized dissipation force
	
	\begin{align*}
		Q_{\theta_1} &= \vec{F}_{D, 1r} \cdot \dfrac{\partial \vec{r}_{1r}}{\partial \theta_1} + \vec{F}_{D, 1b} \cdot \dfrac{\partial \vec{r}_{1b}}{\partial \theta_1}+ \vec{F}_{D, 2r} \cdot \dfrac{\partial \vec{r}_{rod, 2}}{\partial \theta_1} + \vec{F}_{D, 2b} \cdot \dfrac{\partial \vec{r}_{2b}}{\partial \theta_1}
	\end{align*}
	
<b>Bob 1</b><br/>
	The dissipation force applied to the first bob is
	\begin{align*}
		\vec{F}_{D,1b} &= -b_{1b}\vec{v}_{1b} - c_{1b}|\vec{v}_{1b}|\vec{v}_{1b}.
	\end{align*}
	
	Hence the generalized dissipation force, canonical to $\theta_1$, applied to the first bob is
	\begin{align*}
		\vec{F}_{D, 1b} \cdot \dfrac{\partial \vec{r}_{1b}}{\partial \theta_1}&= (-b_{1b} - c_{1b}l_1\dot{\theta}_1)\begin{bmatrix}
			-l_1\dot{\theta}_1 \sin{\theta_1} \\
			l_1\dot{\theta}_1 \cos{\theta_1}
			\end{bmatrix} \cdot \begin{bmatrix}
			-l_1 \sin{\theta_1} \\
			l_1 \cos{\theta_1}
		\end{bmatrix} \\
		&= -(b_{1b}+c_{1b}l_1\dot{\theta}_1) (l_1^2 \dot{\theta}_1 \sin^2{\theta_1}+l_1^2 \dot{\theta}_1 \cos^2{\theta_1}) \\
		&= -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1^2 \dot{\theta}_1.
	\end{align*}
	
<b>Rod 1</b><br/>
	
	Next we will calculate the dissipation force on the first rod. We will use a centre of mass approximation (as otherwise we would likely have to integrate over the rod, which would drastically complicate the calculation)
	
	\begin{align*}
		\vec{F}_{D,1r} &= -\left(b_{1r} + c_{1r}|\vec{v}_{1r}|\right)\vec{v}_{1r} \\
		\vec{F}_{D,1r} &= -\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right)\dfrac{l_1 \dot{\theta}_1}{2} \begin{bmatrix}
			-\sin{\theta_1} \\
			\cos{\theta_1}
		\end{bmatrix}.
	\end{align*}
	
	Hence the generalized dissipation force, canonical to $\theta_1$, on the first rod is
	\begin{align*}
		\vec{F}_{D,1r} \cdot \dfrac{\partial \vec{r}_{1r}}{\partial \theta_1} &= -\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right)\dfrac{l_1 \dot{\theta}_1}{2} \begin{bmatrix}
			-\sin{\theta_1} \\
			\cos{\theta_1}
		\end{bmatrix} \cdot \dfrac{l_1}{2}\begin{bmatrix}
			-\sin{\theta_1} \\
			\cos{\theta_1}
		\end{bmatrix} \\
		&= -\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1^2 \dot{\theta}_1}{4}.
	\end{align*}

<b>Bob 2</b><br/>
	The dissipation force on the second bob would be
	
	\begin{align*}
		\vec{F}_{D,2b} &= -b_{2b}\vec{v}_{2b} - c_{2b}|\vec{v}_{2b}|\vec{v}_{2b} \\
		&= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right) \begin{bmatrix}
			-l_1\dot{\theta}_1 \sin{\theta_1}-l_2\dot{\theta}_2 \sin{\theta_2} \\
			l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2}
		\end{bmatrix}.
	\end{align*}
	
	Hence the generalized dissipation force, canonical to $\theta_1$, would be
	\begin{align*}
		\vec{F}_{D, 2b} \cdot \dfrac{\partial \vec{r}_{2b}}{\partial \theta_1} &= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right) \begin{bmatrix}
			-l_1\dot{\theta}_1 \sin{\theta_1}-l_2\dot{\theta}_2 \sin{\theta_2} \\
			l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2}
		\end{bmatrix} \cdot \begin{bmatrix}
		-l_1 \sin{\theta_1} \\
		l_1 \cos{\theta_1}
		\end{bmatrix}\\
		&= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 \sin^2{\theta_1}+l_1l_2 \dot{\theta}_2 \sin{\theta_1}\sin{\theta_2} +l_1^2 \dot{\theta_1}\cos^2{\theta_1}+l_1l_2\dot{\theta}_2 \cos{\theta_1}\cos{\theta_2}\right) \\
		&= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 + l_1l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}\right).
	\end{align*}
	
<b>Rod 2</b><br/>
	The friction force applied to the second rod would be
	\begin{align*}
		\vec{F}_{D,2r} &= -\left(b_{2r} + c_{2r}|\vec{v}_{2r}|\right)\vec{v}_{2r} \\
		\vec{F}_{D,2r} &= -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\begin{bmatrix}
			-l_1 \dot{\theta}_1\sin{\theta_1} - \dfrac{l_2\dot{\theta}_2\sin{\theta_2}}{2}\\
			l_1 \dot{\theta}_1 \cos{\theta_1} + \dfrac{l_2\dot{\theta}_2\cos{\theta_2}}{2}
		\end{bmatrix}.
	\end{align*}
	
	Therefore the generalized dissipation force, canonical to $\theta_1$, for the second rod would be 
	\begin{align*}	
		\vec{F}_{D,2r}\cdot \dfrac{\partial \vec{r}_{2r}}{\partial \theta_1} &= -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right) \left(l_1^2 \dot{\theta}_1 \sin^2{\theta_1} + \dfrac{l_1 l_2 \dot{\theta}_2 \sin{\theta_1}\sin{\theta_2}}{2} +l_1^2 \dot{\theta}_1\cos^2{\theta_1} + \dfrac{l_1l_2 \dot{\theta}_2 \cos{\theta_1}\cos{\theta_2}}{2}\right) \\
		&= -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1^2 + \dfrac{l_1 l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right).
	\end{align*}

<b>Final answer</b>
	Hence $Q_{\theta_1}$ is
	\begin{align*}
		Q_{\theta_1} &= -\left(b_{1b} + c_{1b} l_1 \dot{\theta}_1\right)l_1^2 \dot{\theta}_1 -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 + l_1l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}\right) -\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1^2 \dot{\theta}_1}{4} \\
		&-\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1^2 + \dfrac{l_1 l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right).
	\end{align*}
	
 ### Final equation for $\ddot{\theta}_1$
	The equation of motion for $\theta_1$ is therefore
	
	\begin{align*}
		\dot{p}_{\theta_1} - F_{\theta_1} &= Q_{\theta_1}
	\end{align*}
	
	\begin{align*}
		&\dfrac{m_{1r} l_1^2 \ddot{\theta}_1}{12} + (m_{1b}+m_{2b})l_1^2 \ddot{\theta}_1 + m_{2b}l_1 l_2 \ddot{\theta}_2\cos{(\theta_1-\theta_2)} - m_{2b}l_1 l_2 \dot{\theta}_2\left(\dot{\theta}_1 - \dot{\theta}_2\right)\sin{(\theta_1-\theta_2)} - \left(-m_{2b}l_1l_2\dot{\theta}_1\dot{\theta}_2 \sin{(\theta_1-\theta_2)} - \dfrac{m_{1r}gl_1 \cos{\theta_1}}{2} -m_{2r}gl_1 \cos{\theta_1} -m_{1b}gl_1 \cos{\theta_1} \right.\\
		&\left.-m_{2b}gl_1 \cos{\theta_1}\right) = -\left(b_{1b} + c_{1b} l_1 \dot{\theta}_1\right)l_1^2 \dot{\theta}_1 -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 + l_1l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}\right) -\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1^2 \dot{\theta}_1}{4} \\
		& -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 + \dfrac{l_1 l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right) \\
		&\dfrac{m_{1r} l_1^2 \ddot{\theta}_1}{12} + (m_{1b}+m_{2b})l_1^2 \ddot{\theta}_1 + m_{2b}l_1 l_2 \ddot{\theta}_2\cos{(\theta_1-\theta_2)} + m_{2b}l_1 l_2 \dot{\theta}_2^2\sin{(\theta_1-\theta_2)} + \dfrac{m_{1r}gl_1 \cos{\theta_1}}{2} +m_{2r}gl_1 \cos{\theta_1} +m_{1b}gl_1 \cos{\theta_1} +m_{2b}gl_1 \cos{\theta_1} = -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1^2 \dot{\theta}_1 \\
		&-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)(l_1^2 \dot{\theta}_1 + l_1l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}) -\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1^2 \dot{\theta}_1}{4} \\
		&-\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 + \dfrac{l_1 l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right) \\
	\end{align*}
	\begin{align*}
		& \left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)l_1^2 \ddot{\theta}_1 = - m_{2b}l_1 l_2 \left( \ddot{\theta}_2\cos{(\theta_1-\theta_2)} +\dot{\theta}_2^2\sin{(\theta_1-\theta_2)}\right) - gl_1 \cos{\theta_1}\left(\dfrac{m_{1r}}{2} +m_{2r} +m_{1b} + m_{2b}\right) -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1^2 \dot{\theta}_1 \\
		&-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 + l_1l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}\right) -\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1^2 \dot{\theta}_1}{4} \\
		& -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1^2 \dot{\theta}_1 + \dfrac{l_1 l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right) \\
	\end{align*}
	\begin{align}
		\ddot{\theta}_1 = &\dfrac{-1}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)l_1} \left[m_{2b}l_2 ( \ddot{\theta}_2\cos{(\theta_1-\theta_2)} +\dot{\theta}_2^2\sin{(\theta_1-\theta_2)}) 
		+ g \cos{\theta_1}\left(\dfrac{m_{1r}}{2} +m_{2r} +m_{1b} + m_{2b}\right) -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1 \dot{\theta}_1 \right.\\
		&\left.+\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)(l_1 \dot{\theta}_1 + l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}) +\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1 \dot{\theta}_1}{4} \right.\\
		&\left.+\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1 \dot{\theta}_1 + \dfrac{l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right)\right] \label{d2theta11}
	\end{align}
	
## $\theta_2$
Next, we will derive the equations for $\theta_2$
	
### Left-hand side
\begin{align*}
	p_{\theta_2} &= \dfrac{\partial \mathcal{L}}{\partial \dot{\theta}_2} \\
		&= \dfrac{m_{2r}l_2^2 \dot{\theta}_2}{12} + m_{2b}l_2 \left(l_2 \dot{\theta}_2 + l_1 \dot{\theta}_1\cos{(\theta_1-\theta_2)}\right)\\
		\dot{p}_{\theta_2} &= \dfrac{m_{2r}l_2^2 \ddot{\theta}_2}{12} + m_{2b}l_2 \left(l_2 \ddot{\theta}_2 + l_1 \ddot{\theta}_1\cos{(\theta_1-\theta_2)}-l_1\dot{\theta}_1(\dot{\theta}_1-\dot{\theta}_2)\sin{(\theta_1-\theta_2)}\right) \\
		&= l_2^2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}l_1l_2\left(\ddot{\theta}_1\cos{(\theta_1-\theta_2)}-\dot{\theta}_1(\dot{\theta}_1-\dot{\theta}_2)\sin{(\theta_1-\theta_2)}\right)\\
		F_{\theta_2} &= \dfrac{\partial \mathcal{L}}{\partial \theta_2} \\
		&= m_{2b} l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \sin{(\theta_1-\theta_2)}-\dfrac{m_{2r}gl_2 \cos{\theta_2}}{2} - m_{2b}gl_2\cos{\theta_2} \\
		&= m_{2b} l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \sin{(\theta_1-\theta_2)}-gl_2\cos{\theta_2}\left(\dfrac{m_{2r}}{2} + m_{2b}\right).
\end{align*}
The LHS of our equation will be
	
\begin{align*}
		\dot{p}_{\theta_2} - F_{\theta_2} &= l_2^2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}l_1l_2\left(\ddot{\theta}_1\cos{(\theta_1-\theta_2)}-\dot{\theta}_1(\dot{\theta}_1-\dot{\theta}_2)\sin{(\theta_1-\theta_2)}\right) - m_{2b} l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \sin{(\theta_1-\theta_2)}+gl_2\cos{\theta_2}\left(\dfrac{m_{2r}}{2} + m_{2b}\right) \\
		&= l_2^2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}l_1l_2\left(\ddot{\theta}_1\cos{(\theta_1-\theta_2)}-\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}\right)+gl_2\cos{\theta_2}\left(\dfrac{m_{2r}}{2} + m_{2b}\right).
\end{align*}
	
### Right-hand side
\begin{align*}
		Q_{\theta_2} &= \sum_{j} \vec{F}_{D,j} \cdot \dfrac{\partial \vec{r}_j}{\partial \theta_2} \\
		&= \vec{F}_{D,1r} \cdot \dfrac{\partial \vec{r}_{1r}}{\partial \theta_2} + \vec{F}_{D,1b} \cdot \dfrac{\partial \vec{r}_{1b}}{\partial \theta_2} + \vec{F}_{D,2r} \cdot \dfrac{\partial \vec{r}_{2r}}{\partial \theta_2} + \vec{F}_{D,2b} \cdot \dfrac{\partial \vec{r}_{2b}}{\partial \theta_2}.
\end{align*}
	
The terms corresponding to the first rod and first bob will be zero as their position and velocities are independent of $\theta_2$. The remaining terms are below.
	
<b>Rod 2</b><br/>
The dissipation force applied to the second rod is
	\begin{align*}
		\vec{F}_{D,2r} &= -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\begin{bmatrix}
			-l_1 \dot{\theta}_1\sin{\theta_1} - \dfrac{l_2\dot{\theta}_2\sin{\theta_2}}{2}\\
			l_1 \dot{\theta}_1 \cos{\theta_1} + \dfrac{l_2\dot{\theta}_2\cos{\theta_2}}{2}
		\end{bmatrix}.
	\end{align*}
	Hence the generalized dissipation force, canonical to $\theta_2$, applied to the second rod is
	\begin{align*}
		\vec{F}_{D,2r} \cdot \dfrac{\partial \vec{r}_{2r}}{\partial \theta_2} &= -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\begin{bmatrix}
			-l_1 \dot{\theta}_1\sin{\theta_1} - \dfrac{l_2\dot{\theta}_2\sin{\theta_2}}{2}\\
			l_1 \dot{\theta}_1 \cos{\theta_1} + \dfrac{l_2\dot{\theta}_2\cos{\theta_2}}{2}
		\end{bmatrix} \cdot \dfrac{l_2}{2} \begin{bmatrix}
			-\sin{\theta_2} \\
			\cos{\theta_2}
		\end{bmatrix} \\
		&= -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(\dfrac{l_1l_2 \dot{\theta}_1 \sin{\theta_1}\sin{\theta_2}}{2} + \dfrac{l_2^2 \dot{\theta}_2\sin^2{\theta_2}}{4}+\dfrac{l_1l_2\dot{\theta}_1\cos{\theta_1}\cos{\theta_2}}{2}+\dfrac{l_2^2\dot{\theta}_2\cos^2{\theta_2}}{4}\right) \\
		&= -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)(2l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2^2 \dot{\theta}_2).
	\end{align*}

<b>Bob 2</b><br/>
	The dissipation force applied to the second bob is
	\begin{align*}
		\vec{F}_{D,2b} &= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right) \begin{bmatrix}
			-l_1\dot{\theta}_1 \sin{\theta_1}-l_2\dot{\theta}_2 \sin{\theta_2} \\
			l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2}
		\end{bmatrix}.
	\end{align*}
	Therefore the generalized dissipation force, canonical to $\theta_2$, applied to the second bob is
	\begin{align*}
		\vec{F}_{D,2b} \cdot \dfrac{\partial \vec{r}_{2b}}{\partial \theta_2} &= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right) \begin{bmatrix}
			-l_1\dot{\theta}_1 \sin{\theta_1}-l_2\dot{\theta}_2 \sin{\theta_2} \\
			l_1\dot{\theta}_1 \cos{\theta_1}+l_2\dot{\theta}_2 \cos{\theta_2}
		\end{bmatrix} \cdot l_2\begin{bmatrix}
		-\sin{\theta_2} \\
		\cos{\theta_2}
		\end{bmatrix} \\
		&= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right) \left(l_1 l_2 \dot{\theta}_1 \sin{\theta_1}\sin{\theta_2} + l_2^2 \dot{\theta}_2 \sin^2{\theta_2} + l_1l_2 \dot{\theta}_1 \cos{\theta_1}\cos{\theta_2} + l_2^2\dot{\theta}_2 \cos^2{\theta_2}\right) \\
		&= -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2^2 \dot{\theta}_2\right).
	\end{align*}
	
### Right-hand side
\begin{align*}
		Q_{\theta_2} &= -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(2l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2^2 \dot{\theta}_2\right)\\
		& -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2^2 \dot{\theta}_2\right).
	\end{align*}
	
### Final equation
	\begin{align*}
		&l_2^2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}l_1l_2\left(\ddot{\theta}_1\cos{(\theta_1-\theta_2)}-\dot{\theta}_1(\dot{\theta}_1-\dot{\theta}_2)\sin{(\theta_1-\theta_2)}\right) - m_{2b} l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \sin{(\theta_1-\theta_2)} +\dfrac{m_{2r}gl_2 \cos{\theta_2}}{2}+ m_{2b}gl_2\cos{\theta_2}\\
		&= -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)(2l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2^2 \dot{\theta}_2) -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)}\right. \\
		&\left.+ l_2^2 \dot{\theta}_2\right)
	\end{align*}
	Moving everything on the left-hand side except the $\ddot{\theta}_1$ term to the right-hand side yields
	\begin{align*}
		&m_{2b}l_1l_2\cos{(\theta_1-\theta_2)}\ddot{\theta}_1 = -l_2^2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}l_2(l_1\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}-g\cos{\theta_2})
		-\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)(2l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} \\
		&+ l_2^2 \dot{\theta}_2) -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1l_2 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2^2 \dot{\theta}_2\right).
	\end{align*}
	
	Dividing both sides by $m_{2b}l_1l_2 \cos{(\theta_1-\theta_2)}$:
	
	\begin{align*}
		&\ddot{\theta}_1 = \dfrac{\sec{(\theta_1-\theta_2)}}{m_{2b}l_1} \left[-l_2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}(l_1\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}-g\cos{\theta_2})-\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(2l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2 \dot{\theta}_2\right)\right.\\ &\left.-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2 \dot{\theta}_2\right)\right].
	\end{align*}
	
## Rewriting $\ddot{\theta}_2$ to not include $\ddot{\theta}_1$
	Dividing both sides by $m_2l_1l_2\cos{(\theta_1-\theta_2)}$ and replacing $\ddot{\theta}_1$ on the LHS with the right-hand side of Equation \eqref{d2theta11} yields
	
	\begin{align*}
		&\dfrac{1}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)l_1} [-m_{2b}l_2 ( \ddot{\theta}_2\cos{(\theta_1-\theta_2)} +\dot{\theta}_2^2\sin{(\theta_1-\theta_2)}) 
		- g \cos{\theta_1}\left(\dfrac{m_{1r}}{2} +m_{2r} +m_{1b} + m_{2b}\right) -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1 \dot{\theta}_1 \\
		&-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)(l_1 \dot{\theta}_1 + l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)})-\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1 \dot{\theta}_1}{4} -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)(l_1 \dot{\theta}_1 \\
		&+ \dfrac{l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2})] =\dfrac{\sec{(\theta_1-\theta_2)}}{m_{2b}l_1} [-l_2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}(l_1\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}-g\cos{\theta_2}) -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\\
		&(2l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)}+ l_2 \dot{\theta}_2) -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)(l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2 \dot{\theta}_2)]\\
	\end{align*}
	
	Multiply both sides by $l_1m_{2b}\cos{(\theta_1-\theta_2)}$
	
	\begin{align*}
		&\dfrac{m_{2b}\cos{(\theta_1-\theta_2)}}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)} \left[-m_{2b}l_2 \left( \ddot{\theta}_2\cos{(\theta_1-\theta_2)} +\dot{\theta}_2^2\sin{(\theta_1-\theta_2)}\right) 
		- g \cos{\theta_1}\left(\dfrac{m_{1r}}{2} +m_{2r} +m_{1b} + m_{2b}\right) -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1 \dot{\theta}_1\right. \\
		&\left.-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)(l_1 \dot{\theta}_1 + l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)})-\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1 \dot{\theta}_1}{4} -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1 \dot{\theta}_1\right.\right. \\
		&\left.\left.+ \dfrac{l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right)\right] =-l_2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2 + m_{2b}(l_1\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}-g\cos{\theta_2}) -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\\
		&\left(2l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)}+ l_2 \dot{\theta}_2\right) -\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2 \dot{\theta}_2\right)\\
	\end{align*}
	
	Adding $l_2 \left(\dfrac{m_{2r}}{12} + m_{2b}\right)\ddot{\theta}_2$ to both sides
	
	\begin{align*}
		&\ddot{\theta}_2\left(\left(\dfrac{m_{2r}}{12} + m_{2b}\right)l_2 - \dfrac{m_{2b}^2l_2\cos^2{(\theta_1-\theta_2)}}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)}\right) + \dfrac{m_{2b}\cos{(\theta_1-\theta_2)}}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)}\left[-m_{2b}l_2\dot{\theta}_2^2\sin{(\theta_1-\theta_2)} - g \cos{\theta_1}\left(\dfrac{m_{1r}}{2} +m_{2r} +m_{1b} + m_{2b}\right) -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1 \dot{\theta}_1\right. \\
		&\left.-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1 \dot{\theta}_1 + l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}\right)-\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1 \dot{\theta}_1}{4} -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1 \dot{\theta}_1\right.\right. \\
		&\left.\left.+ \dfrac{l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right)\right] =   m_{2b}\left(l_1\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}-g\cos{\theta_2}\right) -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)(2l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)}+ l_2 \dot{\theta}_2)\\
		&-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2 \dot{\theta}_2\right)\\\\\\
		&\ddot{\theta}_2\left(\left(\dfrac{m_{2r}}{12} + m_{2b}\right)l_2 - \dfrac{m_{2b}^2l_2\cos^2{(\theta_1-\theta_2)}}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)}\right) = -\dfrac{m_{2b}\cos{(\theta_1-\theta_2)}}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)}\left[-m_{2b}l_2\dot{\theta}_2^2\sin{(\theta_1-\theta_2)} - g \cos{\theta_1}\left(\dfrac{m_{1r}}{2} +m_{2r} +m_{1b} + m_{2b}\right) -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1 \dot{\theta}_1 \right.\\
		&\left.-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1 \dot{\theta}_1 + l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}\right)-\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1 \dot{\theta}_1}{4} -\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1 \dot{\theta}_1 \right.\right.\\
		&\left.\left.+ \dfrac{l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right)\right] + m_{2b}\left(l_1\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}-g\cos{\theta_2}\right) -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(2l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)}+ l_2 \dot{\theta}_2\right)\\
		&-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)\left(l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2 \dot{\theta}_2\right)\\\\\\
	\end{align*}

	\begin{align}
		&\ddot{\theta}_2 = \dfrac{1}{\left(\dfrac{m_{2r}}{12} + m_{2b}\right)l_2 - \dfrac{m_{2b}^2l_2\cos^2{(\theta_1-\theta_2)}}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)}}\left[\dfrac{m_{2b}\cos{(\theta_1-\theta_2)}}{\left(\dfrac{m_{1r}}{12} + m_{1b}+m_{2b}\right)}\left[m_{2b}l_2\dot{\theta}_2^2\sin{(\theta_1-\theta_2)} + g \cos{\theta_1}\left(\dfrac{m_{1r}}{2} +m_{2r} +m_{1b} + m_{2b}\right) -(b_{1b} + c_{1b} l_1 \dot{\theta}_1)l_1 \dot{\theta}_1 \right.\right.\\
		&\quad\left.\left.+\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1 l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)(l_1 \dot{\theta}_1 + l_2 \dot{\theta}_2 \cos{(\theta_1-\theta_2)})+\left(b_{1r} + \dfrac{c_{1r}l_1 \dot{\theta}_1}{2}\right) \dfrac{l_1 \dot{\theta}_1}{4} +\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)\left(l_1 \dot{\theta}_1 \right.\right.\right. \\
		&\quad\left.\left.\left.+ \dfrac{l_2\dot{\theta}_2 \cos{\left(\theta_1 - \theta_2\right)}}{2}\right)\right] + m_{2b}(l_1\dot{\theta}_1^2\sin{(\theta_1-\theta_2)}-g\cos{\theta_2}) -\dfrac{1}{4}\left(b_{2r} + c_{2r}\sqrt{l_1^2 \dot{\theta}_1^2 + \dfrac{l_2^2 \dot{\theta}_2^2}{4} + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1 -\theta_2)}}\right)(2l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)}+ l_2 \dot{\theta}_2)\right.\\
		&\quad\left.-\left(b_{2b}+c_{2b}\sqrt{l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 +2l_1l_2\dot{\theta}_1 \dot{\theta}_2 \cos{(\theta_1-\theta_2)}}\right)(l_1 \dot{\theta}_1 \cos{(\theta_1-\theta_2)} + l_2 \dot{\theta}_2)\right].\label{d2theta2}
	\end{align}
	
	At this point, we could rewrite one of our $\ddot{\theta}_1$ equations with $\ddot{\theta}_2$ replaced with the right-hand side of Equation \eqref{d2theta2}. We will not do this as this will add even more potential for errors to creep in, and Equations \eqref{d2theta11} and \eqref{d2theta2} are already suitable for numerical integration. 
