@def hassim=false;
@def title="Modelling downhill descent with particular reference to Castle Hill"
@def tag=["Mathematics", "Physics"]
@def mintoclevel=1


~~~
<figure style="float: left; margin-right: 20px; margin-bottom: 15px; width: 600px;">
<img src="https://fusion809.github.io/images/Castle%20Hill/20261012_063443.jpg" width="100%"/>
<figcaption>View of Castle Hill, the principal hill modelled in this article, from Jezzine Barracks.</figcaption>
</figure>
~~~

~~~
<figure style="float: right; margin-right: 20px; margin-bottom: 15px; width: 400px;">
<img src="https://fusion809.github.io/images/Castle%20Hill/20261018_043344.jpg" width="100%"/>
<figcaption>My tricycle at the top of Castle Hill on 18 October 2026.</figcaption>
</figure>
~~~

This article derives the equations of motion of a tricycle going downhill after brake failure. This article definitely does not advocate sabotaging one's brakes, as modelling the descent clearly shows that potentially lethal and uncontrollable speeds are achieved rather quickly. Further, on Castle Hill, where many of these models are specifically applied to, the velocities modelled would be illegal given it has a 40 km/hr speed limit. For those unfamiliar, Castle Hill is a 286-metre granite monolith in Townsville, Australia. 

See [this notebook](https://github.com/fusion809/SageMath/blob/master/ChatGPT_playground.ipynb) showing some simulations, conducted with the help of ChatGPT, based on the equations of motion formulated in this article. 

\toc
	
# No rolling resistance and air resistance losses
In the event of there being no non-conservative forces like rolling resistance and air resistance acting on the trike, the total energy of the system is given by:
\begin{align*}
	E &= T + V
\end{align*}

Where $T$ is kinetic energy (in joules) and $V$ is potential energy (also in joules). In turn, they are defined as:
\begin{align}
	T=\dfrac{mv^2}{2} \label{T}
\end{align}

and

\begin{align}
	V=mgh(s). \label{Vc}
\end{align}

Where $v$ is velocity, $m$ is the mass of the tricycle (in kilograms), $g$ is gravitational acceleration (in metres per second squared ($\mathrm{m}\cdot \mathrm{s}^{-2}$)) and $h(s)$ is the height (in $m$) of the base of the tricycle above sea level. The $mgh$ value for $V$ assumes constant gravitational acceleration. 
	
Alternatively, if gravity is not constant:
\begin{align*}
	V &= \dfrac{GMm}{r_{\mathrm{base}}} - \dfrac{GMm}{r_{\mathrm{summit}}}\\
      &=GMm\left(r_{\mathrm{base}}^{-1}-r_{\mathrm{summit}}^{-1}\right).
\end{align*}
Where $G$ is the gravitational constant, $M$ is the mass of the Earth and $r_{\mathrm{summit}}$ and $r_{\mathrm{base}}$ are the distances of the summit and base of the hill/mountain from the centre of the Earth. 
	
If we are at the top of the hill and the tricycle is stationary, then $E=V$. So, when the trike moves downhill, $T$ will eventually equal $V$ at the bottom of the hill or mountain. Hence to determine the maximum of $v$ we merely need to solve:
	\begin{align*}
		\dfrac{mv^2}{2} &= mgh.
	\end{align*}
Multiplying both sides by $\dfrac{2}{m}$ yields:
	\begin{align*}    
		v^2 &= 2gh.
	\end{align*}
Whereas square rooting both sides yields:
	\begin{align}
		v &= \sqrt{2gh}. \label{vconst}
	\end{align}	
We deliberately ignored the negative solution for $v$ as we have defined $s$ to be increasing as we go downhill and hence velocity is positive in this direction. This also means that so long as motion is in the direction of downhill, which it always will be when it is moving solely under the force of gravity, velocity will always be positive. 
	
Or, if gravity cannot be assumed constant:
\begin{align*}
		\dfrac{mv^2}{2} &= GMm\left(r_{\mathrm{base}}^{-1}-r_{\mathrm{summit}}^{-1}\right)  \\
		v^2 &= 2GM \left(r_{\mathrm{base}}^{-1}-r_{\mathrm{summit}}^{-1}\right).
\end{align*}
Square rooting both sides yields:
\begin{align}
		v &= \sqrt{2GM\left(r_{\mathrm{base}}^{-1}-r_{\mathrm{summit}}^{-1}\right)}. \label{vnonconst}
\end{align}
Keeping in mind that the tallest point of Castle Hill is not reachable by road, so maximum altitude would be approximately somewhere in the vicinity of 265.9 metres above sea level --- at least based on what Google Fit says of my altitude when I journey there. According to Equations \eqref{vconst} and \eqref{vnonconst}, the maximum velocity of a trike travelling down Castle Hill would be approximately 260 kilometres per hour (kph). 
	
# With rolling resistance and air resistance
More realistically, however, we would expect there to be non-conservative forces acting on the trike, such as rolling resistance and aerodynamic drag. When these forces act upon the trike a simple analytical solution becomes less easy to obtain, at least if the angle of descent is allowed to not remain constant. To parametize the problem we will use an independent variable $s$ that will be distance, in metres, from the summit of the hill or mountain. We will assume this mountain has a road that we are travelling down and $s$ is the distance from the summit along this road. 
	
## Constant gravity
We will first assume that the variation in gravitational acceleration between the base and summit of the hill or mountain is negligible and derive the equation of motion based on this assumption. 
	
### Formulating the Lagrangian
The first step in deriving the equation of motion is to determine the Lagrangian of the system. The Lagrangian contains both the kinetic and potential energy of the system and hence both will be calculated next.

#### Kinetic energy
Substituting:
\begin{align}
v=\dot{s} \label{v}
\end{align}
where $\dot{s}=\dfrac{ds}{dt}$ is the trike's velocity in $\mathrm{m}\cdot \mathrm{s}^{-1}$ into Equation \eqref{T} yields:
\begin{align}
	T=\dfrac{m\dot{s}^2}{2}. \label{Ts}
\end{align}
	
#### Potential energy
The potential energy term will only include the conservative force of gravitation. To compute it, we need to be able to formulate the height of the trike above sea level in terms of $s$. If we have inclination/declination angle data for the hill in question, let us define $\theta(s)$ as the declination angle (positive when downhill) at $s$ metres from the summit. If we form an infinitesimal triangle on the hill/mountain, we will see that $dh = -d s \sin{\theta}$ and hence:
\begin{align}
	\dfrac{dh}{ds} &= -\sin{\theta(s)}.
\end{align}
Integrating from the summit, we obtain that:
\begin{align}
	h(s) &= h_{\mathrm{summit}} - \int_{0}^s \sin{\theta(\sigma)}d\sigma. \label{hs}
\end{align}
Where $h_{\mathrm{summit}}$ is the height of the hill/mountain, or rather its tallest point accessible by road, and we use $\sigma$ as a dummy (integration) variable for distance from the summit. 
	
Hence if gravitational acceleration can be assumed to be constant:
\begin{align}
	V &= mg\left(h_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right). \label{Vcf}
\end{align}

#### Final Lagrangian
The Lagrangian is, by definition:
\begin{align}
	L &= T - V. \label{L}
\end{align}
	
Substituting Equations \eqref{T} and \eqref{Vcf} into Equation \eqref{L}. Therefore our Lagrangian is:
\begin{align}
	L &= \dfrac{m\dot{s}^2}{2} - mg\left(h_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right). \label{Lc}
\end{align}
	
### Non-conservative forces
It is important to note, however, that nonconservative forces --- such as air resistance and rolling resistance --- are a major factor in determining the speed of a tricycle going downhill. To account for these effects, we need to calculate generalized forces. Fortunately, since our problem only has a single independent variable --- that of distance from the summit --- calculating the generalized forces should be straightforward. 
	
Rolling resistance is proportional to the normal force. The normal force could be calculated as:
\begin{align}
N=mg\cos{\theta(s)}. \label{N}
\end{align}
The rolling resistance force is defined as:
\begin{align*}
	F_r &= C_{rr}N. \label{FrN}
\end{align*}
Where $C_{rr}$ is the coefficient of rolling resistance. Substituting Equation \eqref{N} into Equation \eqref{FrN} yields:
\begin{align}
	F_r &= mgC_{rr}\cos{\theta(s)}. \label{Fr}
\end{align}
This force directly opposes motion in the direction of $s$. 
	
As for air resistance (or drag):
\begin{align}
	F_D &= \dfrac{\rho C_D A v^2}{2}.\label{Fdv}
\end{align}
Substituting Equation \eqref{v} into Equation \eqref{Fdv} yields:
\begin{align}
	F_D &= \dfrac{\rho C_D A \dot{s}^2}{2}. \label{Fd}
\end{align}
Where $\rho$ is air density in $\mathrm{kg} \cdot \mathrm{m}^{-3}$, $C_D$ is a dimensionless drag coefficient, $A$ is the cross sectional surface area (in $m^2$) that the area acts on. This is also in the opposite direction to motion. Hence the generalized non-conservative forces acting on our tricycle would be:
\begin{align*}
	Q_s &= -mgC_{rr} \cos{\theta(s)} - \dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
	
### Euler-Lagrange equation with non-conservative forces
The Euler-Lagrange equation with dissipative forces in the $s$ coordinate is:
\begin{align}
	\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{s}} - \dfrac{\partial L}{\partial s} &= Q_s. \label{ELENCF}
\end{align}
First we shall calculate $\dfrac{\partial L}{\partial \dot{s}}$, the canonical momentum:
\begin{align*}
	\dfrac{\partial L}{\partial \dot{s}} &= m\dot{s}.
\end{align*}
Hence:
\begin{align}
	\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{s}} &= m\ddot{s}. \label{ma}
\end{align}
As for the second term in Equation \eqref{ELENCF}, it is:
\begin{align*}
	\dfrac{\partial L}{\partial s} &= mg\sin{\theta(s)}.
\end{align*}
Hence Equation \eqref{ELENCF} becomes:
\begin{align*}
	m\ddot{s} - mg\sin{\theta(s)} &= -mgC_{rr} \cos{\theta(s)} - \dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
Rearranging and factorizing yields:
\begin{align*}
	m\ddot{s} &= mg(\sin{\theta(s)}+C_{rr}\cos{\theta(s)})-\dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
Dividing both sides by $m$ yields:
\begin{align}
	\ddot{s} &= g(\sin{\theta(s)}+C_{rr}\cos{\theta(s)})-\dfrac{\rho C_D A\dot{s}^2}{2m}\label{eq1}
\end{align}
assuming constant gravitational acceleration.
	
### Simplification to first order
\begin{align*}
	\ddot{s} &= \dfrac{d\dot{s}}{ds} \dfrac{ds}{dt} \\
	&= \dfrac{ds}{dt} \dfrac{d\dot{s}}{ds} \\
	&= \dot{s} \dfrac{d\dot{s}}{ds}.
\end{align*}
Hence Equation \eqref{eq1} becomes, after dividing by $m$:
\begin{align*}
	\dot{s} \dfrac{d\dot{s}}{ds} &= g(\sin{\theta(s)} - C_{rr} \cos{\theta(s)}) - \dfrac{\rho C_D A\dot{s}^2}{2m}.
\end{align*}
Let $u=\dot{s}^2$, then:
\begin{align*}
	\dfrac{du}{ds} &= \dfrac{d(\dot{s}^2)}{ds}
\end{align*}
Hence by the chain rule:
\begin{align*}
	\dfrac{du}{ds} &= 2\dot{s}\dfrac{d\dot{s}}{ds}.
\end{align*}
Consequently:
\begin{align*}
	\dfrac{1}{2}\dfrac{du}{ds} &= g(\sin{\theta(s)} - C_{rr} \cos{\theta(s)})  - \dfrac{\rho C_D A u}{2m}.    
\end{align*}
Multiplying both sides by 2 yields:
\begin{align*}
	\dfrac{du}{ds} &= 2g(\sin{\theta(s)} - C_{rr} \cos{\theta(s)})  - \dfrac{\rho C_D A u}{m}.
\end{align*}
### Solution
\begin{align}
	u(s) = C_1 e^{-\frac{C_D As}{m}} + 2ge^{-\frac{C_D As}{m}}\int_0^s e^{\frac{C_D A\sigma}{m}}(\sin{\theta(\sigma)}-C_{rr} \cos{\theta(\sigma)})d \sigma.
\end{align}
The trike starts at the summit with zero velocity, consequently $u(0)=0$ and hence $C_1=0$. This simplifies $u(s)$ to:
\begin{align}
	u(s) = 2ge^{-\frac{C_D As}{m}}\int_0^s e^{\frac{C_D A\sigma}{m}}(\sin{\theta(\sigma)}-C_{rr} \cos{\theta(\sigma)})d \sigma.
\end{align}
	
Square rooting we obtain the velocity:
\begin{align}
    v(s) &= \sqrt{2ge^{-\frac{C_D As}{m}}\int_0^s e^{\frac{C_D A\sigma}{m}}(\sin{\theta(\sigma)}-C_{rr} \cos{\theta(\sigma)})d \sigma}.
\end{align}
We neglected the negative possibility because velocity is positive in the downhill direction, which is the only direction the trike will go when under the action of gravity. 

Below are the velocities achieved in the above model for various $C_D A$ and $C_{rr}$ values when descending Castle Hill. The simulation did not end at 0 metres above sea level; it ended at about 36 metres above sea level, which is approximately at the high end of Gregory Street. $\theta(s)$ values were calculated from QldDem data on elevation. 

~~~
<figure>
<img src="https://fusion809.github.io/images/SageMath/Castle%20hill%20descent%20speed%20constant%20gravity.svg" width="100%">
<figcaption>Velocity while descending Castle Hill, assuming constant gravity.</figcaption>
</figure>
~~~

## Non-constant gravity
Gravitational acceleration depends on the distance from the centre of the source of gravity, so there is some variation in gravitational acceleration between the base and summit of hills and mountains although it is minor. In this section, the equation of motion will be derived while accounting for this variability.

### Formulating the Lagrangian
The kinetic energy we use when studying non-constant gravity is the same as before. Consequently, Equation \eqref{ma} is also true for this Lagrangian.
	
#### Potential energy
If we cannot assume constant gravity, then:
\begin{align}
	V &=\dfrac{GMm}{(r_{\mathrm{base}} + h(s))} - \dfrac{GMm}{r_{\mathrm{summit}}}. 
\end{align}
Factorizing:
\begin{align}
	V &= GMm\left((r_{\mathrm{base}} + h(s))^{-1} - r_{\mathrm{summit}}^{-1}\right).\label{Vnc}
\end{align}
Where $G$ is the gravitational constant and $M$ is the mass of the Earth.
	
Substituting Equation \eqref{hs} into Equation \eqref{Vnc} yields:
\begin{align*}
	V &= GMm\left[\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^{-1}-r_{\mathrm{summit}}^{-1}\right].
\end{align*}
#### Final Lagrangian
Or, if gravitational acceleration cannot be assumed to be constant:
\begin{align*}
	L &= \dfrac{m\dot{s}^2}{2} - GMm\left[\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^{-1}-r_{\mathrm{summit}}^{-1}\right].
\end{align*}
### Non-conservative forces
If gravity cannot be assumed to be constant then:
\begin{align*}
	g=\dfrac{GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2}. \label{g}
\end{align*}
Substituting Equation \eqref{g} into Equation \eqref{Fr} yields:
\begin{align}
	F_r &= \dfrac{GMmC_{rr}\cos{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2}
\end{align}
### Computing Equation \eqref{ELENCF}
The conservative force term of the Euler-Lagrange equation is therefore:
\begin{align}
	\dfrac{\partial L}{\partial s} &= \dfrac{GMm\sin{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2}. \label{Fnc}
\end{align}
The time derivative of generalized momentum is the same as above. As for $Q_s$ it is:
\begin{align}
	Q_s &= -\dfrac{GMmC_{rr}\cos{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} - \dfrac{\rho C_D A \dot{s}^2}{2}. \label{Qs}
\end{align}
Substituting Equations \eqref{Fnc} and \eqref{Qs} into Equation \eqref{ELENCF} yields:
\begin{align*}
	m\ddot{s} - \dfrac{GMm\sin{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} &= -\dfrac{GMmC_{rr}\cos{\theta(s)}}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} - \dfrac{\rho C_D A\dot{s}^2}{2}.
\end{align*}
Rearranging and factorizing:
\begin{align}
	m\ddot{s} &= \dfrac{GMm}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right] - \dfrac{\rho C_D A\dot{s}^2}{2} \\
	\ddot{s} &= \dfrac{GM}{\left(r_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right)^2} \left[\sin{\theta(s)}-C_{rr}\cos{\theta(s)}\right] - \dfrac{\rho C_D A\dot{s}^2}{2m}.\label{eq2}
\end{align}
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
### Solution
In the denominator, we will use a new dummy variable of integration called $\xi$. Hence:
\begin{align}
	u(s) = C_1 e^{-\frac{C_D As}{m}} + 2GMe^{-\frac{C_D As}{m}}\int_0^s e^{\frac{C_D A\sigma}{m}}\dfrac{\sin{\theta(\sigma)}-C_{rr} \cos{\theta(\sigma)}}{\left(r_{\mathrm{summit}}-\int_0^\sigma \sin{\theta(\xi)} d\xi\right)^2} d \sigma.
\end{align}
Yet again $u(0)=0$ means $C_1=0$ and hence:
\begin{align}
	u(s) &= 2GMe^{-\frac{C_D As}{m}}\int_0^s e^{\frac{C_D A\sigma}{m}}\dfrac{\sin{\theta(\sigma)}-C_{rr} \cos{\theta(\sigma)}}{\left(r_{\mathrm{summit}}-\int_0^\sigma \sin{\theta(\xi)} d\xi\right)^2} d \sigma.
\end{align}
Square rooting $u$ to obtain velocity yields:
\begin{align}
    v(s) &= \sqrt{2GMe^{-\frac{C_D As}{m}}\int_0^s e^{\frac{C_D A\sigma}{m}}\dfrac{\sin{\theta(\sigma)}-C_{rr} \cos{\theta(\sigma)}}{\left(r_{\mathrm{summit}}-\int_0^\sigma \sin{\theta(\xi)} d\xi\right)^2} d \sigma}.
\end{align}
Below are the velocities achieved in the above model for various $C_D A$ and $C_{rr}$ values when descending Castle Hill. The simulation did not end at 0 metres above sea level; it ended at about 36 metres above sea level, which is approximately at the high end of Gregory Street. 

~~~
<figure>
<img src="https://fusion809.github.io/images/SageMath/Castle%20hill%20descent%20speed%20non-constant%20gravity.svg" width="100%">
<figcaption>Velocity while descending Castle Hill, assuming non-constant gravity.</figcaption>
</figure>
~~~