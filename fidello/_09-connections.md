# Connections to other branches of physics
## Special relativity
In special relativity, we can make relativistic corrections to kinetic energy in our Lagrangian and Hamiltonians. Namely, we can define it as:

\begin{align*}
T &= \dfrac{mc^2}{\sqrt{1-\dfrac{v^2}{c^2}}} - mc^2 \\
&= mc^2\left(\dfrac{1}{\sqrt{1-\dfrac{v^2}{c^2}}}-1\right).
\end{align*}

## General relativity
In general relativity, we can use the Einstein-Hilbert action to derive the Einstein field equations (EFE). The EFE are a series of partial differential equations that can be used to calculate the components of the metric tensor, $g_{\mu\nu}$. The metric tensor is a type of matrix with sixteen total components that allow us to define proper time in curved spacetime. It has sixteen components because $\mu$ and $\nu$ can take on four different values each, corresponding to the four dimensions of spacetime (e.g. in Cartesian coordinates they are $t$, $x$, $y$ and $z$), so that gives us $4\times 4=16$ combinations of possible $\mu$ and $\nu$ values. Proper time is analogous to distance in spacetime and is the time that passes locally for an observer. 

I should mention, in case you are unfamiliar, that in relativity, the passage of time itself is relative. Objects travelling close to the speed of light or in powerful gravitational fields have time pass substantially slower for them than those at lower velocities and in weaker gravitational fields. 

The Einstein-Hilbert action is:
\begin{align}
S &= \dfrac{c^4}{16\pi G} \int_{\mathbb{R}^4} R \sqrt{\left|\det{g_{\mu\nu}}\right|} d^4 x.
\end{align}
One major difference you may notice between the classical action and the Einstein-Hilbert action is that the classical action is integrated over time whereas the Einstein-Hilbert action is integrated over spacetime. $R$ is the Ricci scalar which is calculated based on $g_{\mu\nu}$. 

Another application of action in general relativity is that it can be used to calculate the path an object will take in curve spacetime. As the path the object follows will be the path that maximizes its proper time. The action for this calculation is:

\begin{align}
S &= \int_{\lambda_0}^{\lambda_f} \sqrt{\Big| g_{\mu\nu} \dfrac{dx^{\mu}}{d\lambda} \dfrac{dx^{\nu}}{d\lambda}\Big|}d\lambda.
\end{align}

Where $\lambda$ is an arbitrary parameter. Often we choose it to be proper time. Setting $\delta S=0$ yields the geodesic equations which are a series of ordinary differential equations that maximize proper time. 

## Quantum mechanics
We have seen in this article that there are multiple different formulations of classical mechanics. Similarly, there are multiple different formulations of quantum mechanics. One formulation of quantum mechanics is the path-integral formulation developed by Richard Feynmann. In it, we expect our particle or wave to explore all possible paths and we integrate over all those paths in order to compute the probability amplitude of our particle. The exploration of all possible paths is a consequence of the randomness of quantum mechanics. Specifically, we compute the integral:

\begin{align}
K(\mathbf{x}_f, t_f; \mathbf{x}_i, t_i) &= \int_{\mathbf{x}(t_i)=\mathbf{x}_i}^{\mathbf{x}(t_f)=\mathbf{x}_f} e^{\frac{iS(\mathbf{x})}{\hbar}} \mathcal{D}\mathbf{x}.
\end{align}

Where $K$ is the propagator, which specifies the probability amplitude for travelling from $\mathbf{x}_i$ at time $t_i$ to $\mathbf{x}_f$ at time $t_f$. $S$ is the action of the particle, $\hbar$ is the reduced Planck's constant, $i$ is the imaginary number and $\mathcal{D}\mathbf{x}$ refers to the fact we are integrating over all possible paths. The probability amplitude is a quantity that, if you take the absolute value of, and square this quantity you get the probability of the particle travelling from $\mathbf{x}_i$ at time $t_i$ to $\mathbf{x}_f$ at time $t_f$. 

One major strength of the path-integral formulation of quantum mechanics is that it provides an explanation for how quantum mechanics gives rise to classical mechanics on larger scales. 

I am going to have to assume you are familiar with complex numbers. In which case, you will know that $e^{\frac{iS}{\hbar}}$ is like a 2D unit vector as it points in particular direction in the $xy$ plane and has unit length. Specifically, it points in the direction $\dfrac{S}{\hbar}$ radians counterclockwise from the positive $x$ axis. Each integer multiple of $2\pi$, which corresponds to a complete revolution around the Cartesian plane, that is present in $\dfrac{S}{\hbar}$ makes no difference to its overall direction. 

Interestingly, when we integrate over all possible paths on classical mechanics scales we find that the $e^{\frac{iS}{\hbar}}$, for most paths, cancel each other out. To see why, we must first start with the fact that $\hbar$ is a very tiny number &mdash; it is approximately $1.055\times 10^{-34}J\cdot s$. Consequently, $\dfrac{S}{\hbar}$ is a massive number when we are analysing classical mechanical systems. As typically in these scenarios $S$ is $\gt0.01J\cdot s$. Consequently, even a relatively small difference in $S$, can lead to a majorly different direction for $e^{\frac{iS}{\hbar}}$.

This leads to the actions of nearby paths often being in opposite directions and cancelling each other when we calculate the path integral. The one exception is the path corresponding to the stationary action. Remember the equation $\delta S=0$ that we solve to find the path of our classical mechanical system? It says that a tiny change to the stationary action path should make no change (to first order) to the action and consequently, when we integrate $e^{\frac{iS}{\hbar}}$ for that path and its neighbours, they do not cancel each other out.

Another connection between what we have covered and quantum mechanics is in the operator formalism of quantum mechanics. In it, we represent the state of our physical system using a function called the wavefunction which is denoted as $\psi$, which can also thought of as a vector. We can calculate the wavefunction using Schrödinger's equation:

\begin{align}
i\hbar \dfrac{\partial \psi}{\partial t} &= \left(-\dfrac{\hbar^2}{2m}\Delta + V(\mathbf{x})\right)\psi.\label{Schrod}
\end{align}

$V(\mathbf{x})$ is the potential energy and $\Delta$ is the Laplacian which is essentially just $\dfrac{\partial^2}{\partial x^2} + \dfrac{\partial^2}{\partial y^2}+\dfrac{\partial^2}{\partial z^2}$. The right-hand side of Equation \eqref{Schrod} can also be represented by $\hat{H}\psi$ where $\hat{H}$ is the Hamiltonian operator and calculates, from the wavefunction, the total energy of our system.

## Quantum field theory
Lagrangians, Hamiltonians and actions also appear in quantum field theory, including The Standard Model of particle physics. In fact, if you want to see a dauntingly long and complex Lagrangian, Google The Standard Model Lagrangian.