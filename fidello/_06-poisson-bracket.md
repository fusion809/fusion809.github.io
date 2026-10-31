# Poisson brackets
Poisson brackets provide us an elegant and powerful way to rewrite Hamilton's equations. If we have $n$ degrees of freedom and two functions $f$ and $g$ that depend on $p_i$ and $q_i$:[^1] 

\begin{align}
\{f, g\} &= \sum_{i=1}^n \left( \dfrac{\partial f}{\partial q_i} \dfrac{\partial g}{\partial p_i} - \dfrac{\partial f}{\partial p_i}\dfrac{\partial g}{\partial q_i} \right). \label{Poissonfg}
\end{align}

Poisson brackets also have several interesting properties, such as antisymmetry:
\begin{align}
\{A, B\} &= -\{B, A\}. \label{asym}
\end{align}
Further, assuming each of the coordinates of the phase space are independent, it follows that:
\begin{align}
\{q_k, p_k\} &= \sum_{i=1}^n \left( \dfrac{\partial q_k}{\partial q_i} \dfrac{\partial p_k}{\partial p_i} - \dfrac{\partial q_k}{\partial p_i} \dfrac{\partial p_k}{\partial q_i} \right) \nonumber\\
&= \sum_{i=1}^n \delta_{ik}\delta_{ik} \nonumber\\
&= 1. \label{qkpk}
\end{align}
Where $\delta_{ik}$ is the Kronecker delta function. It is equal to zero if $i\neq k$ and one if $i=k$. 
\begin{align}
\{q_j, q_k\} &= \sum_{i=1}^n \left( \dfrac{\partial q_j}{\partial q_i} \dfrac{\partial q_k}{\partial p_i} - \dfrac{\partial q_j}{\partial p_i} \dfrac{\partial q_k}{\partial q_i} \right) \nonumber\\
&= 0 \label{qjqk}\\
\{p_j, p_k\} &= \sum_{i=1}^n \left(\dfrac{\partial p_j}{\partial q_i}\dfrac{\partial p_k}{\partial p_i} - \dfrac{\partial p_j}{\partial p_i}\dfrac{\partial p_k}{\partial q_i}\right) \nonumber\\
&= 0. \label{pjpk}
\end{align}
The quality of two functions, say $f$ and $g$, having a Poisson bracket that equals zero is phrased as $f$ and $g$  **commuting** under the Poisson bracket. 

If a function $f(q_i, p_i)$ is not explicitly dependent on $t$ and only depends on $t$ indirectly through $q_i$ and $p_i$, then, via the chain rule:
\begin{align}
\dfrac{df}{dt} &= \sum_{i=1}^n \dfrac{\partial f}{\partial q_i}\dfrac{d q_i}{dt} + \dfrac{\partial f}{\partial p_i}\dfrac{dp_i}{dt}.
\end{align}
We can use Hamilton's equations to replace $\dfrac{dq_i}{dt}$ and $\dfrac{dp_i}{dt}$:
\begin{align}
\dfrac{df}{dt} &= \sum_{i=1}^n \dfrac{\partial f}{\partial q_i}\dfrac{\partial H}{\partial p_i} - \dfrac{\partial f}{\partial p_i}\dfrac{\partial H}{\partial q_i}.
\end{align}
Interestingly, the right-hand side is equal to the Poisson bracket of $f$ with $H$, in other words:
\begin{align}
\dfrac{df}{dt} &= \{f, H\}.
\end{align}
Hence commuting with the Hamiltonian under the Poisson bracket is a sufficient condition for a quantity that does not explicitly depend on time to have a zero time derivative. A zero time derivative means the quantity is static over time which means it is conserved.

This quality makes it unsurprising that Hamilton's equations can be rewritten as:

\begin{align}
\dot{q}_k &= \{q_k, H\}\label{qdotk}\\
\dot{p}_k &= -\{p_k, H\}. \label{pdotk}
\end{align}