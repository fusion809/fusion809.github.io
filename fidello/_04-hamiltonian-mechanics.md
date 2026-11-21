# Hamiltonian mechanics
In Hamiltonian mechanics, we can *usually* define the Hamiltonian as:
\begin{align}
H &= T + V.\label{H}
\end{align}
Which is the total energy of our system. More generally we define the Hamiltonian as:
\begin{align}
H &\stackrel{\text{def}}{=} \sum_{i=1}^n p_i \dot{q}_i - L. \label{HLeq}
\end{align}
Equation \eqref{HLeq} is called the Legendre transformation of the Lagrangian. So when are Equations \eqref{H} and \eqref{HLeq} not equivalent? The answer is when our potential energy depends, at least in part, on the generalized velocities. Another difference from Lagrangian mechanics is that we define $H$ in terms of $p_i$ and $q_i$ instead of $q_i$ and $\dot{q}_i$. 

So what is the Legendre transformation? It takes a function, say $f(x)$, and it uses its tangent lines and their $y$ intercepts to create a new function, $g(p)$ that incorporates all the information that $f(x)$ contained. [Here](https://www.youtube.com/watch?v=P2HZelQm7Lw) is a great video that explains the Legendre transformation and its derivation. 

The equations of motion with no non-conservative forces are:
\begin{align}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} \label{Heq1}\\
\dot{q}_i &= \dfrac{\partial H}{\partial p_i}.\label{Heq2}
\end{align}