# Hamiltonian mechanics
In Hamiltonian mechanics, we define the Hamiltonian as:
\begin{align}
H &= T + V.\label{H}
\end{align}
Which is the total energy of our system. I should say that this is *usually* the form of the Hamiltonian, more generally it is defined as:
\begin{align}
H &\triangleq \sum_{i=1}^n p_i \dot{q}_i - L. \label{HLeq}
\end{align}
Equation \eqref{HLeq} is called the Legendre transformation of the Lagrangian. Another difference from Lagrangian mechanics is that we define $H$ in terms of $p_i$ and $q_i$ instead of $q_i$ and $\dot{q}_i$. The equations of motion with no non-conservative forces are:
\begin{align}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} \label{Heq1}\\
\dot{q}_i &= \dfrac{\partial H}{\partial p_i}.\label{Heq2}
\end{align}
Equations \eqref{Heq1} and \eqref{Heq2} are often called Hamilton's equations. To account for non-conservative forces, we merely add the corresponding generalized force terms to the right-hand side of Equation \eqref{Heq1}. In other words:
\begin{align}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} + Q_i. \label{HNConseq1}
\end{align}
This creates a new kind of space that we can visualize our solution on that has dimensions $q_i$ and $p_i$. This type of space is known as phase space.