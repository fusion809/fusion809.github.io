# Hamiltonian mechanics
In Hamiltonian mechanics, we define the Hamiltonian as:
\begin{align}
H &= T + V.\label{H}
\end{align}
Which is the total energy of our system. I should say that this is *usually* the form of the Hamiltonian, more generally it is defined as:
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
We can derive these equations of motion by considering an infinitesimal change to our path and its effect on $H$, which we will denote as $dH$:
\begin{align}
dH &= \sum_{i=1}^n d(p_i\dot{q}_i) - dL \label{dH1}\\
&= \sum_{i=1}^n dp_i \dot{q}_i + p_i d\dot{q}_i - dL \label{dH2}\\
&= \sum_{i=1}^n \dot{q}_i dp_i + p_i d\dot{q}_i - \dfrac{\partial L}{\partial q_i} dq_i - \dfrac{\partial L}{\partial \dot{q}_i}d\dot{q}_i.\nonumber
\end{align}
To go from Equations \eqref{dH1} to \eqref{dH2}, we used the product rule of differentiation, in case you were wondering. Next we use the definition of canonical momentum:
\begin{align*}
dH &= \sum_{i=1}^n \dot{q}_i dp_i + p_i d\dot{q}_i - \dfrac{\partial L}{\partial q_i} dq_i - p_i d\dot{q}_i \\
&= \sum_{i=1}^n \dot{q}_i dp_i - \dfrac{\partial L}{\partial q_i} dq_i \\
&= \sum_{i=1}^n - \dfrac{\partial L}{\partial q_i} dq_i + \dot{q}_i dp_i.
\end{align*}
Next we can use Equations \eqref{ELE} to simplify our $\dfrac{\partial L}{\partial q_i}$ term. As remember this, according to Equations \eqref{ELE} is equal to $\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{q}_i}$ which can also be written as $\dot{p}_i$:
\begin{align}
dH &= \sum_{i=1}^n \dot{q}_i dp_i - \dot{p}_i dq_i. \label{dH3}
\end{align}
We can also write $dH$ as:
\begin{align}
dH &= \sum_{i=1}^n \dfrac{\partial H}{\partial q_i} dq_i + \dfrac{\partial H}{\partial p_i}dp_i. \label{dH4}
\end{align}
Equating like terms in the right-hand sides of Equations \eqref{dH3} and \eqref{dH4} yields:
\begin{align*}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} \\
\dot{q}_i &= \dfrac{\partial H}{\partial p_i}.
\end{align*}
Which are Equations \eqref{Heq1} and \eqref{Heq2}, which are known as Hamilton's equations. I should mention this is not the only way to derive Hamilton's equations; it is also possible to derive them directly from the stationary action principle, for instance. 

To account for non-conservative forces, we merely add the corresponding generalized force terms to the right-hand side of Equation \eqref{Heq1}. In other words:
\begin{align}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} + Q_i. \label{HNConseq1}
\end{align}
This creates a new kind of space that we can visualize our solution on that has dimensions $q_i$ and $p_i$. This type of space is known as phase space.