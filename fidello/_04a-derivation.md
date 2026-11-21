## Derivation
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
Assuming our Hamiltonian does not explicitly depend on $t$, we can also write $dH$ as:
\begin{align}
dH &= \sum_{i=1}^n \dfrac{\partial H}{\partial q_i} dq_i + \dfrac{\partial H}{\partial p_i}dp_i. \label{dH4}
\end{align}
As this is just the change in $H$ caused directly by the change in $q_i$ plus the change in $H$ that is directly caused by the change in $p_i$. 
Equating like terms in the right-hand sides of Equations \eqref{dH3} and \eqref{dH4} yields:
\begin{align*}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} \\
\dot{q}_i &= \dfrac{\partial H}{\partial p_i}.
\end{align*}
Which are Equations \eqref{Heq1} and \eqref{Heq2}, which are known as Hamilton's equations. I should mention this is not the only way to derive Hamilton's equations; it is also possible to derive them directly from the stationary action principle, for instance. 
