# Noether's theorem
Noether's theorem essentially states that for every continuous symmetry in the Lagrangian, there is a corresponding conservation law. For instance, time translation symmetry of the Lagrangian yields the conservation of energy. To show this, we will assume our Lagrangian does not explicitly depend on $t$. That is, $t$ itself does not appear directly in our Lagrangian. The only time dependence $L$ has on $t$ is through the dependence of $q_i$ and $\dot{q}_i$ on time. If the Lagrangian did explicitly depend on time, we would expect it not to have time translation symmetry. Say we increase time according to the transformation: $t\rightarrow t+\epsilon$. 

\begin{align}
L(q_i(t+\epsilon), \dot{q}_i(t+\epsilon)) &\approx L(q_i(t), \dot{q}_i(t)) + \dfrac{dL}{dt} \epsilon \\
&= L(q_i(t), \dot{q}_i(t)) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \dfrac{dq_i}{dt} \epsilon + \dfrac{\partial L}{\partial \dot{q}_i} \dfrac{d\dot{q}_i}{dt} \epsilon \label{Lshift}..
\end{align}
Subtracting $L(q_i(t), \dot{q}_i(t))$ from both sides of Equation \eqref{Lshift} yields:
\begin{align}
L(q_i(t+\epsilon), \dot{q}_i(t+\epsilon)) - L(q_i(t), \dot{q}_i(t)) &\approx \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \dfrac{dq_i}{dt} \epsilon + \dfrac{\partial L}{\partial \dot{q}_i} \dfrac{d\dot{q}_i}{dt} \epsilon. \label{diffL}
\end{align}
We will call the left-hand side of Equation \eqref{diffL} $\Delta L$, hence Equation \eqref{diffL} becomes:
\begin{align}
\Delta L &\approx \sum_{i=1}^n \dfrac{\partial L}{\partial q_i} \dfrac{dq_i}{dt} \epsilon + \dfrac{\partial L}{\partial \dot{q}_i} \dfrac{d\dot{q}_i}{dt} \epsilon. \label{DeltaL1}
\end{align}
Rearranging Equations \eqref{ELE} yields that:
\begin{align}
\dfrac{\partial L}{\partial q_i} = \dfrac{d}{dt}\left(\dfrac{\partial L}{\partial \dot{q}_i}\right). \label{ELEnh}
\end{align}
If you are wondering what about if we have non-conservative forces? Would that not mean Equations \eqref{ELEnh} would not hold as there would be a missing generalized force term? Yes and in that event the Lagrangian would indeed change in time and conservation of energy would not hold, as our system would be losing energy to those non-conservative forces. 

Substituting Equations \eqref{ELEnh} into Equation \eqref{DeltaL1} yields:
\begin{align}
\Delta L &\approx \sum_{i=1}^n \epsilon\left(\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{q}_i}\dot{q}_i + \dfrac{\partial L}{\partial \dot{q}_i}\ddot{q}_i\right). \label{DeltaL2}
\end{align}
Substituting Equations \eqref{canmom} into Equation \eqref{DeltaL2} yields:
\begin{align}
\Delta L &\approx \epsilon\sum_{i=1}^n \left(\dfrac{dp_i}{dt} \dot{q}_i + p_i \ddot{q}_i\right)\nonumber\\
&= \epsilon\sum_{i=1}^n \left(\dot{p}_i \dot{q}_i + p_i \ddot{q}_i\right). \label{DeltaL3}
\end{align}
We can rewrite the bracketed expression on the right-hand side of Equation \eqref{DeltaL3} using the product rule of differentiation as:
\begin{align*}
\Delta L &\approx \epsilon\sum_{i=1}^n \dfrac{d}{dt}\left(p_i \dot{q}_i\right) \\
&= \epsilon\dfrac{d}{dt}\sum_{i=1}^n p_i \dot{q}_i.
\end{align*}
If $L$ does not change in time, $\Delta L=0$, and hence:
\begin{align}
\epsilon \dfrac{d}{dt} \sum_{i=1}^n p_i \dot{q}_i &= 0. \label{epsdTdt}
\end{align}
Dividing both sides of Equations \eqref{epsdTdt} by $\epsilon$ yields:
\begin{align*}
\dfrac{d}{dt} \sum_{i=1}^n p_i \dot{q}_i &= 0.
\end{align*}
$L$ does not change, so $\dfrac{dL}{dt}=0$, too. Hence:
\begin{align}
\dfrac{d}{dt} \left(\sum_{i=1}^n p_i \dot{q}_i - L\right) &= 0.\label{constE}
\end{align}
Interestingly, $\displaystyle\sum_{i=1}^n p_i \dot{q}_i=2T$, that is it is twice the kinetic energy of our system. Hence if we substitute Equation \eqref{L} into Equation \eqref{constE}, we obtain:
\begin{align}
\dfrac{d}{dt} \left(2T-(T-V)\right) &= 0 \nonumber\\
\dfrac{d}{dt} \left(T+V\right) &= 0.\label{constE2}
\end{align}
Hence substituting Equation \eqref{H} into Equation \eqref{constE2} yields:
\begin{align}
\dfrac{dH}{dt} &= 0.\label{constE3}
\end{align}
In other words, we obtain that the Hamiltonian &mdash; the total energy of the system &mdash; does not change, hence it is conserved.