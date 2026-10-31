## Non-conservative forces
I guess next I should explain how we account for non-conservative forces. Essentially we replace the $0$ in Equations \eqref{ELE} with the generalized force $Q_i$. It is defined as:
\begin{align}
    Q_i &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \textbf{r}_j}{\partial q_i}. \label{Qir}
\end{align}
Where $\textbf{F}_j$ represents the vector of the $j$th non-conservative force we wish to account for in our equations of motion. As you can see by the summation bounds, we have $m$ non-conservative forces we wish to account for. $\textbf{r}_j$ is the position vector for force $j$. $\textbf{F}_j\cdot\dfrac{\partial \textbf{r}_j}{\partial q_i}$ denotes the dot product between $\textbf{F}_j$ and $\dfrac{\partial \textbf{r}_j}{\partial q_i}$. Alternatively, to calculate $Q_i$ we can use velocities instead of positions:
\begin{align}
    Q_i &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \dot{\textbf{r}}_j}{\partial \dot{q}_i}.\label{Qirdot}
\end{align}
$\dot{\textbf{r}}_j$ is the time derivative of the position vector (that is, the velocity vector) for force $j$. Hence the more generalized form of Equations \eqref{ELE} is:
\begin{align}
    \hspace{4cm}\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i} - \dfrac{\partial L}{\partial q_i} &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \textbf{r}_j}{\partial q_i}, \hspace{2cm} i\in\{1,2,3,...,n\}. \label{ELENCF}
\end{align}
The velocity-based definition of $Q_i$ can also be used, which yields:
\begin{align}
    \hspace{4cm}\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i} - \dfrac{\partial L}{\partial q_i} &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \dot{\textbf{r}}_j}{\partial \dot{q}_i}, \hspace{2cm} i\in\{1,2,3,...,n\}. \label{ELENCFrdot}
\end{align}