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