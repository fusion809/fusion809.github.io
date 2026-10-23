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

Where: $v$ is velocity; $m$ is the mass of the tricycle (in kilograms), including the mass of anything on it (e.g. its rider); $g$ is gravitational acceleration (in metres per second squared ($\mathrm{m}\cdot \mathrm{s}^{-2}$)); and $h(s)$ is the height (in $m$) of the base of the tricycle above sea level. The $mgh$ value for $V$ assumes constant gravitational acceleration. 
	
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
In Equation \eqref{vnonconst}, we deliberately ignored the negative possibility as the velocity as we only care about the magnitude of the velocity not its direction. Keeping in mind that the tallest point of Castle Hill is not reachable by road, so maximum altitude would be approximately somewhere in the vicinity of 265.9 metres above sea level --- at least based on what Google Fit says of my altitude when I journey there. According to Equations \eqref{vconst} and \eqref{vnonconst}, the maximum velocity of a trike travelling down Castle Hill would be approximately 260 kilometres per hour (kph). 