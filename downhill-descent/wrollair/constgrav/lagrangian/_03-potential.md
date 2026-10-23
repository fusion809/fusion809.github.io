#### Potential energy
The potential energy term will only include the conservative force of gravitation. To compute it, we need to be able to formulate the height of the trike above sea level in terms of $s$. If we have inclination/declination angle data for the hill in question, let us define $\theta(s)$ as the declination angle (positive when downhill) at $s$ metres from the summit. If we form an infinitesimal triangle on the hill/mountain, we will see that $dh = -d s \sin{\theta}$ and hence:
\begin{align}
	\dfrac{dh}{ds} &= -\sin{\theta(s)}.
\end{align}
Integrating from the summit, we obtain that:
\begin{align}
	h(s) &= h_{\mathrm{summit}} - \int_{0}^s \sin{\theta(\sigma)}d\sigma. \label{hs}
\end{align}
Where $h_{\mathrm{summit}}$ is the height of the hill/mountain, or rather its tallest point accessible by road, and we use $\sigma$ as a dummy (integration) variable for distance from the summit. 
	
Hence if gravitational acceleration can be assumed to be constant:
\begin{align}
	V &= mg\left(h_{\mathrm{summit}} - \int_0^s \sin{\theta(\sigma)}d\sigma\right). \label{Vcf}
\end{align}