## General relativity
In general relativity, we can use the Einstein-Hilbert action to derive the Einstein field equations (EFE). The EFE are a series of partial differential equations that can be used to calculate the components of the metric tensor, $g_{\mu\nu}$. The metric tensor is a type of matrix with sixteen total components that allow us to define proper time in curved spacetime. It has sixteen components because $\mu$ and $\nu$ can take on four different values each, corresponding to the four dimensions of spacetime (e.g. in Cartesian coordinates they are $t$, $x$, $y$ and $z$), so that gives us $4\times 4=16$ combinations of possible $\mu$ and $\nu$ values. Proper time is analogous to distance in spacetime and is the time that passes locally for an observer. 

I should mention, in case you are unfamiliar, that in relativity, the passage of time itself is relative. Objects travelling close to the speed of light or in powerful gravitational fields have time pass substantially slower for them than those at lower velocities and in weaker gravitational fields. 

The Einstein-Hilbert action is:
\begin{align}
S &= \dfrac{c^4}{16\pi G} \int_{\mathbb{R}^4} R \sqrt{\left|\det{g_{\mu\nu}}\right|} d^4 x.
\end{align}
One major difference you may notice between the classical action and the Einstein-Hilbert action is that the classical action is integrated over time whereas the Einstein-Hilbert action is integrated over spacetime (that is what the $d^x$ and $\mathbb{R}^4$ is about). $R$ is the Ricci scalar which is calculated based on $g_{\mu\nu}$. 

Another application of action in general relativity is that it can be used to calculate the path an object will take in curved spacetime. As the path the object follows will be the path that maximizes its proper time. The action for this calculation is:

\begin{align}
S &= \int_{\lambda_0}^{\lambda_f} \sqrt{\Big| g_{\mu\nu} \dfrac{dx^{\mu}}{d\lambda} \dfrac{dx^{\nu}}{d\lambda}\Big|}d\lambda.
\end{align}

Where $\lambda$ is an arbitrary parameter. Often we choose it to be proper time. Setting $\delta S=0$ yields the geodesic equations which are a series of ordinary differential equations that maximize proper time. 