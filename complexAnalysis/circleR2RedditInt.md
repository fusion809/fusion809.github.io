@def title="Contour integration with pole on the contour itself"
@def tags="Contour integration"

This post is in response to [this thread on Reddit](https://www.reddit.com/r/MathHelp/comments/ioku9v/contour_integration_with_pole_on_the_contour/). 

The original question was essentially to evaluate the integral:

\[
    \oint_C \dfrac{e^{2z}}{(z-2)(z-3)} dz
\]

on the circle defined by $|z| = 2$. For brevity's sake, let's call our integrand $f(z)$. This presents an obvious challenge, as a pole exists on this curve, which means that Cauchy's integral formula and the residue theorem cannot be used as is.

The way around this problem is to modify the curve so that the pole no longer occurs on it. Below is one such modified curve:

~~~
<p>
<img src="/assets/Circle%20with%20radius%202%20and%20enclave.svg" width="80%" alt>
</p>
~~~

As no pole occurs within C' and $f(z)$ is defined along and in C', we know by the Cauchy integral theorem that:

\[
    \oint_{C'} \dfrac{e^{2z}}{(z-2)(z-3)} dz = 0
\]

and splitting $\displaystyle \oint_{C'} f(z) dz$ into subintegrals from A to B and B to A, we get:

\[
    \int_{AB} \dfrac{e^{2z}}{(z-2)(z-3)} dz + \int_{BA} \dfrac{e^{2z}}{(z-2)(z-3)} dz = 0.
\]

Along AB:

$z = 2e^{i\theta}$, $dz = 2ie^{i\theta} d\theta$ and $\delta \leq \theta \leq 2\pi - \delta$, where $\delta = \sin^{-1}\left(\dfrac{\epsilon}{2}\right)$.

Along BA:

$z = 2 + \epsilon e^{i\theta}$, $dz = 2i\epsilon e^{i\theta} d\theta$ and $\theta$ goes from $\dfrac{3\pi}{2}$ to $\dfrac{\pi}{2}$. 

Therefore (3) becomes:

\[
    \int_{\delta}^{2\pi - \delta} \dfrac{e^{2 \cdot 2 e^{i\theta}}}{(2e^{i\theta}-2)(2e^{i\theta}-3)} 2i e^{i\theta} d\theta + \int_{\dfrac{3\pi}{2}}^{\dfrac{\pi}{2}} \dfrac{e^{2 \cdot (2+\epsilon e^{i\theta})}}{(2 + \epsilon e^{i\theta}-2)(2 + \epsilon e^{i\theta}-3)} i \epsilon e^{i\theta} d\theta = 0.
\]

As $\epsilon \rightarrow 0$, the first of these integrals becomes the same as $\displaystyle \oint_C f(z) dz$ and the second integral becomes:

\[
    \begin{aligned}
    \lim_{\epsilon \rightarrow 0} \int_{\dfrac{3\pi}{2}}^{\dfrac{\pi}{2}} \dfrac{e^{2 \cdot (2+\epsilon e^{i\theta})}}{(2 + \epsilon e^{i\theta}-2)(2 + \epsilon e^{i\theta}-3)} i \epsilon e^{i\theta} d\theta &= \lim_{\epsilon \rightarrow 0} ie^4 \int_{\dfrac{3\pi}{2}}^{\dfrac{\pi}{2}} \dfrac{e^{2 \epsilon e^{i\theta}}}{(\epsilon e^{i\theta})(\epsilon e^{i\theta}-1)} \epsilon e^{i\theta} d\theta \\
    &= \lim_{\epsilon \rightarrow 0} ie^4 \int_{\dfrac{3\pi}{2}}^{\dfrac{\pi}{2}} \dfrac{e^{2 \epsilon e^{i\theta}}}{\epsilon e^{i\theta}-1} d\theta \\
    &= ie^4 \int_{\dfrac{3\pi}{2}}^{\dfrac{\pi}{2}} \dfrac{1}{-1} d\theta \\
    &= -ie^4 \int_{\dfrac{3\pi}{2}}^{\dfrac{\pi}{2}} d\theta \\
    &= ie^4 \int_{\dfrac{\pi}{2}}^{\dfrac{3\pi}{2}} d\theta \\
    &= ie^4 \pi.
    \end{aligned}
\]

(4) therefore becomes:

\[
    \begin{aligned}
        \oint_C \dfrac{e^{2z}}{(z-2)(z-3)} dz + ie^4 \pi &= 0.
    \end{aligned}
\]

Hence our original integral is:

\[
    \oint_C \dfrac{e^{2z}}{(z-2)(z-3)} dz = - ie^4 \pi.
\]