@def title="Evaluating the integral of sine mx over x times x to the eighth plus a to the eighth"
@def tags=["Contour integration"]

The integral we wish to evaluate is $I = \displaystyle \int_{-\infty}^{\infty} \dfrac{\sin{mx}}{x(x^8+a^8)} dx$, where $a, m \gt 0$ and are real.

\toc

# Curve and integrand hunting
The best curve to use to evaluate this integral is likely this one:

~~~
<img src="/assets/Infinite semicircle with infinitesimal interior semicircle.svg" width="60%">
~~~

We will use this contour integral: $\displaystyle \oint_C \dfrac{e^{imz}}{z(z^8+a^8)} dz$. 

# Pole hunting and residue evaluation
The integrand of which has poles where:

\[
    \begin{aligned}
    z^8 + a^8 &= 0 \\
    z^8 &= -a^8 \\
    z^8 &= a^8 \cdot (-1) \\
    z^8_n &= a^8 e^{\pi i (1+2n)},\hspace{0.1cm} n\in \mathbb{Z} \\
    z_n &= ae^{\dfrac{\pi i}{8} (1+2n)}.
    \end{aligned}
\]

$n=0,1,2,3,4,5,6,7$ gives us our only unique poles (which are all first order). Of these, only $n=0,1,2,3$ is within the region enclosed by C. Let $r_n$ be the residue of our integrand at $z_n$. 

\[
    \begin{aligned}
    r_n &= \lim_{z\rightarrow z_n} \dfrac{e^{imz} (z-z_0)}{z(z^8+a^8)} \\
    &= \lim_{z\rightarrow z_n} \dfrac{e^{imz}}{z\cdot 8z^7} \\
    &= \dfrac{e^{imz_n}}{8z^8_n} \\
    &= -\dfrac{e^{imz_n}}{8a^8}.
    \end{aligned}
\]

Let us do some working to simplify $r_0 = -\dfrac{e^{imz_0}}{8a^8}.$ $imz_0 = ame^{\dfrac{5\pi i}{8}} = -am\sin{\dfrac{\pi}{8}}+ami\cos{\dfrac{\pi}{8}}$. Let $p=\sin{\dfrac{\pi}{8}}$ and $q=\cos{\dfrac{\pi}{8}}$. Therefore:

\[
    \begin{aligned}
    r_0 &= -\dfrac{e^{-amp+amiq}}{8a^8} \\
    &= -\dfrac{e^{-amp}e^{amiq}}{8a^8}.\\
    r_1 &= -\dfrac{e^{imz_1}}{8a^8}.
    \end{aligned}
\]

Let us do some working to simplify this last line. $imz_1 = ame^{\dfrac{7\pi i}{8}} = -am\cos{\dfrac{\pi}{8}}+ami\sin{\dfrac{\pi}{8}} = -amq+amip$. Therefore:

\[
    \begin{aligned}
    r_1 &= -\dfrac{e^{-amq+amip}}{8a^8} \\
    &= -\dfrac{e^{-amq}e^{amip}}{8a^8}.\\
    r_2 &= -\dfrac{e^{imz_2}}{8a^8}.
    \end{aligned}
\]

Working on this last line. $imz_2 = ame^{\dfrac{9\pi i}{8}} = -amq-amip$. Therefore:

\[
    \begin{aligned}
    r_2 &= -\dfrac{e^{-amq-amip}}{8a^8} \\
    &= -\dfrac{e^{-amq}e^{-amip}}{8a^8}.\\
    r_3 &= -\dfrac{e^{imz_3}}{8a^8}.
    \end{aligned}
\]

Working on this last line. $imz_3 = ame^{\dfrac{11\pi i}{8}} = -amp-amiq$. Therefore:

\[
    \begin{aligned}
    r_3 &= -\dfrac{e^{-amp-amiq}}{8a^8} \\
    &= -\dfrac{e^{-amp}e^{-amiq}}{8a^8}.
    \end{aligned}
\]

Thus the sum of the residues of the poles within the region enclosed by C is:

\[
  \begin{aligned}
  \sum_{n=0}^3 r_n &= -\dfrac{1}{8a^8} \left(e^{-amp}e^{amiq} + e^{-amq}e^{amip} + e^{-amq}e^{-amip} + e^{-amp}e^{-amiq}\right) \\
  &= -\dfrac{1}{8a^8} \left(e^{-amp}\left(e^{amiq}+e^{-amiq}\right)+e^{-amq}\left(e^{amip}+e^{-amip}\right)\right) \\
  &= -\dfrac{1}{8a^8} \left(2e^{-amp}\cos{amq} + 2e^{-amq}\cos{amp}\right) \\
  &= -\dfrac{1}{4a^8} \left(e^{-amp}\cos{amq} + e^{-amq}\cos{amp}\right).
  \end{aligned}  
\]

# Contour integral evaluation
Therefore our contour integral is:

\[
    \begin{aligned}
    \oint_C \dfrac{e^{imz}}{z(z^8+a^8)} dz &= 2\pi i \cdot -\dfrac{1}{4a^8} \left(e^{-amp}\cos{amq} + e^{-amq}\cos{amp}\right) \\
    &= -\dfrac{\pi i}{2a^8} \left(e^{-amp}\cos{amq} + e^{-amq}\cos{amp}\right).
    \end{aligned}
\]

# Subintegral evaluation
Splitting our contour integral into subintegrals, we obtain:

\[
    \begin{aligned}
        \oint_C \dfrac{e^{imz}}{z(z^8+a^8)} dz &= \int_{AB} + \int_{BD} + \int_{DE} + \int_{EA}
    \end{aligned}
\]

where each integral on the right-hand side of course shares the same integrand as that on the left.

Along AB, $z=x$, $dz=dx$ and $x$ goes from $-R$ to $-\epsilon$. 

\[
    \begin{aligned}
    \int_{AB} = \int_{-R}^{-\epsilon} \dfrac{e^{imx}}{x(x^8+a^8)} dx.
    \end{aligned}
\]

When $\epsilon \rightarrow 0$ and $R\rightarrow \infty$:

\[
    \int_{AB} = \int_{-\infty}^0 \dfrac{e^{imx}}{x(x^8+a^8)} dx.
\]

Along BD, $z=\epsilon e^{i\theta}$, $dz=i\epsilon e^{i\theta} d\theta$ and $\theta$ goes from $\pi$ to $0$.

\[
    \begin{aligned}
    \int_{BD} &= \int_{\pi}^0 \dfrac{e^{im\epsilon e^{i\theta}}}{\epsilon e^{i\theta}(\epsilon^8 e^{8i\theta}+a^8)} i\epsilon e^{i\theta} d\theta \\
    &= -i \int_0^{\pi} \dfrac{e^{im\epsilon e^{i\theta}}}{\epsilon^8 e^{8i\theta}+a^8} d\theta.
    \end{aligned}
\]

When $\epsilon \rightarrow 0$:

\[
    \begin{aligned}
    \int_{BD} &= -i \int_0^{\pi} \dfrac{d\theta}{a^8} \\
    &= -\dfrac{\pi i}{a^8}.
    \end{aligned}
\]

Along DE: $z=x$, $dz=dx$ and $x$ goes from $\epsilon$ to $R$.

\[
    \begin{aligned}
    \int_{DE} &= \int_{\epsilon}^R \dfrac{e^{imx}}{x(x^8+a^8)} dx.
    \end{aligned}
\]

When $\epsilon \rightarrow 0$ and $R\rightarrow \infty$:

\[
    \int_{DE} = \int_{0}^{\infty} \dfrac{e^{imx}}{x(x^8+a^8)} dx.
\]

Finally, along EA, $z=Re^{i\theta}$, $dz=iRe^{i\theta}$ and $\theta$ goes from $0$ to $pi$. Therefore:

\[
    \begin{aligned}
    \int_{EA} &= \int_{0}^{\pi} \dfrac{e^{imR e^{i\theta}}}{R e^{i\theta}(R^8 e^{8i\theta}+a^8)} iR e^{i\theta} d\theta \\
    &= i \int_0^{\pi} \dfrac{e^{imR \cos{\theta}-mR\sin{\theta}}}{R^8 e^{8i\theta}+a^8} d\theta \\
    &= i \int_0^{\pi} \dfrac{e^{imR \cos{\theta}}e^{-mR\sin{\theta}}}{R^8 e^{8i\theta}+a^8} d\theta.
    \end{aligned}
\]

For $0\leq \theta \leq \pi$ (our bounds of integration), $e^{-mR\sin{\theta}} \leq 1$, and $|e^{imR\cos{\theta}}| = 1$ provided $m \geq 0$. Therefore an upper bound on our integral is:

\[
    \begin{aligned}
    \int_{EA} \leq i \int_0^{\pi} \dfrac{d\theta}{R^8 e^{8i\theta}+a^8}
    \end{aligned}
\]

The magnitude of that denominator is:

\[
    \begin{aligned}
    |R^8e^{8i\theta} + a^8| &= \sqrt{(R^8 \cos{8\theta}+a^8)^2 + R^{16} \sin^2{8\theta}} \\
    &= \sqrt{R^{16} + a^{16} + 2R^8a^8 \cos{8\theta}} \\
    &\geq \sqrt{R^{16}+a^{16} - 2R^8a^8} \\
    &= |R^8-a^8| \\
    \therefore \int_{EA} &\leq i\int_0^{\pi} \dfrac{d\theta}{|R^8-a^8|} \\
    &= \dfrac{\pi i}{|R^8-a^8|}.
    \end{aligned}
\]

As $R\rightarrow \infty$, it is obvious therefore that $\displaystyle \int_{EA} \rightarrow 0$.

# Putting the pieces together 
Therefore the sum of our subintegrals (which is equal to our contour integral) is:

\[
    \begin{aligned}
    \oint_C \dfrac{e^{imz}}{z(z^8+a^8)} dz &= \int_{-\infty}^{\infty} \dfrac{e^{imx}}{x(x^8+a^8)} dx -\dfrac{\pi i}{a^8}
    \end{aligned}
\]

Replacing the left-hand side with what the residue theorem tells us this contour integral equals yields:

\[
    \begin{aligned}
    -\dfrac{\pi i}{2a^8} \left(e^{-amp}\cos{amq} + e^{-amq}\cos{amp}\right) &= \int_{-\infty}^{\infty} \dfrac{e^{imx}}{x(x^8+a^8)} dx -\dfrac{\pi i}{a^8} \\
    \therefore \int_{-\infty}^{\infty} \dfrac{e^{imx}}{x(x^8+a^8)} dx &= \dfrac{\pi i}{2a^8} \left(2 - e^{-amp}\cos{amq} - e^{-amq}\cos{amp}\right) \\
    \int_{-\infty}^{\infty} \dfrac{\cos{mx}+i\sin{mx}}{x(x^8+a^8)} dx &= \dfrac{\pi i}{2a^8} \left(2 - e^{-amp}\cos{amq} - e^{-amq}\cos{amp}\right) \\
    \int_{-\infty}^{\infty} \dfrac{\cos{mx}}{x(x^8+a^8)} dx + i\int_{-\infty}^{\infty} \dfrac{\sin{mx}}{x(x^8+a^8)} dx &= \dfrac{\pi i}{2a^8} \left(2 - e^{-amp}\cos{amq} - e^{-amq}\cos{amp}\right).
    \end{aligned}
\]

Equating the imaginary part of this expression yields:

\[
    \begin{aligned}
    \int_{-\infty}^{\infty} \dfrac{\sin{mx}}{x(x^8+a^8)} dx &= \dfrac{\pi}{2a^8} \left(2 - e^{-amp}\cos{amq} - e^{-amq}\cos{amp}\right) \\
    &= \dfrac{\pi}{2a^8} \left(2 - e^{-am\sin{\dfrac{\pi}{8}}}\cos{\left(am\cos{\dfrac{\pi}{8}}\right)} - e^{-am\cos{\dfrac{\pi}{8}}}\cos{\left(am\sin{\dfrac{\pi}{8}}\right)}\right).
    \end{aligned}
\]

Which is what we were required to determine.