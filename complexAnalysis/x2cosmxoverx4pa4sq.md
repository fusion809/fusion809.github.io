@def title="Evaluating example definite integral on the infinite domain using contour integration"
@def tags="Contour integration"

In this article, we shall evaluate the following integral using contour integration methods.

\[
    I = \int_{-\infty}^{\infty} \dfrac{x^2\cos{mx}}{(x^4+a^4)^2} dx
\]

\toc

## Contour and integrand hunting
Our first couple of questions, which are impossible to separate from one another as they are closely connected, is which curve and which integrand should we use to determine the value of $I$? The contour and integrand should satisfy the following criteria:

1. The contour must have a component whose line integral evaluates to either $I$, or $\displaystyle \int_{0}^{\infty} \dfrac{x^2\cos{mx}}{(x^4+a^4)^2} dx = \dfrac{1}{2} I$ (thanks to the symmetry of our integrand). It is also OK if $I$ or half $I$ is only stored in the real or imaginary part of the integral along this segment, and there's another value in the other part of the integral for this segment.
2. All other components of the contour must either only contribute to a separate part of the result (e.g. if $I$ is stored in the real part of the solution, the other parts of the integral may only contribute to the imaginary part of the solution), or equal 0.
3. There must be no singularities of our integrand along this contour.
4. The only singularities allowed within the region enclosed by the contour are poles.

One curve and integrand that would satisfy these criteria is the infinite semi-circle that covers the entire positive $y$-axis, but does not cover any of the negative $y$-axis and is centred on the origin. Namely, this curve:

~~~
<img src="/assets/Infinite%20semicircle.svg" width="60%">
~~~

with this contour integral:

\[
    \oint_C \dfrac{z^2 e^{imz}}{(z^4+a^4)^2} dz.
\]

## Evaluating the contour integral
The residue theorem tells us that our contour integral will equal $2\pi i \sum \mathrm{residues}$. These residues will correspond to when the denominator of our integrand equals zero, that is when:

\[
\begin{aligned}
    z^4+a^4 &= 0 \\
    z^4 &= -a^4 \\
    z^4 &= a^4 \cdot (-1) \\
    z^4_n &= a^4 e^{\pi i + 2\pi i n}, \hspace{0.1cm} n\in \mathbb{Z} \\
    z_n &= a e^{\dfrac{\pi i(1+2n)}{4}}.
\end{aligned}
\]

$n=0, 1, 2, 3$ gives us the only unique roots of the denominator and hence our only poles. It is also important to note that $z^4+a^4 = (z-z_0)(z-z_1)(z-z_2)(z-z_3)$. Of these poles, only $z_0$ ($z_n$ for $n=0$) and $z_1$ lie within the region enclosed by C. Both these poles will be second order, as to remove the singularity, our integrand will need to be multiplied by $(z-z_n)^2$. 

The general formula for calculating residues is (where $p$ is a pole of order $k$ in the function, $f(z)$):

\[
    \begin{aligned}
    \mathrm{residue} = \lim_{z\rightarrow p} \dfrac{1}{(k-1)!} \dfrac{d^{k-1}}{dz^{k-1}} \left((z-p)^k f(z)\right)
    \end{aligned}
\]

For $z_0$, the residue $r_0$ is:
\[
    \begin{aligned}
    r_0 &= \lim_{z\rightarrow z_0} \dfrac{1}{(2-1)!} \dfrac{d^{2-1}}{dz^{2-1}} \left((z-z_0)^2 \dfrac{z^2e^{imz}}{(z^4+a^4)^2}\right) \\
    &= \lim_{z\rightarrow z_0} \dfrac{d}{dz} \left(\dfrac{z^2 e^{imz}}{(z-z_1)^2(z-z_2)^2(z-z_3)^2}\right) \\
    &= \lim_{z\rightarrow z_0} \dfrac{2ze^{imz}}{(z-z_1)^2(z-z_2)^2(z-z_3)^2} + \dfrac{imz^2 e^{imz}}{(z-z_1)^2(z-z_2)^2(z-z_3)^2} - \dfrac{2z^2e^{imz}}{(z-z_1)^2(z-z_2)^2(z-z_3)^2} \\\left(\dfrac{1}{z-z_1}+\dfrac{1}{z-z_2}+\dfrac{1}{z-z_3}\right) \\
    &= \lim_{z\rightarrow z_0} \dfrac{ze^{imz}}{(z-z_1)^2(z-z_2)^2(z-z_3)^2} \left(2+imz-2z\left(\dfrac{1}{z-z_1}+\dfrac{1}{z-z_2}+\dfrac{1}{z-z_3}\right)\right) \\
    &= \dfrac{z_0 e^{imz_0}}{(z_0-z_1)^2(z_0-z_2)^2(z_0-z_3)^2} \left(2 + imz_0 - 2z_0 \left(\dfrac{1}{z_0-z_1}+\dfrac{1}{z_0-z_2}+\dfrac{1}{z_0-z_3}\right)\right)
    \end{aligned}
\]

To simplify our above equation, we need to simplify our denominators:

\[
    \begin{aligned}
    z_0-z_1 &= a \left(e^{\dfrac{\pi i}{4}} - e^{\dfrac{3\pi i}{4}}\right) \\
    &= ae^{\dfrac{\pi i}{4}} \left(1-e^{\dfrac{\pi i}{2}}\right) \\
    &= ae^{\dfrac{\pi i}{4}} (1-i) \\
    &= ae^{\dfrac{\pi i}{4}} \sqrt{2} e^{-\dfrac{\pi i}{4}} \\
    &= a\sqrt{2}.\\
    z_0 - z_2 &= a \left(e^{\dfrac{\pi i}{4}} - e^{\dfrac{5\pi i}{4}}\right) \\
    &= ae^{\dfrac{\pi i}{4}} \left(1-e^{\pi i}\right) \\
    &= ae^{\dfrac{\pi i}{4}} (1 - (-1)) \\
    &= 2ae^{\dfrac{\pi i}{4}}.\\
    z_0 - z_3 &= a \left(e^{\dfrac{\pi i}{4}} - e^{\dfrac{7\pi i}{4}}\right) \\
    &= ae^{\dfrac{\pi i}{4}} \left(1-e^{\dfrac{3\pi i}{2}}\right) \\
    &= ae^{\dfrac{\pi i}{4}} (1 - (-i))\\
    &= ae^{\dfrac{\pi i}{4}} (1+i) \\
    &= ae^{\dfrac{\pi i}{4}} \sqrt{2} e^{\dfrac{\pi i}{4}} \\
    &= a\sqrt{2} e^{\dfrac{\pi i}{2}} \\
    &= a\sqrt{2} i.
    \end{aligned}
\]

Likewise, $iz_0 = e^{\dfrac{\pi i}{2}} \cdot ae^{\dfrac{\pi i}{4}} = ae^{\dfrac{3\pi i}{4}} = \dfrac{a}{\sqrt{2}} (-1+i)$.

\[
    \begin{aligned}
    r_0 &= \dfrac{ae^{\dfrac{\pi i}{4}} e^{\dfrac{am}{\sqrt{2}}(-1+i)}}{2a^2 \cdot 4a^2e^{\dfrac{\pi i}{2}} \cdot 2a^2 e^{\pi i}} \left(2+\dfrac{am}{\sqrt{2}}(-1+i) - 2ae^{\dfrac{\pi i}{4}}\left(\dfrac{1}{a\sqrt{2}} + \dfrac{1}{2ae^{\dfrac{\pi i}{4}}} + \dfrac{1}{a\sqrt{2}e^{\dfrac{\pi i}{2}}}\right)\right) \\
    &= \dfrac{e^{-\dfrac{5\pi i}{4}} e^{\dfrac{am}{\sqrt{2}}(-1+i)}}{16a^5} \left(2+\dfrac{am}{\sqrt{2}}(-1+i) - \sqrt{2}e^{\dfrac{\pi i}{4}} - 1 - \sqrt{2}e^{-\dfrac{\pi i}{4}}\right) \\
    \end{aligned}
\]

Using: $e^{\dfrac{\pi i}{4}} + e^{-\dfrac{\pi i}{4}} = \sqrt{2}$ and $e^{-\dfrac{5\pi i}{4}} = e^{\dfrac{3\pi i}{4}}$, we can rewrite the above expression as:

\[
    \begin{aligned}
    r_0 &= \dfrac{e^{\dfrac{3\pi i}{4}} e^{\dfrac{am}{\sqrt{2}}(-1+i)}}{16a^5} \left(\dfrac{am}{\sqrt{2}}(-1+i) -1\right) \\
    &= \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1+i)}}{16a^5\sqrt{2}}\left(\dfrac{am}{\sqrt{2}}(-1+i)^2 - (-1+i)\right) \\
    &= \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1+i)}}{16a^5\sqrt{2}}\left(\dfrac{am}{\sqrt{2}}\cdot -2i+1-i\right) \\
    &= \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1+i)}}{16a^5\sqrt{2}}\left(-i\left(am\sqrt{2}+1\right)+1\right)
    \end{aligned}
\]

Next, we must find $r_1$. 

\[
    \begin{aligned}
    r_1 &= \lim_{z\rightarrow z_1} \dfrac{d}{dz} \left((z-z_1)^2 \dfrac{z^2e^{imz}}{(z^4+a^4)^2}\right) \\
    &= \lim_{z\rightarrow z_1} \dfrac{d}{dz} \left(\dfrac{z^2 e^{imz}}{(z-z_0)^2(z-z_2)^2(z-z_3)^2}\right) \\
    &= \lim_{z\rightarrow z_1} \dfrac{2ze^{imz}}{(z-z_0)^2(z-z_2)^2(z-z_3)^2} + \dfrac{imz^2 e^{imz}}{(z-z_0)^2(z-z_2)^2(z-z_3)^2} - \dfrac{2z^2e^{imz}}{(z-z_0)^2(z-z_2)^2(z-z_3)^2} \\\left(\dfrac{1}{z-z_0}+\dfrac{1}{z-z_2}+\dfrac{1}{z-z_3}\right) \\
    &= \lim_{z\rightarrow z_1} \dfrac{ze^{imz}}{(z-z_0)^2(z-z_2)^2(z-z_3)^2} \left(2+imz-2z\left(\dfrac{1}{z-z_0}+\dfrac{1}{z-z_2}+\dfrac{1}{z-z_3}\right)\right) \\
    &= \dfrac{z_1 e^{imz_1}}{(z_1-z_0)^2(z_1-z_2)^2(z_1-z_3)^2} \left(2 + imz_1 - 2z_1 \left(\dfrac{1}{z_1-z_0}+\dfrac{1}{z_1-z_2}+\dfrac{1}{z_1-z_3}\right)\right).
    \end{aligned}
\]

First, we shall find $z_1-z_0$, $z_1-z_2$ and $z_1-z_3$ as these formulae will simplify the equations we have:

\[
    \begin{aligned}
    z_1 - z_0 &= -(z_0-z_1) \\
    &= -a\sqrt{2}.\\
    z_1 - z_2 &= a \left(e^{\dfrac{3\pi i}{4}}-e^{\dfrac{5\pi i}{4}}\right) \\
    &= ae^{\dfrac{3\pi i}{4}}(1-e^{\dfrac{\pi i}{2}}) \\
    &= ae^{\dfrac{3\pi i}{4}}(1-i) \\
    &= ae^{\dfrac{3\pi i}{4}} \cdot \sqrt{2} e^{-\dfrac{\pi i}{4}} \\
    &= a\sqrt{2} e^{\dfrac{\pi i}{2}}.\\
    z_1 - z_3 &= a \left(e^{\dfrac{3\pi i}{4}}-e^{\dfrac{7\pi i}{4}}\right) \\
    &= ae^{\dfrac{3\pi i}{4}} (1-e^{\pi i}) \\
    &= ae^{\dfrac{3\pi i}{4}} (1 - (-1)) \\
    &= 2a e^{\dfrac{3\pi i}{4}}. \\
    iz_1 &= e^{\dfrac{\pi i}{2}} ae^{\dfrac{3\pi i}{4}} \\
    &= ae^{\dfrac{5\pi i}{4}} \\
    &= \dfrac{a}{\sqrt{2}} (-1-i).
    \end{aligned}
\]

\[
    \begin{aligned}
    r_1 &= \dfrac{ae^{\dfrac{3\pi i}{4}}e^{\dfrac{am}{\sqrt{2}}(-1-i)}}{2a^2 \cdot 2a^2 e^{\pi i} \cdot 4a^2 e^{\dfrac{3\pi i}{2}}} \left(2 + \dfrac{am}{\sqrt{2}}(-1-i)-2ae^{\dfrac{3\pi i}{4}}\left(-\dfrac{1}{a\sqrt{2}}+\dfrac{1}{a\sqrt{2}e^{\dfrac{\pi i}{2}}}+\dfrac{1}{2ae^{\dfrac{3\pi i}{4}}}\right)\right) \\
    &= \dfrac{e^{-\dfrac{7\pi i}{4}e^{\frac{am}{\sqrt{2}}(-1-i)}}}{16a^5} (2-\dfrac{am}{\sqrt{2}}(1+i)+\sqrt{2}e^{\dfrac{3\pi i}{4}} - \sqrt{2} e^{\dfrac{\pi i}{4}}-1)
    \end{aligned}
\]

To simplify the above expression, let us do some complex arithmetic:

\[
    \begin{aligned}
    \sqrt{2} e^{\dfrac{3\pi i}{4}} - \sqrt{2}e^{\dfrac{\pi i}{4}} &= \sqrt{2} \cdot \dfrac{1}{\sqrt{2}} (-1+i - 1 - i) \\
    &= -2.\\
    e^{-\dfrac{7\pi i}{4}} &= e^{\dfrac{\pi i}{4}} \\
    &= \dfrac{1}{\sqrt{2}}(1+i).
    \end{aligned}
\]

\[
    \begin{aligned}
    r_1 &= \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1-i)}(1+i)}{16a^5 \sqrt{2}}\left(-\dfrac{am}{\sqrt{2}}(1+i)-1\right) \\
    &= \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1-i)}}{16a^5 \sqrt{2}}\left(-am\sqrt{2}i - 1 - i\right) \\
    &= \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1-i)}}{16a^5 \sqrt{2}}\left(-i(1+am\sqrt{2})-1\right)
    \end{aligned}
\]

The sum of the residues is thus:

\[
    \begin{aligned}
    \sum \mathrm{residues} &= \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1+i)}}{16a^5\sqrt{2}}\left(-i\left(am\sqrt{2}+1\right)+1\right) + \dfrac{e^{\dfrac{am}{\sqrt{2}}(-1-i)}}{16a^5 \sqrt{2}}\left(-i(1+am\sqrt{2})-1\right)\\
    &= \dfrac{e^{-\dfrac{am}{\sqrt{2}}}}{16a^5\sqrt{2}} \left(e^{\dfrac{ami}{\sqrt{2}}}-e^{-\dfrac{ami}{\sqrt{2}}} -i (am\sqrt{2}+1)\left(e^{\dfrac{ami}{\sqrt{2}}}+e^{-\dfrac{ami}{\sqrt{2}}}\right)\right) \\
    &= \dfrac{e^{-\dfrac{am}{\sqrt{2}}}}{16a^5\sqrt{2}} \left(2i\sin{\dfrac{am}{\sqrt{2}}} -2i (am\sqrt{2}+1)\cos{\dfrac{am}{\sqrt{2}}}\right) \\
    &= \dfrac{e^{-\dfrac{am}{\sqrt{2}}}i}{8a^5\sqrt{2}}\left(\sin{\dfrac{am}{\sqrt{2}}} -(am\sqrt{2}+1)\cos{\dfrac{am}{\sqrt{2}}}\right) \\
    \end{aligned}
\]

Therefore our contour integral equals:

\[
    \begin{aligned}
    \oint_C \dfrac{z^2 e^{imz}}{(z^4+a^4)^2} dz &= 2 \pi i \sum \mathrm{residues} \\
    &= \dfrac{e^{-\dfrac{am}{\sqrt{2}}}\pi}{4a^5\sqrt{2}}\left((am\sqrt{2}+1)\cos{\dfrac{am}{\sqrt{2}}}-\sin{\dfrac{am}{\sqrt{2}}}\right).\hspace{5cm}\tag{1}
    \end{aligned}
\]

## Evaluating integral components
Breaking our contour integral into its components:

\[
    \begin{aligned}
    \oint_C \dfrac{z^2 e^{imz}}{(z^4+a^4)^2} dz &= \int_{AB} \dfrac{z^2 e^{imz}}{(z^4+a^4)^2} dz + \int_{BA} \dfrac{z^2 e^{imz}}{(z^4+a^4)^2} dz
    \end{aligned}
\]

Along AB: $z=x$, $dz=dx$ and $x$ goes from $-R$ to $R$.
Along BA: $z=Re^{i\theta}$, $dz=iRe^{i\theta}d\theta$ and $\theta$ goes from $0$ to $\pi$.

\[
    \begin{aligned}
    \int_{AB} \dfrac{z^2 e^{imz}}{(z^4+a^4)^2} dz + \int_{BA} \dfrac{z^2 e^{imz}}{(z^4+a^4)^2} dz &= \int_{-R}^R \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + \int_{0}^{\pi} \dfrac{R^2 e^{2i\theta} e^{imRe^{i\theta}}}{(R^4e^{4i\theta}+a^4)^2} iRe^{i\theta} d\theta \\
    &= \int_{-R}^R \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + i\int_{0}^{\pi} \dfrac{R^3 e^{3i\theta} e^{imRe^{i\theta}}}{(R^4e^{4i\theta}+a^4)^2} d\theta. \\
    \end{aligned}
\]

Taking the limit as $R\rightarrow\infty$:

\[
\begin{aligned}
\lim_{R\rightarrow \infty} \int_{-R}^R \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + i\int_{0}^{\pi} \dfrac{R^3 e^{3i\theta} e^{imRe^{i\theta}}}{(R^4e^{4i\theta}+a^4)^2} d\theta &= \int_{-\infty}^{\infty} \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + \lim_{R\rightarrow \infty} i\int_{0}^{\pi} \dfrac{R^3 e^{3i\theta} e^{imRe^{i\theta}}}{(R^4e^{4i\theta}+a^4)^2} d\theta \\
&= \int_{-\infty}^{\infty} \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + \lim_{R\rightarrow \infty} i\int_{0}^{\pi} \dfrac{R^3 e^{3i\theta} e^{imRe^{i\theta}}}{R^8e^{8i\theta}+a^8+2a^4R^4e^{4i\theta}} d\theta \\
&= \int_{-\infty}^{\infty} \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + \lim_{R\rightarrow \infty} i\int_{0}^{\pi} \dfrac{e^{imRe^{i\theta}}}{R^5e^{5i\theta}+a^8e^{-3i\theta}R^{-3}+2a^4Re^{i\theta}} d\theta \\
&= \int_{-\infty}^{\infty} \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + \lim_{R\rightarrow \infty} i\int_{0}^{\pi} \dfrac{e^{imR\cos{\theta}-mR\sin{\theta}}}{R^5e^{5i\theta}+a^8e^{-3i\theta}R^{-3}+2a^4Re^{i\theta}} d\theta \\
&= \int_{-\infty}^{\infty} \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx + \lim_{R\rightarrow \infty} i\int_{0}^{\pi} \dfrac{e^{imR\cos{\theta}}e^{-mR\sin{\theta}}}{R^5e^{5i\theta}+a^8e^{-3i\theta}R^{-3}+2a^4Re^{i\theta}} d\theta
\end{aligned}
\]

Within the domain of integration for the second integral, $\sin{\theta} \geq 0$ and therefore $e^{-mR\sin{\theta}} \leq 1$ (assuming $m\geq0$). As the denominator contains some positive powers of $R$, as $R\rightarrow\infty$ the denominator goes to infinity and therefore the integrand goes to zero. The integral of 0 is 0 and therefore our second integral becomes zero. 

## Conclusion
Therefore our contour integral becomes:

\[
    \begin{aligned}
    \int_{-\infty}^{\infty} \dfrac{x^2 e^{imx}}{(x^4+a^4)^2} dx &= \int_{-\infty}^{\infty} \dfrac{x^2 (\cos{mx}+i\sin{mx})}{(x^4+a^4)^2} dx \\
    &= \int_{-\infty}^{\infty} \dfrac{x^2 \cos{mx}}{(x^4+a^4)^2} dx + i \int_{-\infty}^{\infty} \dfrac{x^2 \sin{mx}}{(x^4+a^4)^2} dx
    \end{aligned}
\]

Equating with the right-hand side of (1):

\[
    \begin{aligned}
    \int_{-\infty}^{\infty} \dfrac{x^2 \cos{mx}}{(x^4+a^4)^2} dx + i \int_{-\infty}^{\infty} \dfrac{x^2 \sin{mx}}{(x^4+a^4)^2} dx &= \dfrac{e^{-\dfrac{am}{\sqrt{2}}}\pi}{4a^5\sqrt{2}}\left((am\sqrt{2}+1)\cos{\dfrac{am}{\sqrt{2}}}-\sin{\dfrac{am}{\sqrt{2}}}\right)
    \end{aligned}
\]

Equating real and imaginary parts:

\[
    \begin{aligned}
    \int_{-\infty}^{\infty} \dfrac{x^2 \cos{mx}}{(x^4+a^4)^2} dx &= \dfrac{e^{-\dfrac{am}{\sqrt{2}}}\pi}{4a^5\sqrt{2}}\left((am\sqrt{2}+1)\cos{\dfrac{am}{\sqrt{2}}}-\sin{\dfrac{am}{\sqrt{2}}}\right) \\
    \int_{-\infty}^{\infty} \dfrac{x^2 \sin{mx}}{(x^4+a^4)^2} dx &= 0.
    \end{aligned}
\]

It is important to note, however, that because $\cos{(-mx)}=\cos{mx}$ and $(-a)^4=a^4$, the sign of $a$ and $m$ should not matter, therefore provided $a$ and $m$ are real:

\[
\int_{-\infty}^{\infty} \dfrac{x^2 \cos{mx}}{(x^4+a^4)^2} dx = \dfrac{e^{-\dfrac{|a||m|}{\sqrt{2}}}\pi}{4|a|^5\sqrt{2}}\left((|a||m|\sqrt{2}+1)\cos{\dfrac{am}{\sqrt{2}}}-\sin{\dfrac{|a||m|}{\sqrt{2}}}\right)
\]