### Complex analysis
\figure{contour}{right}{700px}{/assets/Infinite%20semicircle.svg}{An example of a curve we might calculate the contour integral of in complex analysis.}
I have definitely not covered everything within the field of complex numbers. There's a field of maths known as complex analysis where we combine calculus and complex numbers. In it, we learn, among other things, how to integrate along curves (or contours) in the complex plane and that there is an extension of Taylor series to the complex plane called the Laurent series. In it, we write:
\begin{align}
f(z) &= \sum_{n=-\infty}^{\infty} a_n (z-z_0)^n
\end{align}
where:
\begin{align}
a_n &= \dfrac{1}{2\pi i} \oint_C \dfrac{f(z)}{(z-z_0)^{n+1}} dz.
\end{align}
This is a contour integral over a curve $C$ that includes the point $z=z_0$. $z_0$ is analogous to the $a$ in our Taylor series as it is arbitrary. There's a formula called Cauchy's differentiation formula that says that:
\begin{align}
f^{(n)}(z_0) &= \dfrac{n!}{2\pi i} \oint_C \dfrac{f(z)}{(z-z_0)^{n+1}} dz, \hspace{0.1cm}n \in \{0, 1, 2, ...\}.
\end{align}
Where $f^{(n)}(z_0)$ is the $n$th derivative of $f(z)$ evaluated at $z=z_0$. Which helps us see how the Laurent series is connected to the Taylor series as for non-negative $n$, $a_n = \dfrac{f^{(n)}(z_0)}{n!}$. As for negative $n$, $a_n=0$ if $f(z)$ is defined everywhere in the complex plane. Why do I mention this? Well, I personally consider complex analysis one of the most beautiful fields of mathematics and I thought I would mention it in case you are ever in a position to take a complex analysis class while at university. 