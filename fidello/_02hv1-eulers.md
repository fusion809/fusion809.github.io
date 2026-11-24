#### Euler's method
Numerical solutions are actually quite interesting, because they typically involve applications of the Taylor series that we have covered. For instance, say we wanted to numerically integrate the differential equation $y'=y$ with initial condition $y(0)=1$. Then if we wanted to approximate $y(h)$ for some small $h$, we could use the Maclaurin series of $y$:
\begin{align*}
y(h) &= y(0) + hy'(0) + \dfrac{h^2}{2}y''(0) + \dfrac{h^3}{6}y'''(0) + ... \\
\end{align*}
Now our differential equation makes it very easy to approximate any arbitrary-order derivative of $y$ at $x=0$, as we can just keep differentiating our differential equation and find that $\dfrac{d^n y}{dx^n} = y$ for any $n\in \mathbb{Z}^{+}$. But say we did not have that luxury. Then what we could do is just approximate $y(h)$ using the first two terms of the Maclaurin series:
\begin{align*}
y(h) &\approx y(0) + hy'(0) \\
&= 1 + h \cdot 1 \\
&= 1 + h.
\end{align*}
We can then approximate $y(2h)$ using the same method:
\begin{align*}
y(2h) &\approx y(h) + hy'(h) \\
&\approx (1+h) + h(1+h) \\
&= 1 + 2h + h^2.
\end{align*}
And continue on from there. This method, where we truncate the Taylor series at the first two terms and numerically approximate our solution by applying this approximation at each step, is called Euler's method. As we ignore so many terms in the Taylor series, this method is rather inaccurate. Specifically, we call it a first-order method. 