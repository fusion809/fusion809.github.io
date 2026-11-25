### Rearranging and integrating
For Equation \eqref{firstDE}, there would be a few ways we could solve this differential equation. As luckily each of the adjectives we used to describe it implies it should be easy to solve. One option is we could divide both sides by $y$ and integrate on $x$:
\begin{align*}
\int \dfrac{\dfrac{dy}{dx}}{y} dx &= \int dx \\
\int \dfrac{dy}{y} &= \int dx \\
\ln{y} &= x + C
\end{align*}
Exponentiating both sides yields:
\begin{align*}
y &= e^{x+C}.
\end{align*}
Note that we have one unknown constant, $C$. This is a product of how derivatives delete constants and is why we need initial, boundary or other conditions to uniquely determine our solution. For instance, if $y(0)=1$, then we can substitute this into our solution to find $C$.
\begin{align*}
1 &= e^{0+C} \\
e^C &= 1 \\
\implies C &= 0.
\end{align*}
Hence:
\begin{align*}
y &= e^x.
\end{align*}