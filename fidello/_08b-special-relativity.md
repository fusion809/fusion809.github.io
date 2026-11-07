## Special relativity
As previously mentioned, in special relativity we can make relativistic corrections to kinetic energy in our Lagrangian and Hamiltonians. Namely, we can define it as:

\begin{align*}
T(v) &= \dfrac{mc^2}{\sqrt{1-\dfrac{v^2}{c^2}}} - mc^2 \\
&= mc^2\left(\dfrac{1}{\sqrt{1-\dfrac{v^2}{c^2}}}-1\right).
\end{align*}

### Classical limit of special relativity's kinetic energy
Using the Maclaurin series, we are going to show that $T(v)$ is well approximated by $T=\dfrac{mv^2}{2}$ when $v << c$. To make things simpler for ourselves, we will consider the Maclaurin series of the function:

\begin{align*}
f(x) &= \dfrac{1}{\sqrt{1-x}} \\
&= (1-x)^{-1/2}.
\end{align*}

(This is as opposed to directly calculating the Maclaurin series of our $T(v)$. If we directly applied the Maclaurin series formula to $T(v)$ the derivatives we would need to calculate would get a bit unruly.) In this case:
\begin{align*}
f'(x) &= -\dfrac{1}{2}(1-x)^{-3/2} \cdot -1 \\
&= \dfrac{(1-x)^{-3/2}}{2} \\
f''(x) &= -\dfrac{3(1-x)^{-5/2}}{4} \cdot -1 \\
&= \dfrac{3(1-x)^{-5/2}}{4} \\
f^{(3)}(x) &= -\dfrac{15(1-x)^{-7/2}}{8} \cdot -1 \\
&= \dfrac{15(1-x)^{-7/2}}{8}
\end{align*}
Where $f^{(3)}(x)$ denotes the 3rd derivative of $f$ with respect to $x$. Note how for the first derivative the coefficient is just $\dfrac{1}{2}$, for the second derivative it is $\dfrac{3}{4}$ and for the third it is $\dfrac{15}{8}$. Seems to me that the denominator is just $2^n$ where $n$ is the order of the derivative and the numerator of the coefficient is $(2n-1) \times (2n-3) \times (2n-5)\times ... \times 1$. Hence:
\begin{align*}
\dfrac{d^nf}{dx^n} &= \dfrac{\left(\displaystyle \prod_{i=1}^n 2i-1\right)(1-x)^{-(2n+1)/2}}{2^{n}} \hspace{0.1cm}\forall n\in \mathbb{Z}^{+}.
\end{align*}
$\displaystyle \prod_{i=1}^n 2i-1$ is a way of writing a multi-term product that is somewhat analogous to the summation notation $\sum_{i=1}^n 2i-1$. In this sum, we add up 2i-1 for $i=1$, $i=2$, $i=3$, and so forth, all the way up until $i=n$. Whereas in the product $\displaystyle \prod_{i=1}^n 2i-1$ we would multiply 2i-1 for $i=1$, $i=2$, $i=3$, and so forth up until $i=n$. Likewise $\mathbb{Z}^{+}$ denotes the set of all positive integers (so $\{1, 2, 3, ...\}$). For instance if $n=3$:
\begin{align*}
\displaystyle \prod_{i=1}^3 2i-1 &= (2(1)-1) \times (2(2)-1) \times (2(3)-1) \\
&= 1 \times 3 \times 5 \\
&= 15.
\end{align*}

Substituting in $x=0$ into $f(x)$ and our derivatives yields:
\begin{align*}
f(0) &= 1 \\
f'(0) &= \dfrac{1}{2} \\
\dfrac{d^nf}{dx^n}\Big|_{x=0} &= \dfrac{\displaystyle \prod_{i=1}^n 2i-1}{2^{n}}.
\end{align*}

In this case, $x=\dfrac{v^2}{c^2}$ and hence:
\begin{align*}
T(v) &= mc^2 \left(1 + \sum_{n=1}^{\infty} \dfrac{\left(\displaystyle \prod_{i=1}^n 2i-1\right)\left(\dfrac{v^2}{c^2}\right)^n}{2^{n}n!}-1\right)\\
&= mc^2 \sum_{n=1}^{\infty} \dfrac{\left(\displaystyle \prod_{i=1}^n 2i-1\right)\left(\dfrac{v^2}{c^2}\right)^n}{2^{n}n!}.
\end{align*}
Hence the first few terms of our Maclaurin series are:
\begin{align*}
T(v) &= mc^2 \left(\dfrac{v^2}{2c^2}+\dfrac{3v^4}{8c^4}+\dfrac{15v^6}{48c^6}+\dots\right) \\
&= \dfrac{mv^2}{2} + \dfrac{3mv^4}{8c^2} + \dfrac{15mv^6}{48c^4} + \dots.
\end{align*}
So what happens when $v$ is much smaller than $c$? Well, we would expect every term with a non-zero power of $c$ in the denominator to be negligible in its contribution to $T$ and this leaves only the first term. Hence we obtain that for $v << c$:
\begin{align*}
T(v) &\approx \dfrac{mv^2}{2}.
\end{align*}
