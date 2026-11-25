### Power series method
Another method to solve this and many other linear ODEs is to assume a power series expansion of $y$:
\begin{align*}
y(x) &= \sum_{n=0}^{\infty} a_n x^n \\
\therefore y'(x) &= \sum_{n=1}^{\infty} a_n nx^{n-1}.
\end{align*}
We ommitted $n=0$ from the $y'(x)$ power series because the derivative of $x^0=1$ is zero. We can rewrite the derivative power series as:
\begin{align*}
y'(x) &= \sum_{n=0}^{\infty} a_{n+1} (n+1) x^n.
\end{align*}
Hence our differential equation becomes:
\begin{align*}
\sum_{n=0}^{\infty} a_{n+1} (n+1) x^n &= \sum_{n=0}^{\infty} a_n x^n.
\end{align*}
Equating the coefficients of $x^n$ yields:
\begin{align*}
a_{n+1}(n+1) &= a_n.
\end{align*}
Dividing by $n+1$ yields:
\begin{align*}
a_{n+1} &= \dfrac{a_n}{n+1}
\end{align*}
Applying this recursive relationship over and over yields:
\begin{align*}
a_{n+1} &= \dfrac{a_{n-1}}{(n+1)n} \\
&= \dfrac{a_{n-2}}{(n+1)n(n-1)} \\
&= \dfrac{a_0}{(n+1)\times n\times(n-1)\times...\times1} \\
&= \dfrac{a_0}{(n+1)!}.
\end{align*}
Hence the coefficients of our power series are:
\begin{align*}
a_n &= \dfrac{a_0}{n!}.
\end{align*}
Then using our knowledge of Maclaurin series, we know that $a_0$ must be $y(0)$. $y(0)=1$ and hence:
\begin{align*}
a_n &= \dfrac{1}{n!}.
\end{align*}
Hence:
\begin{align*}
y &= \sum_{n=0}^{\infty} \dfrac{x^n}{n!}.
\end{align*}
Which we can recognize as the Maclaurin series of $e^x$. I guess we could have more directly used our knowledge of Maclaurin series and the fact that $y'=y$ implies that:
\begin{align*}
\dfrac{d^n y}{dx^n} = y.
\end{align*}
Hence $y(0)=1$ implies that:
\begin{align*}
\dfrac{d^n y}{dx^n} \Big|_{x=0} &= 1
\end{align*}
and hence:
\begin{align*}
y(x) &= \sum_{n=0}^{\infty} \dfrac{x^n}{n!}.
\end{align*}