### Laplace transform
The Laplace transform provides a method for converting linear ordinary differential equations into algebraic equations. It is defined as:
\begin{align*}
\mathcal{L}(f) &= \int_0^{\infty} f(t) e^{-st} dt.
\end{align*}
Where $s$ is a complex number. It is invertible, although the inverse Laplace transform requires some knowledge of complex analysis to define, so I am going to skip the definition. The Lapace transform of derivatives is also well-known:
\begin{align*}
\mathcal{L}(f^{n}) &= s^n \mathcal{L}(f) - \sum_{k=1}^{n} s^{n-k} f^{(k-1)}(0).
\end{align*}
$f^{(k-1)}(0-)=f^{(k-1)}(0)$ typically. If $f^{(k-1)}(0)$ is not clearly defined, we are required to determine $f^{(k-1)}(0-)$ through the limit:
\begin{align*}
f^{(k-1)}(0-) &= \lim_{x\rightarrow 0^{-}} f^{(k-1)}(x).
\end{align*}
$x\rightarrow 0^{-}$ means that we approach $0$ from the negative side of the number line. To see the power of this approach, let us consider the equation:
\begin{align*}
y'' - y &= \cos{x}.
\end{align*}
The Laplace transform of $y''$ is $s^2 \mathcal{L}(y) - sy(0) - y'(0)$. As for the $\cos{x}$ its Laplace transform is $\dfrac{s}{s^2+1}$. Hence our differential equation becomes
\begin{align*}
s^2 \mathcal{L}(y) - sy(0) - y'(0) - \mathcal{L}(y) &= \dfrac{s}{s^2+1} \\
(s^2-1) \mathcal{L}(y) &= sy(0) + y'(0) + \dfrac{s}{s^2+1} \\
\mathcal{L}(y) &= \dfrac{sy(0)+y'(0)}{s^2-1} + \dfrac{s}{(s^2-1)(s^2+1)}.
\end{align*}
Splitting the second right-hand side term through partial fraction decomposition:
\begin{align*}
\dfrac{s}{(s^2-1)(s^2+1)} &= \dfrac{As+B}{s^2-1} + \dfrac{Cs+D}{s^2+1} \\
s &= (As+B)(s^2+1) + (Cs+D)(s^2-1) \\
&= As^3+As + Bs^2 + B + Cs^3 - Cs + Ds^2 - D \\
&= (A+C)s^3 + (B+D)s^2 + (A-C)s + B-D.
\end{align*}
Equating coefficients yields:
\begin{align*}
A+C &= 0 \\
C &= -A \\
B+D &= 0 \\
D &= -B \\
A-C &= 1 \\
A-(-A) &= 1 \\
2A &= 1 \\
A &= \dfrac{1}{2} \\
C &= -\dfrac{1}{2} \\
B-D &= 0 \\
B-(-B) &= 0 \\
2B &= 0 \\
B &= 0 \\
\implies D &= 0.
\end{align*}
Hence:
\begin{align*}
\dfrac{s}{(s^2-1)(s^2+1)} &= \dfrac{1}{2} \dfrac{s}{s^2-1} - \dfrac{1}{2} \dfrac{s}{s^2+1}.
\end{align*}
Inverse Laplace transforming this yields:
\begin{align*}
\mathcal{L}^{-1}\left(\dfrac{s}{(s^2-1)(s^2+1)}\right) &= \dfrac{1}{2}\left(\cosh{x} - \cos{x}\right).
\end{align*}
As for the inverse Laplace transform of $\dfrac{sy(0)+y'(0)}{s^2-1}$ yields:
\begin{align*}
\mathcal{L}^{-1}\left(\dfrac{sy(0)+y'(0)}{s^2-1}\right) &= y(0) \cosh{x} + y'(0)\sinh{x}.
\end{align*}
Hence:
\begin{align*}
y(x) &= y(0) \cosh{x} + y'(0)\sinh{x} + \dfrac{1}{2}\left(\cosh{x} - \cos{x}\right) \\
&= \left(y(0)+\dfrac{1}{2}\right)\cosh{x} + y'(0)\sinh{x} - \dfrac{1}{2}\cos{x}.
\end{align*}