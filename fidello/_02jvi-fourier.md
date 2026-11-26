### Fourier transform
The Fourier transform is an integral transform that takes, as its argument, one function and returns another. It is of the form:
\begin{align*}
\hat{f}(\varphi) &= \int_{-\infty}^{\infty} f(x) e^{-2\pi i \varphi x} dx.
\end{align*}
Its inverse is:
\begin{align*}
f(x) &= \int_{-\infty}^{\infty} \hat{f}(\varphi) e^{2\pi i \varphi x} d\varphi.
\end{align*}
We denote the Fourier transform of a specified function as $\mathcal{F}(f)$. It has broad applications, including in signal processing (as it identifies frequencies present in sounds) and in solving differential equations. This latter use is the result of the relationship between $\mathcal{F}(f^{(n)})$ and $f$ being a simple multiplication by $i^n \varphi^n$.

It also has a discrete counterpart called the discrete Fourier transform. If you have a sequence of numbers $\{x_j\}=\{x_0, x_1, x_2, ..., x_{n-1}\}$, the discrete Fourier transform yields:
\begin{align*}
X_j &= \sum_{k=0}^{n-1} x_k e^{-\frac{2\pi i jk}{n}}.
\end{align*}
Its inverse is defined as:
\begin{align*}
x_k &= \dfrac{1}{n} \sum_{j=0}^{n-1} X_j e^{\frac{2\pi i jk}{n}}.
\end{align*}
The fast Fourier transform is an especially efficient implementation of the discrete Fourier transform. 