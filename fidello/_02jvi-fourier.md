### Fourier transform
The Fourier transform is an integral transform that takes, as its argument, one function and returns another. It is of the form:
\begin{align*}
\hat{f}(k) &= \mathcal{F}(f) \\
&= \int_{-\infty}^{\infty} f(x) e^{-i k x} dx.
\end{align*}
Its inverse is:
\begin{align*}
f(x) &= \mathcal{F}^{-1}(\hat{f}) \\
&= \int_{-\infty}^{\infty} \hat{f}(k) e^{ikx} dk.
\end{align*}
We denote the Fourier transform of a specified function as $\mathcal{F}(f)$. It has broad applications, including in signal processing (as it identifies frequencies present in sounds) and in solving differential equations. This latter use is the result of the relationship between $\mathcal{F}(f^{(m)})$ and $f$ being a simple multiplication by $(ik)^m$.

It also has a discrete counterpart called the discrete Fourier transform. If you have a sequence of numbers corresponding to your values of $f(x)$ at $\{x_j\}=\{h, 2h, 3h, ..., 2\pi\}$ and $k\in \{-\frac{n}{2}+1, -\dfrac{n}{2}+2, -\dfrac{n}{2}+3, ..., \frac{n}{2}\}$ where $h=\dfrac{2\pi}{n}$ \cite{trefethen_spectral_2000}{18}, the discrete Fourier transform yields \cite{trefethen_spectral_2000}{18}:
\begin{align*}
\hat{f}_k &= h\sum_{j=1}^{n} f(x_j) e^{-ikx_j}.
\end{align*}
Its inverse is defined as \cite{trefethen_spectral_2000}{19}:
\begin{align*}
f(x_j) &= \dfrac{1}{2\pi} \sum_{k=-n/2+1}^{n/2} \hat{f}_k e^{ikx_j}.
\end{align*}
The fast Fourier transform is an especially efficient implementation of the discrete Fourier transform. This simply means that when you compute it, say in a programming language like Julia, it runs faster than you might otherwise expect while still returning the correct result. The DFT of the $m$th derivative of $f$ should be $(ik)^m\hat{f}_k$, except when $k=\frac{n}{2}$ and $m$ is odd in which case it would be $0$ \cite{trefethen_spectral_2000}{23}.