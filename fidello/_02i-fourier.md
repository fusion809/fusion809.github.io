## Fourier series
Fourier series are a way of rewriting periodic functions in terms of sine and cosine waves. Specifically, if we wish to write $f(x)$ as a Fourier series over the domain $x \in [-L,L]$, we would write (keep in mind the series coefficients $a_n$ and $b_n$ are not the same as those we have been using so far):
\begin{align*}
f(x) &= a_0 + \sum_{n=1}^{\infty} a_n \cos{\left(\dfrac{n\pi}{L}x\right)} + b_n \sin{\left(\dfrac{n\pi}{L}x\right)}. 
\end{align*}
Where:
\begin{align}
a_0 &= \dfrac{1}{2L} \int_{-L}^{L} f(x) dx \label{Fouriera0} \\
a_n &= \dfrac{1}{L} \int_{-L}^{L} f(x) \cos{\left(\dfrac{2\pi nx}{L}\right)} dx \label{Fourieran}\\
b_n &= \dfrac{1}{L} \int_{-L}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx \label{Fourierbn}.
\end{align}