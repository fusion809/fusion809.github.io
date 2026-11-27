##### Chebyshev spectral methods
Many different orthogonal functions can be used when applying spectral methods. Among the most popular are Chebyshev polynomials of the first kind and Fourier series functions. Chebyshev polynomials of the first kind, $T_i(x)$, are defined on $x\in [-1,1]$ as:
\begin{align*}
T_j(x) = \cos(j\cos^{-1}{x}).
\end{align*}
I would imagine that you are looking at this expression and wondering, how can $T_j$ be a polynomial when it is defined in terms of cosine and its inverse? Say $j=0$, we get:
\begin{align*}
T_0 &= \cos{0} \\
&=1.
\end{align*} 
If $j=1$, we get:
\begin{align*}
T_1 &=\cos{\left(\cos^{-1}{x}\right)} \\
&=x.
\end{align*} 
If $j_2$, we get:
\begin{align*}
T_2 &= \cos{\left(2\cos^{-1}{x}\right)} \\
&= \cos^2{\left(\cos^{-1}{x}\right)}-\sin^2{\left(\cos^{-1}{x}\right)} \\
&= \cos^2{\left(\cos^{-1}{x}\right)}- (1-\cos^2{\left(\cos^{-1}{x}\right)}) \\
&= x^2 - (1-x^2) \\
&= 2x^2 - 1.
\end{align*}
This pattern continues on from here, too. Further, this is just one way to define Chebyshev polynomials of the first kind and is the most convenient way for when you are using them to approximate the solution to differential equations. That being said, they can also be defined using a recursion relation that makes the fact that they are polynomials all the more clear. Specifically this is the recursion relation \cite{boyd_chebyshev_2000}{323}:
\begin{align*}
T_0 &= 1, \hspace{0.2cm} T_1 &&= x \\
T_{j+1}(x) &= 2xT_j(x) - T_{j-1}(x) && \forall j \geq 1.
\end{align*}

The grid points normally used are:
\begin{align*}
x_i &= -\cos{\left(\dfrac{\pi i}{n}\right)}.
\end{align*}
These points extremize the $n$th Chebyshev polynomial of the first kind and this makes the matrices we get easier to work with, as these matrices will be better conditioned (that is, have a condition number close to one). 

I should mention that while our basis functions are limited to $x\in [-1,1]$, this does not mean we are limited to $x\in[-1,1]$ for the differential equations we can solve. As we can apply any number of different transformations to transform $x_i$ from $[-1,1]$ to other domains. The most straight-forward is a linear transformation. Say we want to transform to $y\in [a,b]$ we would use:
\begin{align*}
y &= \dfrac{b-a}{2}x + \dfrac{b+a}{2}.
\end{align*}

We can also define a differentiation matrix. If $\mathbf{dT} = (T'_j(x_i))$ and $\mathbf{T}=(T_j(x_i))$ we can define the differentiation matrix:
\begin{align*}
\mathbf{D} &= \mathbf{dT} \mathbf{T}^{-1}.
\end{align*}
This matrix, when run on a vector of function values at $x_i$ computes its first derivative. If you multiply the vector by $\mathbf{D}$ twice, you will get the second derivative. Let $\mathbf{D}^{(2)}$ represent the second derivative matrix; it will equal $\mathbf{D}^2=\mathbf{D}\mathbf{D}$. It can also be calculated using $\mathbf{d^2T}=(T''_j(x_i))$ as $\mathbf{D}^{(2)}=\mathbf{d^2T}\mathbf{T}^{-1}$. Luckily, there are also more efficient ways of calculating $\mathbf{D}$, such as using \cite{trefethen_spectral_2000}{53}:
\begin{align*}
D_{ij} = \begin{cases}
    -\dfrac{2n^2+1}{6}, & i=j=0 \\\\
    \dfrac{2n^2+1}{6}, & i=j=n \\\\
    -\dfrac{x_j}{2(1-x_j^2)}, & i=j, i\neq 0, n \\\\
    \dfrac{(1+\delta_{i0}+\delta_{in})(-1)^{i+j}}{(1+\delta_{j0}+\delta_{jn})(x_i-x_j)}, & i\neq j.
\end{cases}
\end{align*}