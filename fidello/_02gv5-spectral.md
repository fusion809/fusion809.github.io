#### Spectral methods
Sturm-Liouville problems also give rise to infinite sets of orthogonal functions that are useful when numerically approximating the solution to other differential equations.[^2] If $\{\phi_i\}$, $i\in \{0,1,2,...\}$ is our orthogonal function set, we assume our solution, $u$, to some ordinary differential equation can be written as:
\begin{align}
u(x) &= \sum_{j=0}^{\infty} a_j \phi_j(x). \label{spectral} 
\end{align}
Then we discretize our grid into points $x_i$, truncate the above series expansion, and can write a matrix equation for $\mathbf{a}$ (a vector containing $\{a_j\}$, the series coefficients). Then this equation would be:
\begin{align}
\mathbf{A}\mathbf{a} &= \mathbf{b}. \label{speca}
\end{align}
Let $\boldsymbol{\phi}$ be a matrix containing $\phi_j(x_i)$ and $\mathbf{u}$ be a vector of $u(x_i)$ values. Then we can write:
\begin{align}
\mathbf{a} &=\boldsymbol{\phi}^{-1} \mathbf{u}.\label{afromu}
\end{align}
We can substitute Equation \eqref{afromu} into Equation \eqref{speca}, yielding:
\begin{align*}
\mathbf{A}\boldsymbol{\phi}^{-1}\mathbf{u} &= \mathbf{b}.
\end{align*}
In this case, we would be solving for $\mathbf{u}$ instead of $\mathbf{a}$. Methods that approximate the solution to a differential equation over the entire domain of the problem using a single sum of basis functions like this are called spectral methods. 

If the domain is split up into smaller sections with a separate function of the form Equation \eqref{spectral} defined on each, this is called a spectral element method. Finite element methods are very similar as they do split the domain up and define a series-based solution on each part of the domain. The difference is that the series usually has fewer terms and achieves lower accuracy than their spectral element counterparts.

Many different orthogonal functions can be used when applying spectral methods. Among the most popular are Chebyshev polynomials of the first kind and Fourier series functions. Chebyshev polynomials of the first kind, $T_i(x)$, are defined on $x\in [-1,1]$ as:
\begin{align*}
T_j(x) = \cos(j\cos^{-1}{x}).
\end{align*}
I would imagine that you are looking at this expression and wondering, how can $T_j$ be a polynomial when it is defined in terms of cosine and its inverse? Well, this is just one way to define Chebyshev polynomials of the first kind and is the most convenient way for when you are using them to approximate the solution to differential equations. That being said, they can also be defined using a recursion relation that makes the fact that they are polynomials all the more clear. Specifically this is the recursion relation:[^2]
\begin{align*}
T_0 &= 1, \hspace{0.2cm} T_1 &&= x \\
T_{j+1}(x) &= 2xT_j(x) - T_{j-1}(x) && \forall j \geq 1.
\end{align*}

The grid points normally used are:
\begin{align*}
x_i &= \cos{\dfrac{\pi i}{n}}.
\end{align*}
These points extremize the $n$th Chebyshev polynomial of the first kind and this makes the matrices we get easier to work with. 

I should mention that while our basis functions are limited to $x\in [-1,1]$, this does not mean we are limited to $x\in[-1,1]$ for the differential equations we can solve. As we can apply any number of different transformations to transform $x_i$ from $[-1,1]$ to other domains. 

We can also define a differentiation matrix. If $\mathbf{dT} = (T'_j(x_i))$ and $\mathbf{T}=(T_j(x_i))$ we can define the differentiation matrix:
\begin{align*}
\mathbf{D} &= \mathbf{dT} \mathbf{T}^{-1}.
\end{align*}
This matrix, when run on a vector of function values at $x_i$ computes its first derivative. If you multiply the vector by $\mathbf{D}$ twice, you will get the second derivative. Luckily, there are also more efficient ways of calculating $\mathbf{D}$. Like with the finite difference method, spectral methods are also easiest applied to linear differential equations. Interestingly, spectral methods are of order that is proportional to the number of grid points used. Which means that while increasing the number of grid points does improve the accuracy of every method one can use for solving differential equations, spectral methods take this to a whole other level. They are often orders of magnitude more accurate than other methods at a given $n$. 

Say we wish to solve $u''-xu=0$ on $x\in[0,\infty]$. The exact solution to this differential equation is:
\begin{align*}
u(x) &= C_1 \mathrm{Ai}(x) + C_2 \mathrm{Bi}(x).
\end{align*}
Where $\mathrm{Ai}(x)$ and $\mathrm{Bi}(x)$ are Airy functions, which are actually defined as the solution to this exact ODE. Say we have the boundary conditions $u(0)=\mathrm{Ai}(0)$ and $u(\infty)=0$. Then the solution would be $u(x)=\mathrm{Ai}(x)$. This is because:
\begin{align*}
\lim_{x\rightarrow \infty} \mathrm{Bi}(x) = \infty
\end{align*}
and:
\begin{align*}
\lim_{x\rightarrow \infty} \mathrm{Ai}(x)=0.
\end{align*}
Then to transform our grid points from $[-1,1]$ to $[0,\infty]$ we will use the transformation (where $y$ are our Chebyshev grid points and $x$ is our semi-infinite domain grid points):
\begin{align*}
x &= L \dfrac{1+y}{1-y}.
\end{align*}
Where $L$ is a parameter we might need to fine tune to get the best approximation. Hence:
\begin{align*}
\dfrac{d}{dx} &= \dfrac{1}{x'(y)} \dfrac{d}{dy} \\
&= \dfrac{1}{\dfrac{(1-y)+(1+y)}{(1-y)^2}} \dfrac{d}{dy} \\
&= \dfrac{1}{\dfrac{2L}{(1-y)^2}} \dfrac{d}{dy} \\
&= \dfrac{(1-y)^2}{2L} \dfrac{d}{dy}.
\end{align*}
Which allows us to calculate our differentiation matrix for $x$. It will simply be:
\begin{align*}
\mathbf{D}_x &= \mathbf{Diag}\left(\dfrac{(1-y_i)^2}{2L}\right)\mathbf{D}_y.
\end{align*}
$\mathbf{D}_x$ is the differentiation matrix with respect to $x$ and $\mathbf{D}_y$ is the differentiation matrix with respect to $y$. $\mathbf{Diag}\left(\dfrac{(1-y_i)^2}{2L}\right)$ represents a matrix with $\dfrac{(1-y_i)^2}{2L}$ on its diagonal and $0$ in all other entries. To clarify what I mean, I will show you what this matrix looks like:
\begin{align*}
\mathbf{Diag}\left(\dfrac{(1-y_i)^2}{2L}\right) &= \begin{bmatrix}
\dfrac{(1-y_0)^2}{2L} & 0 & 0 & 0 & ... & 0 \\
0 & \dfrac{(1-y_1)^2}{2L} & 0 & 0 & ... & 0 \\
0 & 0 & \dfrac{(1-y_2)^2}{2L} & 0 & ... & 0 \\
... & ... & ... & ... & ... & ... \\
0 & 0 & 0 & 0 & ... & \dfrac{(1-y_n)^2}{2L}
\end{bmatrix}.
\end{align*}
Why do we left-multiply $\mathbf{D}_y$ by this matrix? Well, if you manually calculate the matrix product between these matrices you will find that each row of $\mathbf{D}_y$ (which corresponds to a different value of $y_i$) is multiplied by the corresponding value of $\dfrac{(1-y_i)^2}{2L}$. This is the matrix analogue to our expression for $\dfrac{d}{dx}$ in terms of $y$.

Hence we will be solving the matrix equation:
\begin{align*}
\left(\mathbf{D}_x^2 - \mathbf{Diag}{(x_i)}\right)\mathbf{u} &= \mathbf{0}.
\end{align*}
Where $\mathbf{0}$ denotes a column vector of just zeros. Although, we will need to remove two rows from this matrix equation and replace them with rows to impose our boundary conditions. Say these rows are rows 1 and $n+1$, which correspond to $x_0$ and $x_n$. Say we call $\mathbf{A}=\mathbf{D}_x^2 - \mathbf{Diag}{(x_j)}$. Then the edits we make to impose our boundary conditions will be setting $A_{1,1} = 1$ and all remaining entries in the first row of $\mathbf{A}$ will be $0$. As for $A_{n+1,n+1}=1$ and all other entries in the $n+1$th row of $\mathbf{A}$ will be zero. 

Code Listing \coderef{chebyshev} was used to solve our Airy equation ODE.

\codelisting{chebyshev}{Solution to the Airy equation via Chebyshev spectral method.}{./code/chebyshev}{fidello/code/chebyshev.jl}
\output{./code/chebyshev}

An important issue I should probably address is how do we solve Sturm-Liouville and other differential equations with eigenvalues? Well, the answer is that we construct an eigenvector equation of the form:
\begin{align*}
\mathbf{A} \mathbf{u} = \lambda \mathbf{u}.
\end{align*}
Where $\mathbf{u}$ is our vector of values of our function at our grid values. Typically, we omit the endpoints from our matrices to impose our boundary conditions. Then we use established techniques for finding eigenvalues and eigenvectors. 