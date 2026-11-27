##### Example: Airy equation
Say we wish to solve:
\begin{align*}
u''-xu &= 0
\end{align*}
on $x\in[0,\infty]$. The exact solution to this differential equation is:
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
\dfrac{d}{dx} &= \dfrac{dy}{dx} \dfrac{d}{dy} \\
&= \dfrac{1}{\dfrac{dx}{dy}} \dfrac{d}{dy} \\
&= \dfrac{1}{x'(y)} \dfrac{d}{dy} \\
&= \dfrac{1}{\dfrac{L[(1-y)+(1+y)]}{(1-y)^2}} \dfrac{d}{dy} \\
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