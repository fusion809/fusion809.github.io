## Newton's method
Given we do later cover the Newton-Kantorovich method, which is closely related to Newton's method, I thought I should probably cover Newton's method here. Newton's method is a method for solving nonlinear algebraic equations, including systems of nonlinear algebraic equations. For instance, say we wanted to solve:
\begin{align*}
x = \cos^2{x}.
\end{align*}
Then Newton's method would provide us a way to do so. Newton's method essentially finds the root of a function by looking at where its tangent line intersects the $x$ axis. Another way to think about it is that it uses the first two terms (the constant term and the linear term) of the function's Taylor series to determine its $x$ intercept. 

For our example, we would rearrange it to make it a function whose root we are finding, namely:
\begin{align*}
f(x) &= x -\cos^2{x} \\
&= 0.
\end{align*}
Then Newton's method says we can apply this scheme to find its zero:
\begin{align*}
x_{n+1}=x_n - \dfrac{f(x_n)}{f'(x_n)}.
\end{align*}
Where $x_n$ is our $n$th approximation to $x$. $x_0$ will be our first guess as to $f(x)$'s root. For instance, I would guess $x\in[0,1]$ would be where we would find the root for $f(x)$. $\cos^2{x}$ is always between $0$ and $1$ for all real $x$, hence $x\in[0,1]$ sounds like a viable solution. Let us try $x=1$ as our first guess. Differentiating $f(x)$ yields:
\begin{align*}
f'(x) &= 1+2\sin{x}\cos{x} \\
&= 1+\sin{2x}.
\end{align*}
Hence Newton's method would come down to:
\begin{align*}
x_{n+1} &= x_n - \dfrac{x_n-\cos^2{x_n}}{1+\sin{2x_n}}.
\end{align*}
For $n=0$ we get:
\begin{align*}
x_1 &= 1 - \dfrac{1-\cos^2{1}}{1+\sin{2}} \\
&\approx 0.6291445.
\end{align*}
For $n=1$ and $n=2$, we get:
\begin{align*}
x_2 &= x_1 - \dfrac{x_1-\cos^2{x_1}}{1+\sin{2x_1}} \\
&= 0.6417386 \\
x_3 &= x_2 - \dfrac{x_2-\cos^2{x_2}}{1+\sin{2x_2}} \\
&= 0.6417144.
\end{align*}
If we continue further iterations of Newton's method, we will find that $x_n$ stays $0.6417144$ when rounding to seven decimal places. 

How do we extend this method to solving a system of algebraic equations? We use \cite{boyd_chebyshev_2000}{527}:
\begin{align*}
\mathbf{x}_{n+1} &= \mathbf{x}_n - \mathbf{J}_n^{-1}\mathbf{f}_n.
\end{align*}
$\mathbf{x}_{n}$ is a column vector of our $n$th approximation to the independent variables that our system of equations depends on. $\mathbf{f}_n$ is our system of equations written as a vector of functions we are simultaneously setting to zero evaluated at $\mathbf{x}=\mathbf{x}_n$. $\mathbf{J}_n$ is the Jacobian matrix evaluated at $\mathbf{x}=\mathbf{x}_n$. Let the elements of $\mathbf{x}_n$ and $\mathbf{f}_n$ be numbered $1$ to $m$ (where $m$ is the number of equations we are solving). What is the Jacobian matrix? Well, it is a matrix with elements \cite{boyd_chebyshev_2000}{526}:
\begin{align*}
\mathbf{J} &= \begin{bmatrix}
\dfrac{\partial f_1}{\partial x_1} & \dfrac{\partial f_1}{\partial x_2} & \dfrac{\partial f_1}{\partial x_3} & ... & \dfrac{\partial f_1}{\partial x_m} \\
\dfrac{\partial f_2}{\partial x_1} & \dfrac{\partial f_2}{\partial x_2} & \dfrac{\partial f_2}{\partial x_3} & ... & \dfrac{\partial f_2}{\partial x_m} \\
\dfrac{\partial f_3}{\partial x_1} & \dfrac{\partial f_3}{\partial x_2} & \dfrac{\partial f_3}{\partial x_3} & ... & \dfrac{\partial f_3}{\partial x_m} \\
... & ... & ... & ... & ... \\
\dfrac{\partial f_m}{\partial x_1} & \dfrac{\partial f_m}{\partial x_2} & \dfrac{\partial f_m}{\partial x_3} & ... & \dfrac{\partial f_m}{\partial x_m} \\
\end{bmatrix}.
\end{align*}
Which can be written more succinctly as:
\begin{align*}
J_{ij} = \left(\dfrac{\partial f_i}{\partial x_j}\right).
\end{align*}

