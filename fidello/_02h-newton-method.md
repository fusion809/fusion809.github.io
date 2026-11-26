## Newton's method
Given we do later cover the [Newton-Kantorovich method](#newton-kantorovich_method), which is closely related to Newton's method, I thought I should probably cover Newton's method here. Newton's method, which is also known as the Newton-Raphson method, is a method for solving nonlinear algebraic equations, including systems of nonlinear algebraic equations. For instance, say we wanted to solve:
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
Where $x_n$ is our $n$th approximation to $x$. $x_0$ will be our first guess as to $f(x)$'s root. Newton's method will not always converge towards a solution if our initial guess is bad enough. Getting back to our example root-finding problem, I would guess $x\in[0,1]$ would be where we would find the root for $f(x)$. $\cos^2{x}$ is always between $0$ and $1$ for all real $x$, hence $x\in[0,1]$ sounds like a viable solution. Let us try $x=1$ as our first guess. Differentiating $f(x)$ yields:
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

For example, say we have the equations:
\begin{align*}
xy &= \log{(x+y)} \\
x^2 + y^2 = 9.
\end{align*}
First we would rewrite these equations to set one side of them to zero:
\begin{align*}
f_1(x,y) &= xy - \log{(x+y)} \\
&= 0 \\
f_2(x,y) &= x^2 + y^2 - 9 \\
&= 0.
\end{align*}
Hence $\mathbf{f}$ is:
\begin{align*}
\mathbf{f} &= \begin{bmatrix}
xy - \log{(x+y)} \\
x^2 + y^2 - 9
\end{bmatrix}.
\end{align*}
And our Jacobian is:
\begin{align*}
\mathbf{J} &= \begin{bmatrix}
\dfrac{\partial f_1}{\partial x} & \dfrac{\partial f_1}{\partial y} \\
\dfrac{\partial f_2}{\partial x} & \dfrac{\partial f_2}{\partial y}
\end{bmatrix} \\
&= \begin{bmatrix}
y - \dfrac{1}{x_n+y_n} & x- \dfrac{1}{x_n+y_n} \\
2x & 2y
\end{bmatrix}.
\end{align*}
Hence Newton's method boils down to applying:
\begin{align*}
\mathbf{x}_{n+1} &= \mathbf{x}_n - \mathbf{J}_n^{-1}\mathbf{f} \\
\begin{bmatrix}
x_{n+1} \\
y_{n+1}
\end{bmatrix} &= \begin{bmatrix}
x_n \\
y_n
\end{bmatrix} - \begin{bmatrix}
y_n - \dfrac{1}{x_n+y_n} & x_n- \dfrac{1}{x_n+y_n} \\
2x_n & 2y_n
\end{bmatrix}^{-1}\begin{bmatrix}
x_n y_n - \log{(x_n+y_n)} \\
x_n^2 + y_n^2 - 9
\end{bmatrix}.
\end{align*}
We can rewrite $y$ in terms of $x$ in the second equation, yielding:
\begin{align*}
y &= \pm \sqrt{9-x^2}.
\end{align*}
Hence $x\in [-3,3]$ and $y\in [-3,3]$. $x+y>0$ thanks the logarithm in the first equation, hence $y>-x$. Let us choose $x=2$ and $y=+\sqrt{9-2^2}=\sqrt{5}$ as a first guess. 

\codelisting{newtonsystem}{Newton's method solution to $\mathbf{f}=\mathbf{0}$ with initial guess $x=2$ and $y=\sqrt{5}$.}{./code/newtonsystem}{fidello/code/newton-system.jl}
\output{./code/newtonsystem}

As we can see above, $x=2$ and $y=\sqrt{5}$ are actually pretty poor first guesses, but our solution still luckily converged. 