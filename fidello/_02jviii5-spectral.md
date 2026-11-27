#### Spectral methods
Sturm-Liouville problems also give rise to infinite sets of orthogonal functions that are useful when numerically approximating the solution to other differential equations \cite{boyd_chebyshev_2000}{66-67}. If $\{\phi_i\}$, $i\in \{0,1,2,...\}$ is our orthogonal function set, we assume our solution, $u$, to some ordinary differential equation can be written as:
\begin{align}
u(x) &= \sum_{j=0}^{\infty} a_j \phi_j(x). \label{spectral} 
\end{align}
Then we discretize our grid into points $x_i$ ($i \in \{0, 1, 2,...,n\}$) and truncate the above series expansion to $j\in \{0, 1, 2, ..., n\}$. What we do next depends on the precise method we are using. The simplest method we can use is to write a matrix equation for $\mathbf{a}$ (a vector containing $\{a_j\}$, the series coefficients). Then this equation would be:
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
In this case, we would be solving for $\mathbf{u}$ instead of $\mathbf{a}$. Methods that approximate the solution to a differential equation over the entire domain of the problem using a single sum of basis functions whose number of terms is equal to the number of grid points like this are called spectral methods. 

Why do we choose orthogonal function as our basis set for spectral methods? Well, orthogonality means our basis functions will be very different to one another and we want our basis functions to be as different from one another as possible to ensure that they can, as accurately as possible, represent all the variation in our solution. This difference helps spectral methods achieve exponential convergence instead of polynomial convergence. 

What does exponential and polynomial convergence mean? Remember how the truncation error for Euler's method, Runge-Kutta methods, Runge-Kutta-Fehlberg method, finite difference method, etc. were all of the form $\mathcal{O}(h^k)$ where $k$ is some fixed number that depends on the precise method used and $h$ is the step size? Well, this is what we call polynomial convergence. As $h \propto n^{-1}$ (where $n+1$ is our number of grid points), this can be rewritten as $\mathcal{O}(n^{-k})$. Spectral methods, however, have a truncation error of $O(e^{-qn^{r}})$ where $q>0$ and $r>0$. Why did we switch to using $n$ in our truncation error? As we will soon see, it is often convenient not to use a fixed step size $h$ for spectral methods and to instead use a set of $n+1$ grid points specific to the orthogonal function set we are using. 

As we can represent our spectral method solution to a differential equation as a function (the sum in Equation \eqref{spectral}) that approximates the true solution, spectral methods are sometimes called semi-analytical methods instead of true numerical methods which give us only approximated values of our solution.

If the domain is split up into smaller sections with a separate function of the form Equation \eqref{spectral} defined on each, this is called a spectral element method. Finite element methods are very similar as they do split the domain up and define a series-based solution on each part of the domain. The difference is that the series usually has fewer terms and achieves lower accuracy than their spectral element counterparts.

Like with the finite difference method, spectral methods, spectral element methods and finite element methods are also easiest applied to linear differential equations. This is because they result in a system of equations &mdash; if the differential equation is linear, this system is linear, but if not the system is nonlinear, which requires the application of [Newton's method](#newtons-method). Either that or we can linearize our differential equation using the [Newton-Kantorovich method](#newton-kantorovich_method).