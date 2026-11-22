### Numerical solution techniques
In this article, however, we will encounter differential equations that have no exact solutions and for which numerical solutions are necessary. Numerical solution methods approximate the solution to a differential equations, and as such they are not exact. Although, quite often we can make them accurate to ten or more significant figures. 

#### Euler's method
Numerical solutions are actually quite interesting, because they typically involve applications of the Taylor series that we have covered. For instance, say we wanted to numerically integrate the differential equation $y'=y$ with initial condition $y(0)=1$. Then if we wanted to approximate $y(h)$ for some small $h$, we could use the Maclaurin series of $y$:
\begin{align*}
y(h) &= y(0) + hy'(0) + \dfrac{h^2}{2}y''(0) + \dfrac{h^3}{6}y'''(0) + ... \\
\end{align*}
Now our differential equation makes it very easy to approximate any arbitrary-order derivative of $y$ at $x=0$, as we can just keep differentiating our differential equation and find that $\dfrac{d^n y}{dx^n} = y$ for any $n\in \mathbb{Z}^{+}$. But say we did not have that luxury. Then what we could do is just approximate $y(h)$ using the first two terms of the Maclaurin series:
\begin{align*}
y(h) &\approx y(0) + hy'(0) \\
&= 1 + h \cdot 1 \\
&= 1 + h.
\end{align*}
We can then approximate $y(2h)$ using the same method:
\begin{align*}
y(2h) &\approx y(h) + hy'(h) \\
&\approx (1+h) + h(1+h) \\
&= 1 + 2h + h^2.
\end{align*}
And continue on from there. This method, where we truncate the Taylor series at the first two terms and numerically approximate our solution by applying this approximation at each step, is called Euler's method. As we ignore so many terms in the Taylor series, this method is rather inaccurate. Specifically, we call it a first-order method. 

#### Runge-Kutta methods
Other schemes for numerical integration of ordinary differential equations also exist. Among the most popular are the Runge-Kutta family of methods, such as the Runge-Kutta fourth-order method. These methods do not utilize higher-order derivatives, but they use clever tricks for approximating the solution to higher-order accuracy using only first-order derivatives. I should say that Euler's method and the Runge-Kutta methods are all only applicable when we have initial conditions, not boundary or other conditions. There are other integration schemes for ODEs that work well with these other conditions, however.

To illustrate how we use these methods, let us try to solve the differential equation:

\begin{align*}
y' &= y, \hspace{0.2cm} y(0)=1.
\end{align*}

So our starting point is going to be $x=0$. We will call this point $x_0$; the value of $y$ at $x_0$ we will call $y_0$. The fourth-order Runge-Kutta method, says that we should approximate $y$ at $x_0+h$ for the differential equation:
\begin{align*}
y' &= f(x,y)
\end{align*}
using:
\begin{align*}
k_1 &= hf(x_0, y_0) \\
k_2 &= hf\left(x_0 + \dfrac{h}{2}, y_0 + \dfrac{k_1}{2}\right) \\
k_3 &= hf\left(x_0 + \dfrac{h}{2}, y_0 + \dfrac{k_2}{2}\right) \\
k_4 &= hf(x_0 + h, y_0 + k_3) \\
y(x_0+h) &= y_0 + \dfrac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4).
\end{align*}
And we should repeat this process for all subsequent points we wish to approximate $y(x)$ at. That is why, when we write numerical integration schemes like this we usually use $x_n$ to denote the $n$th $x$ value for which we are approximating $y$. 

In Code Listing \coderef{differential}, I have implemented both Euler's method and the Runge-Kutta fourth-order method to solve this differential equation. The results are displayed, along with the exact solution, in a table below the code. Here we used $h=0.2$ and continued the approximation for five steps, until $x=1$. As you can see, the Runge-Kutta fourth-order method is much more accurate than Euler's method.

\codelisting{differential}{Differential equation solution code.}{./code/differential}{fidello/code/differential.jl}

\output{./code/differential}

#### Runge-Kutta-Fehlberg 4/5th-order method
An arguably better method for numerically solving ordinary differential equations is the Runge-Kutta-Fehlberg fourth/fifth-order method (RKF45), as it adapts our step size $h$ to ensure that the error in our solution never exceeds a certain threshold. In this method, we use:

\begin{align*}
k_1 &= hf(x_n, y_n) \\
k_2 &= hf\left(x_n + \dfrac{h}{4}, y_n + \dfrac{k_1}{4}\right) \\
k_3 &= hf\left(x_n + \dfrac{3h}{8}, y_n + \dfrac{3k_1}{32} + \dfrac{9k_2}{32}\right) \\
k_4 &= hf\left(x_n + \dfrac{12h}{13}, y_n + \dfrac{1932k_1}{2197} - \dfrac{7200k_2}{2197} + \dfrac{7296k_3}{2197}\right) \\
k_5 &= hf\left(x_n + h, y_n + \dfrac{439k_1}{216} - 8k_2 + \dfrac{3680k_3}{513} - \dfrac{845k_4}{4104}\right) \\
k_6 &= hf\left(x_n + \dfrac{h}{2}, y_n - \dfrac{8k_1}{27} + 2k_2 - \dfrac{3544k_3}{2565} + \dfrac{1859k_4}{4104} - \dfrac{11k_5}{40}\right) \\
z_{n+1} &= y_n + \dfrac{16k_1}{135} + \dfrac{6656k_3}{12825} + \dfrac{28561k_4}{56430} - \dfrac{9k_5}{50} + \dfrac{2k_6}{55} \\
y_{n+1} &= y_n + \dfrac{25k_1}{216} + \dfrac{1408k_3}{2565} + \dfrac{2197k_4}{4104} - \dfrac{k_5}{5} \\
\mathrm{TE} &= |y_{n+1} - z_{n+1}| \\
s &= \left(\dfrac{\epsilon}{2\mathrm{TE}}\right)^{1/5}.
\end{align*}
Here $\epsilon$ is the maximum absolute error in our approximation we are willing to accept. $\mathrm{TE}$ is our approximation of the truncation error as our estimate for $y_{n+1}$ is fourth-order and $z_{n+1}$ is the fifth-order estimate. If $\mathrm{TE} \geq \epsilon$, then we rescale $h$ using $h=h\cdot s$ and repeat the approximation again. If not, then we accept our approximation to the next $y$ value and rescale $h$ using $h=h\cdot s$. This method is very efficient and is the numerical technique I use for integrating ordinary differential equations in all of my JavaScript-based ordinary differential equations solver webpages. 

\codelisting{RKF45}{Runge-Kutta-Fehlberg fourth/fifth-order method solution to $y'=y, y(0)=1$.}{./code/RKF45}{fidello/code/RKF45.jl}

\output{./code/RKF45}

As you can see RKF45 achieved more accuracy than both Euler's method and the Runge-Kutta fourth-order method, while using an, on average, larger step size. Also notice how we set the error tolerance of $10^{-5}$ and none of our errors exceeded this threshold. This is not always the case, however, as this error tolerance is purely local. This means that it is the maximum error we are willing to accept in each individual step of our solution. As these errors accumulate over multiple steps, it is possible for our solution to eventually deviate enough from the exact solution to exceed this threshold. 

Its main drawback is that because it is an adaptive method, the number of steps is unknown until you actually compute the solution and this can make implementation in code a little more challenging. As instead of dealing with vectors or arrays (containing your solution) of a set, predictable size, you will need to append new entries to the end of your solution vector/array at each step.

Further, the adaptive stepsize means that if you were hoping for an evenly spaced set of $x$ values, you may be disappointed by this method. That being said, it is straight forward to fix this issue as you can use interpolation to find the $y$ values on an evenly spaced $x$ grid. 

An interesting feature of higher-order ODEs is that we can rewrite them as a coupled system of first-order ODEs. Say, for instance, have the second-order ODE:

\begin{align*}
y'' = -y.
\end{align*}
Cosine and sine functions both satisfy this equation. The solution would be:
\begin{align*}
y &= C_1 \cos{x} + C_2 \sin{x}.
\end{align*}
So we would need two conditions to uniquely determine $y$. We could rewrite our second-order ODE as the system:
\begin{align*}
u &= y' \\
u' &= -y.
\end{align*}
How would we numerically integrate this system? We could apply RKF45 to determine $y_{n+1}$ (where $y_0$ and $y'_0$ will be our initial conditions) using this scheme:
\begin{align*}
k_1 &= hu_n\\
l_1 &= -hy_n \\
k_2 &= h\left(u_n+\dfrac{l_1}{4}\right) \\
l_2 &= -h\left(y_n+\dfrac{k_1}{4}\right) \\
k_3 &= h\left(u_n+\dfrac{3l_1}{32}+\dfrac{9l_2}{32}\right) \\
l_3 &= -h\left(y_n+\dfrac{3k_1}{32}+\dfrac{9k_2}{32}\right) \\
k_4 &= h\left(u_n+\dfrac{1932l_1}{2197}-\dfrac{7200l_2}{2197}+\dfrac{7296l_3}{2197}\right) \\
l_4 &= -h\left(y_n+\dfrac{1932k_1}{2197}-\dfrac{7200k_2}{2197}+\dfrac{7296k_3}{2197}\right) \\
k_5 &= h\left(u_n+\dfrac{439l_1}{216}-8l_2+\dfrac{3680l_3}{513}-\dfrac{845l_4}{4104}\right) \\
l_5 &= -h\left(y_n+\dfrac{439k_1}{216}-8k_2+\dfrac{3680k_3}{513}-\dfrac{845k_4}{4104}\right) \\
k_6 &= h\left(u_n-\dfrac{8l_1}{27}+2l_2-\dfrac{3544l_3}{2565}+\dfrac{1859l_4}{4104}-\dfrac{11l_5}{40}\right) \\
l_6 &= -h\left(y_n-\dfrac{8k_1}{27}+2k_2-\dfrac{3544k_3}{2565}+\dfrac{1859k_4}{4104}-\dfrac{11k_5}{40}\right) \\
u_{4,n+1} &= u_{n} + \dfrac{25l_1}{216} + \dfrac{1408l_3}{2565} + \dfrac{2197l_4}{4104} - \dfrac{l_5}{5}\\
u_{5,n+1} &= u_{n} + \dfrac{16l_1}{135} + \dfrac{6656l_3}{12825} + \dfrac{28561l_4}{56430} - \dfrac{9l_5}{50} + \dfrac{2l_6}{55}\\
y_{4,n+1} &= y_{n} + \dfrac{25k_1}{216} + \dfrac{1408k_3}{2565} + \dfrac{2197k_4}{4104} - \dfrac{k_5}{5}\\
y_{5,n+1} &= y_{n} + \dfrac{16k_1}{135} + \dfrac{6656k_3}{12825} + \dfrac{28561k_4}{56430} - \dfrac{9k_5}{50} + \dfrac{2k_6}{55}\\
\mathrm{TE} &= \max_{f \in \{u, y\}} \Big|f_{5,n+1}-f_{4,n+1}\Big|\\
s &= \left(\dfrac{\epsilon}{2\mathrm{TE}}\right)^{1/5}.
\end{align*}
I should clarify that we are maximizing over the error estimates in our two functions $u$ and $y$ to calculate $\mathrm{TE}$. So whichever error estimate is bigger out of those for $u$ and $y$, we choose for $\mathrm{TE}$. If $\mathrm{TE} < \epsilon$, we continue onto the next step with our new $y_n$ being the old $y_{5,n+1}$ and the new $u_n$ being the old $u_{5,n+1}$ and obtain our new $h$ by multiplying it by $s$. If $TE \geq \epsilon$, we multiply $h$ by $s$ and repeat the above calculations again.

So we will use the initial conditions $y(0)=0$ and $y'(0)=1$. This gives the solution $y(x)=\sin{x}$ and integrate our ordinary differential equation on the interval $x\in [0, 2\pi]$.

\codelisting{RKF45coup}{Runge-Kutta-Fehlberg 4/5th-order method solution to $y''=-y$ with $y(0)=0$ and $y'(0)=1$}{./code/RKF45coup}{fidello/code/RKF45coup.jl}
\output{./code/RKF45coup}

In Code Listing \coderef{RKF45coup}, we used $\epsilon=10^{-8}$ as our error tolerance. You may notice that our maximum errors in the output table do indeed exceed this. This is possible as our truncation error estimate is solely for that step, it is not the cumulative error for all steps up until the point we are in our calculation. 

#### Finite difference method
So far I have covered how we can numerically integrate initial value problems. Solving boundary value problems is often more challenging, and we usually need to construct and solve matrix equations instead. For instance, if we had the equation $y''=-y$ and $y(0)=y(2\pi)=1$ (which implies $y(x)=\cos{x}$) we could use the finite difference method as it is called. In this, we discretize our grid, which means splitting our domain into steps of size $h$. We will define our grid points as:
\begin{align*}
x_i &= x_0 + hi.
\end{align*}
Where $x_0$ is the start of our domain, which is $0$. If we divide our domain into $N+1$ equally spaced points then $h=\dfrac{2\pi}{N}$. Hence:
\begin{align*}
x_i &= hi \\
&= \dfrac{2\pi i}{N}.
\end{align*}
Let $y_i$ denote $y(x_i)$. Then using the Taylor series expansion at $a=x_i$:
\begin{align*}
y_{i+1} &= y_i + y'_i h + \dfrac{y''_i}{2}h^2 + \dfrac{y'''_i}{6} h^3 + \mathcal{O}(h^4)\\
y_{i-1} &= y_i - y'_i h + \dfrac{y''_i}{2} h^2 - \dfrac{y'''_i}{6} h^3 + \mathcal{O}(h^4)\\
y_{i+1}-y_{i-1} &= y_i + y'_i h + \dfrac{y''_i}{2}h^2  + \dfrac{y'''_i}{6} h^3- \left(y_i - y'_i h + \dfrac{y''_i}{2} h^2 - \dfrac{y'''_i}{6} h^3\right) + \mathcal{O}(h^4) \\
&= 2y'_i h + \dfrac{y'''_i}{3} h^3 + \mathcal{O}(h^4) \\
&= 2y'_i h + \mathcal{O}(h^3).
\end{align*}
Now we need to explain what $\mathcal{O}(h^4)$ is. It essentially refers to terms proportional to $h^4$ and higher powers of $h$. The reason they did not cancel out when we calculated $y_{i+1}-y_{i-1}$ is because they would differ between $y_{i+1}$ and $y_{i-1}$. Hence we can approximate $y'_i$ as:
\begin{align*}
y'_i &= \dfrac{y_{i+1}-y_{i-1}}{2h} + \dfrac{\mathcal{O}(h^3)}{2h} \\
&=\dfrac{y_{i+1}-y_{i-1}}{2h} + \mathcal{O}(h^2).
\end{align*}
If we add $y_{i+1}$ and $y_{i-1}$, we should cancel out the $y'_i h$ and $\dfrac{y'''_i}{6} h^3$ terms. To get rid of $y_i$, we will subtract $2y_i$:
\begin{align*}
y_{i+1}-2y_i+y_{i-1} &= y''_i h^2 + \mathcal{O}(h^4) \\
y''_i &= \dfrac{y_{i+1}-2y_i+y_{i-1}}{h^2} + \dfrac{\mathcal{O}(h^4)}{h^2} \\
&= \dfrac{y_{i+1}-2y_i+y_{i-1}}{h^2} + \mathcal{O}(h^2).
\end{align*}
Hence we have estimates for both $y''_i$ and $y'_i$ that have local truncation errors of $\mathcal{O}(h^2)$. Consequently, they are second-order estimates for $y''_i$ and $y'_i$, respectively. Substituting into our differential equation, ignoring the $\mathcal{O}(h^2)$ yields and adding our boundary conditions:
\begin{align*}
y''_i &= y_i \\
\dfrac{y_{i+1}-2y_i+y_{i-1}}{h^2} &= y_i \\
y_{i+1}-2y_i+y_{i-1} &= h^2 y_i \\
y_{i+1}-(2+h^2)y_i + y_{i-1} &= 0 \\
y_0 &= 1 \\
y_N &= 1.
\end{align*}
We can now write this in matrix form as:
\begin{align*}
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & ... \\
1 & -2+h^2 & 1 & 0 & 0 & ... \\
0 & 1 & -2+h^2 & 1 & 0 & ... \\
... & ... & ... & ... & ... & ... \\
0 & 0 & 0 & 0 & ... & 1
\end{bmatrix} \begin{bmatrix}
y_0 \\
y_1 \\
y_2 \\
...
y_n
\end{bmatrix} &= \begin{bmatrix}
1 \\
0 \\
0 \\
0 \\
... \\
1
\end{bmatrix}.
\end{align*}
Code Listing \coderef{findiff} contains code for implementing the finite difference method for our ODE. 

\codelisting{findiff}{Implementation of the finite difference method for solving $y''=-y$ with $y(0)=y(2\pi)=1$.}{./code/findiff}{fidello/code/findiff.jl}

As we can see from the below output, our implementation of the finite difference method is rather inaccurate. Higher order approximations for the derivative can be achieved by using more grid point values at each step, but this does not make the matrices more complex. 
\output{./code/findiff}

For instance, if we use the fourth-order approximation to $y''_i$:
\begin{align*}
y''_i &= \dfrac{-y_{i+2}+16y_{i+1}-30y_i + 16y_{i-1}-y_{i-2}}{12h^2}.
\end{align*}
For $i=1$ and $i=N-1$, where one of these grid points is not defined, we would instead use:
\begin{align*}
y''_1 &= \dfrac{11y_5 - 56y_4 + 114y_3 - 104y_2+35y_1}{12h^2} \\
y''_{N-1} &= \dfrac{35y_{N-1}-104y_{N-2}+114y_{N-3}-56y_{N-4}+11y_{N-5}}{12h^2}.
\end{align*}
These scheme was used in Code Listing \coderef{findiff4th}.

\codelisting{findiff4th}{Implementation of the fourth-order finite difference method for $y''=-y$, $y(0)=y(2\pi)=1$.}{./code/findiff4th}{fidello/code/findiff4th.jl}
Which produces the output:
\output{./code/findiff4th}
As we can see, it is more accurate, but definitely not perfect. 

It is important to note, however, that so far we have only applied the finite difference method to linear differential equations and this is for good reason. When we try to apply it to nonlinear differential equations, we get a nonlinear coupled system of algebraic equations which are much more challenging to solve. There are ways around this, for instance we can use a method called Newton's method to solve the algebraic equation system that arises, and we can also try linearizing our differential equation using a technique known as the Newton-Kantorovich method. Each of these presents its own challenges though. 


#### Spectral methods
Sturm-Liouville problems also give rise to infinite sets of orthogonal functions that are useful when numerically approximating the solution to other differential equations. If $\{\phi_i\}$, $i\in {0,1,2,...}$ is our orthogonal function set, we assume our solution, $u$, to some ordinary differential equation can be written as:
\begin{align*}
u(x) &= \sum_{j=0}^{\infty} a_j \phi_j(x). 
\end{align*}
Then we discretize our grid into points $x_i$, truncate the above series expansion, and can write a matrix equation for $\mathbf{a}$ (a vector containing $\{a_j\}$, the series coefficients). Then this equation would be:
\begin{align}
\mathbf{A}\mathbf{a} &= \mathbf{b}. \label{speca}
\end{align}
Let $\boldsymbol{\phi}$ be a matrix containing $\phi_j(x_i)$ and $\mathbf{u}$ be a vector of $u(x_i)$ values. Then we can write:
\begin{align*}
\mathbf{a} &=\boldsymbol{\phi}^{-1} \mathbf{u}.\label{afromu}
\end{align*}
We can substitute Equation \eqref{afromu} into Equation \eqref{speca}, yielding:
\begin{align*}
\mathbf{A}\boldsymbol{\phi}^{-1}\mathbf{u} &= \mathbf{b}.
\end{align*}
In this case, we would be solving for $\mathbf{u}$ instead of $\mathbf{a}$. Methods that approximate the solution to a differential equation using a sum of basis functions like this are called spectral methods. 

Many different orthogonal functions can be used; among the most popular are Chebyshev polynomials of the first kind and Fourier basis functions. Chebyshev polynomials of the first kind, $T_i(x)$, are defined on $x\in [-1,1]$ as:
\begin{align*}
T_j(x) = \cos(j\cos^{-1}{x}).
\end{align*}
The grid points normally used as:
\begin{align*}
x_i &= \cos{\dfrac{\pi i}{n}}.
\end{align*}
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
\mathbf{D}_x &= \mathbf{Diag}\left(\dfrac{(1-y_j)^2}{2L}\right)\mathbf{D}_y.
\end{align*}
Where $\mathbf{Diag}\left(\dfrac{(1-y_j)^2}{2L}\right)$ represents a matrix with $\dfrac{(1-y_j)^2}{2L}$ on its diagonal and $0$ in all other entries. $\mathbf{D}_x$ is the differentiation matrix with respect to $x$ and $\mathbf{D}_y$ is the differentiation matrix with respect to $y$. Hence we will be solving the matrix equation:
\begin{align*}
\left(\mathbf{D}_x^2 - \mathbf{Diag}{(x_j)}\right)u &= 0.
\end{align*}
Although, we will need to remove two rows from this matrix equation and replace them with rows to impose our boundary conditions. Say these rows are rows 1 and $n+1$, which correspond to $x_0$ and $x_n$. Say we call $\mathbf{A}=\mathbf{D}_x^2 - \mathbf{Diag}{(x_j)}$. Then the edits we make to impose our boundary conditions will be setting $A_{1,1} = 1$ and all remaining entries in the first row of $\mathbf{A}$ will be $0$. As for $A_{n+1,n+1}=1$ and all other entries in the $n+1$th row of $\mathbf{A}$ will be zero. 

\codelisting{chebyshev}{Solution to the Airy equation via Chebyshev spectral method.}{./code/chebyshev}{fidello/code/chebyshev.jl}
\output{./code/chebyshev}

#### Partial differential equations
The methods for numerically integrating partial differential equations (PDEs) tend to be substantially more complex than those for solving ODEs due to the added complexity of needing to work with more than one independent variable. Although, the finite difference method can be applied directly to boundary value linear PDEs, too. Although the matrices become more complex when we add additional independent variables. Further, some of the other methods used to numerically integrate PDEs utilize ODE methods in some way. 

For instance, the code in Code Listing \coderef{heatan} solves the heat equation by approximating spatial derivatives using something called the fast Fourier transform (FFT) and its inverse and then integrates on time using RKF45. This is one of many different techniques for numerically solving PDEs.
