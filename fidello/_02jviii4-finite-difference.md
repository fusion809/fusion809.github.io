#### Finite difference method
So far I have covered how we can numerically integrate initial value problems. Solving boundary value problems is often more challenging, and we usually need to construct and solve matrix equations instead. For instance, if we had the equation $y''=-y$ and $y(0)=y(2\pi)=1$ (which implies $y(x)=\cos{x}$) we could use the finite difference method as it is called. In this, we discretize our grid, which means splitting our domain into steps of size $h$. We will define our grid points as:
\begin{align*}
x_i &= x_0 + hi.
\end{align*}
Where $x_0$ is the start of our domain, which is $0$. If we divide our domain into $n+1$ equally spaced points then $h=\dfrac{2\pi}{n}$. Hence:
\begin{align*}
x_i &= hi \\
&= \dfrac{2\pi i}{n}.
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
... \\
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

\codelisting{findiff}{Implementation of a second-order finite difference method for solving $y''=-y$ with $y(0)=y(2\pi)=1$.}{./code/findiff}{fidello/code/findiff.jl}

As we can see from the below output, our implementation of the finite difference method is rather inaccurate. Higher order approximations for the derivative can be achieved by using more grid point values at each step, but this does make the matrices more complex. 
\output{./code/findiff}

For instance, if we use the fourth-order approximation to $y''_i$:
\begin{align*}
y''_i &= \dfrac{-y_{i+2}+16y_{i+1}-30y_i + 16y_{i-1}-y_{i-2}}{12h^2}.
\end{align*}
For $i=1$ and $i=n-1$, where one of these grid points is not defined, we would instead use:
\begin{align*}
y''_1 &= \dfrac{11y_5 - 56y_4 + 114y_3 - 104y_2+35y_1}{12h^2} \\
y''_{n-1} &= \dfrac{35y_{n-1}-104y_{n-2}+114y_{n-3}-56y_{n-4}+11y_{n-5}}{12h^2}.
\end{align*}
This scheme was used in Code Listing \coderef{findiff4th}.

\codelisting{findiff4th}{Implementation of the fourth-order finite difference method for $y''=-y$, $y(0)=y(2\pi)=1$.}{./code/findiff4th}{fidello/code/findiff4th.jl}
Which produces the output:
\output{./code/findiff4th}
As we can see, it is more accurate, but definitely not perfect. 

\codelisting{findiff10th}{Implementation of the tenth-order finite difference method for $y''=-y$, $y(0)=y(2\pi)=1$.}{./code/findiff10th}{fidello/code/findiff10th.jl}
Which produces the output:
\output{./code/findiff10th}

As you can see from Code Listing \coderef{findiff10th}, I tried implementing higher-order finite difference estimates of the solution and an interesting trend emerged. My matrices became increasingly poorly conditioned. For the tenth-order scheme, for instance, we obtained a condition number of $6.11\times 10^{10}$. To put this into perspective, the fourth-order method's condition number was $1.51\times 10^6$. The second-order method's condition number was a mere 4,838.

It is important to note, however, that so far we have only applied the finite difference method to linear differential equations and this is for good reason. When we try to apply it to nonlinear differential equations, we get a nonlinear coupled system of algebraic equations which are much more challenging to solve. There are ways around this, for instance we can use Newton's method to solve the algebraic equation system that arises, and we can also try linearizing our differential equation using a technique known as the Newton-Kantorovich method \cite{boyd_chebyshev_2000}{526-535}. Each of these presents its own challenges though. 
