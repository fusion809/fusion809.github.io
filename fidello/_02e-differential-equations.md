## Differential equations
As this article involves a lot of differential equations, I thought I should cover a little on how we categorize and solve differential equations. 

Differential equations are equations that involve the derivatives of a function and typically the aim in solving them is to determine the original function they pertain to. For instance, say you had the uncoupled, linear, homogeneous, first-order ordinary differential equation:
\begin{align*}
y' &= y.
\end{align*}
Our aim would be to determine what $y$ is. 

I also showed you this equation to clarify some of the terminology we use with regards to differential equations. So why is this equation *uncoupled*? Because it is a single equation and not a series of differential equations we need to solve simultaneously. Why is this equation *linear*? Well, it is because there is no term that is nonlinear in $y$ or its derivative(s). For instance, if the equation was $y'=y^2$, that would make it nonlinear because the right-hand side is nonlinear in $y$. Why is it an *ordinary* differential equation (ODE)? Because there are only ordinary derivatives, not partial derivatives. Why is it *homogeneous*? Because there are no terms that do not involve $y$ or its derivative(s). Why is it *first order*? Because the highest order derivative in the equation is the first derivative. 

In this example, there would be a few ways we could solve this differential equation. As luckily each of the adjectives we used to describe it implies it should be easy to solve. One option is we could divide both sides by $y$ and integrate on $x$:
\begin{align*}
\int \dfrac{\dfrac{dy}{dx}}{y} dx &= \int dx \\
\int \dfrac{dy}{y} &= \int dx \\
\ln{y} &= x + C
\end{align*}
Exponentiating both sides yields:
\begin{align*}
y &= e^{x+C}.
\end{align*}
Note that we have one unknown constant, $C$. This is actually a result of our equation being a first-order ordinary differential equation. $n$th-order ordinary differential equations give us $n$ unknown constants and we need to use any initial conditions or other conditions we are given to solve for them. For instance, if $y(0)=1$, then we can substitute this into our solution to find $C$.
\begin{align*}
1 &= e^{0+C} \\
e^C &= 1 \\
\implies C &= 0.
\end{align*}
Hence:
\begin{align*}
y &= e^x.
\end{align*}

### Numerical solution techniques
In this article, however, we will encounter differential equations that have no exact solutions and for which numerical solutions are necessary. Numerical solutions are actually quite interesting, because they typically involve applications of the Taylor series that we have covered. For instance, say we wanted to numerically integrate the differential equation $y'=y$ with initial condition $y(0)=1$. Then if we wanted to approximate $y(h)$ for some small $h$, we could use the Maclaurin series of $y$:
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
And continue on from there. This method, where we truncate the Taylor series at the first two terms and numerically approximate our solution by applying this approximation at each step, is called Euler's method. As we ignore so many terms in the Taylor series, this method is said to be rather inaccurate. Specifically, we call it a first-order method. 

Other schemes for numerical integration of ordinary differential equations also exist. Among the most popular are the Runge-Kutta family of methods, such as the Runge-Kutta fourth-order method. These methods do not utilize higher-order derivatives, but they use clever tricks for approximating the solution to higher-order accuracy using only first-order derivatives.

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

The methods for numerically integrating partial differential equations (PDEs) tend to be substantially more complex than those for solving ODEs due to the added complexity of needing to work with more than one independent variable. Although, sometimes these methods utilize ODE methods in some way. For instance, the code in Code Listing \coderef{heatan} solves the heat equation by approximating spatial derivatives using something called Fast Fourier Transforms (FFTs) and their inverses and then integrating on time using RKF45. This is one of many different techniques for numerically solving PDEs. Although, like ODEs, some PDEs also have analytical solutions and this includes the heat equation. 