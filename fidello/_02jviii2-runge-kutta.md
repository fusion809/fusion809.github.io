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
