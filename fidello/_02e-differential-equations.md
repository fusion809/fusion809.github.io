## Differential equations
As this article involves a lot of differential equations, I thought I should cover a little on how we categorize and solve differential equations. 

Differential equations are equations that involve the derivatives of a function and typically the aim in solving them is to determine the original function they pertain to. For instance, say you had the linear, homogeneous, first-order ordinary differential equation:
\begin{align*}
y' &= y
\end{align*}.
Our aim would be to determine what $y$ is. 

I also showed you this equation to clarify some of the terminology we use with regards to differential equations. So why is this equation *linear*? Well, it is because there is no term that is nonlinear in $y$ or its derivative(s). For instance, if the equation was $y'=y^2$, that would make it nonlinear because the right-hand side is nonlinear in $y$. Why is it an *ordinary* differential equation? Because there are only ordinary derivatives, not partial derivatives. Why is it *homogeneous*? Because there are no terms that do not involve $y$ or its derivative(s). Why is it first-order? Because the highest order derivative in the equation is the first derivative. 

In this example, there would be a few ways we could solve this differential equation. As luckily each of the adjectives we used to describe it, implies it should be easy to solve. One option is we could divide both sides by $y$ and integrate on $x$:
\begin{align*}
\int \dfrac{\dfrac{dy}{dx}}{y} dx &= \int dx \\
\int \dfrac{dy}{y} &= \int dx \\
\ln{y} &= x + C \\
\therefore y &= e^{x+C}.
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

In this article, however, we will encounter that have no exact solutions and for which numerical solutions are necessary. Numerical solutions are actually quite interesting, because they typically involve applications of the Taylor series that we just covered. For instance, say we wanted to numerically integrate the differential equation $y'=y$ with initial condition $y(0)=1$. Then if we wanted to approximate $y(h)$ for some small $h$, we could use the Maclaurin series of $y$:
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
And continue on from there. This method, where we truncate the Taylor series at the first two terms and numerically approximate our solution by applying this approximation at each step is called Euler's method. As we ignore so many terms in the Taylor series, this method is said to be rather inaccurate. Specifically, we call it a first-order method. 

Other schemes for numerical integration of ordinary differential equations also exist. Among the most popular are the Runge-Kutta family of methods, such as the Runge-Kutta fourth-order method. These methods also derive, ultimately, from the Taylor series expansion of our function but they are derived from it in more convoluted ways.