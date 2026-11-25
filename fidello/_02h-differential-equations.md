## Differential equations
As this article involves a lot of differential equations, I thought I should cover a little on how we categorize and solve differential equations. This is definitely not intended to be a comprehensive covering of this topic. As this is a topic so broad that entire textbooks are written about it. 

Differential equations are equations that involve the derivatives of a function and typically the aim in solving them is to determine the original function they pertain to. For instance, say you had the uncoupled, linear, homogeneous, first-order ordinary differential equation:
\begin{align}
y' &= y. \label{firstDE}
\end{align}
Our aim would be to determine what $y$ is. 

I also showed you this equation to clarify some of the terminology we use with regards to differential equations. So why is this equation *uncoupled*? Because it is a single equation and not a series of differential equations we need to solve simultaneously. Why is this equation *linear*? Well, it is because there is no term that is nonlinear in $y$ or its derivative(s). For instance, if the equation was $y'=y^2$, that would make it nonlinear because the right-hand side is nonlinear in $y$. Why is it an *ordinary* differential equation (ODE)? Because there are only ordinary derivatives, not partial derivatives. Why is it *homogeneous*? Because there are no terms that do not involve $y$ or its derivative(s). Why is it *first order*? Because the highest order derivative in the equation is the first derivative. 

I should also tell you that because differentiating deletes constants from our function, in order to obtain a unique solution to a given differential equation of order $n$, we need $n$ additional conditions beyond the differential equation itself. Often these conditions are initial conditions, that is they specify the initial state of our system. In this event, our differential equation is called an initial value problem. For instance, if we wanted to study $\ddot{x}=-g$ over $t \in [0, 10]$ then the appropropriate initial conditions could be $x(0)=50$ and $\dot{x}(0)=0$, as these specify the initial state of our system. 

Other conditions are possible, too. For instance, boundary conditions could also be used for our ODE $\ddot{x}=-g$. For this problem, the appropriate boundary conditions would be the values of $x(0)$ and $x(10)$. In this case, we would call the problem a boundary value problem. We could also have the conditions $x'(2)=50$ and $x(5)=10$, which together do not really fit neatly into the definitions of either initial or boundary conditions. 