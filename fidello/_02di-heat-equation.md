### Heat equation
One of the simplest partial differential equations is the heat equation. In its simplest form it is:
\begin{align}
    \dfrac{\partial u}{\partial t} &= \alpha \dfrac{\partial^2 u}{\partial x^2} \label{heat}
\end{align}
Where $u(x,t)$ is a function describing the temperature of the one-dimensional object (only one dimensional as the above heat equation only has a single spatial dimension, if we want more spatial dimensions, we will have to add second partial derivatives with respect to those other dimensions, too) we are modelling and $\alpha$ is some constant. I know this equation may look confusing at first glance, but it is actually saying something rather intuitive to you. 

Do you remember from your maths class that the second derivative of a function is its concavity? So if the second derivative is negative, that means our function is concave down &mdash; it is leading to some sort of maximum. Whereas if it is positive, the function is concave up, which means it is going up from some minimum. We have the second derivative with respect to space on the right-hand side of Equation \eqref{heat}, so it essentially returns positive if we have a localized colder area and neighbouring hotter points. Whereas it returns negative if we have a really hot point surrounded by cooler points. The left-hand side of Equation \eqref{heat}, however, represents the rate of change of our function in time.

Equation \eqref{heat} therefore essentially says that if we are in a region where there is a really cold point surrounded by hotter points, the temperature will increase in the short-term for that colder point. Whereas the equation says that the temperature will decrease in a hotter region surrounded by colder points. If you have ever brought something hot out of the oven, or put something of room temperature into the fridge or freezer, this should make intuitive sense to you. 

Code Listing \coderef{heat} solves the heat equation in 3D for a $100^{\circ}C$ sphere in a $20^{\circ}C$ room.

\codelisting{heat}{Heat equation simulation code.}{./code/heat}{fidello/code/heat.jl}
\figure{heatan}{right}{800px}{/heat_equation_3d.gif}{Animation showing the heat distribution of an object and how it changes over time.}