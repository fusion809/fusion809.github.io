## Partial derivatives
Partial derivatives are just a straightforward extension of derivatives as you understand them already. The key difference is that they are derivatives of multivariable functions instead of single-variable functions. Further, when you differentiate with respect to one variable, all other variables are kept constant. If you would like to build a graphical image of what the partial derivative is, well recall that the derivative of a single-variable function is just the gradient of that function. Consequently, it is also the gradient of the tangent line to that function. The partial derivative is similar to that. As when we have a function with more than a single independent variable, we have the freedom to choose which variable or variables we wish to determine our function's gradient with respect to. 

If you want to visualize this, imagine a surface defined by $z=f(x,y)$ for example. Because a tangent line is two dimensional and we're in three dimensional space, there can be more than one tangent line to a given point. As when constructing the tangent line, we must answer a question: in which direction do we want our line to be tangential to our surface? It could be tangential in the $x$ direction or the $y$ direction, or even in any arbitrary direction in the $xy$ plane. If we want the tangent line in the $x$ direction, we will need the partial derivative with respect to $x$ to get its gradient, whereas if we want the tangent line in the $y$ direction, we will need the partial derivative with respect to $y$ to get its gradient. More generally, for an arbitrary direction in the $xy$ plane we would use the directional derivative as it is called. 

Say $f$ is a function of two variables &mdash; $x$ and $y$. Then its partial derivative with respect to $x$ is formally defined as:

\begin{align*}
    \dfrac{\partial f}{\partial x} &\stackrel{\text{def}}{=} \lim_{h\rightarrow 0} \dfrac{f(x+h, y)-f(x, y)}{h}.
\end{align*}

This is analogous to the first principles definition of the derivative you would have learned for standard derivatives. Namely, if we are differentiating $f(x)$:
\begin{align*}
    \dfrac{df}{dx} &= \lim_{h\rightarrow 0}\dfrac{f(x+h) - f(x)}{h}.
\end{align*}
The only difference is that now we also have a variable $y$ that our function $f(x,y)$ depends on, but is kept constant when differentiating with respect to $x$.
Say $f(x,y)=2x^2y$. Then $\dfrac{\partial f}{\partial x}$ would be:
\begin{align*}
    \dfrac{\partial f}{\partial x} &= \lim_{h\rightarrow 0} \dfrac{2(x+h)^2y - 2x^2y}{h}.
\end{align*}
Did you notice that $2y$ is a common factor of both terms in our numerator? Let us take it out of our limit altogether, as it is effectively a constant (it does not depend on $h$, after all). 

\begin{align*}
\dfrac{\partial f}{\partial x} &= 2y\lim_{h\rightarrow 0} \dfrac{(x+h)^2-x^2}{h}
\end{align*}
Using the binomial theorem, we know that
\begin{align*}
(x+h)^2 &= x^2 + 2xh + h^2.
\end{align*}

Hence:
\begin{align*}
\dfrac{\partial f}{\partial x} &= 2y\lim_{h\rightarrow 0} \dfrac{x^2+2xh + h^2-x^2}{h} \\
&= 2y\lim_{h\rightarrow 0} \dfrac{2xh + h^2}{h} \\
&= 2y\lim_{h\rightarrow 0} 2x + h \\
&= 2y \cdot 2x \\
&= 4xy.
\end{align*}
Okay, so that is the long way of finding out partial derivatives. Fortunately, there is also a faster way. Which is that we apply normal differentiation rules, but any variable we are not differentiating with respect to is treated like it is a constant. So for $f(x,y)=2x^2y$ we would treat the $2y$ as just a constant number and only differentiate the $x^2$. Applying the power rule, we see that $\dfrac{d}{dx}x^2 = 2x$ and hence:
\begin{align*}
    \dfrac{\partial f}{\partial x} &= 2y \cdot 2x \\
    &= 4xy.
\end{align*}

### Heat equation
One of the simplest partial differential equations is the heat equation. In its simplest form it is:
\begin{align}
    \dfrac{\partial u}{\partial t} &= \alpha \dfrac{\partial^2 u}{\partial x^2} \label{heat}
\end{align}
Where $u(x,t)$ is a function describing the temperature of the one-dimensional object (only one dimensional as the above heat equation only has a single spatial dimension, if we want more spatial dimensions we'll have to add second partial derivatives with respect to those other dimensions, too) we're modelling and $\alpha$ is some constant. I know this equation may look confusing at first glance, but it is actually saying something rather intuitive to you. 

Do you remember from your maths class that the second derivative of a function is its concavity? So if the second derivative is negative, that means our function is concave down &mdash; it is leading to some sort of maximum. Whereas if it is positive, the function is concave up, which means it is going up from some minimum. We have the second derivative with respect to space on the right-hand side of Equation \eqref{heat}, so it essentially returns positive if we have a localized colder area and neighbouring hotter points. Whereas it returns negative if we have a really hot point surrounded by cooler points. The left-hand side of Equation \eqref{heat}, however, represents the rate of change of our function in time.

Equation \eqref{heat} therefore essentially says that if we are in a region where there is a really cold point surrounded by hotter points, the temperature will increase in the short-term for that colder point. Whereas the equation says that the temperature will decrease in a hotter region surrounded by colder points. If you've ever brought something hot out of the oven, or put something of room temperature into the fridge or freezer, this should make intuitive sense to you. 
