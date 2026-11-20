## Derivation of Equations \eqref{ELE}
\figure{stationary-action}{right}{400px}{/fidello/Stationary%20action.svg}{Some of the paths a system could take.}

The derivation of Equations \eqref{ELE} requires knowledge of a field of mathematics known as the calculus of variations as well as a physics principle known as the stationary action principle. I will provide my own derivation here; although, I would suggest checking out [Veritasium's video on this topic](https://www.youtube.com/watch?v=Q10_srZ-pbs) as it provides a more accessible introduction to this topic.

There is a quantity known as the action &mdash; which is denoted as $S$ &mdash; that is defined as:
\begin{align}
    S &\stackrel{\text{def}}{=} \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) dt. \label{act}
\end{align}

One interesting property of the action is that nature consistently sets it to its stationary point, and this is regardless of context. Regardless of whether we are examining the motion of planets with Newton's second law, or studying the behaviour of black holes with general relativity, or checking out the behaviour of subatomic particles, there is always an action that is at its stationary point. If you would like an explanation for how this arises in classical mechanics, I suggest checking out [this excellent video](https://www.youtube.com/watch?v=QbnkIdw0HJQ) by Dr Jorge S Diaz. Traditionally, this stationary point is typically thought of as a minimum (hence why this is also known as the principle of **least** action), but strictly speaking it is not always a minimum. 

Of all the paths that an object could explore, the one that it actually explores is the one that sets the action to its stationary point. Figure \figref{stationary-action} provides some example paths that we examine when applying the principle of stationary action. In the terminology of the calculus of variations, we can express the fact the action is at its stationary point as:
\begin{align}
\delta S &= 0.\label{deltaS}
\end{align}
This is analogous to saying, in regular calculus, that our function's first derivative is zero. Which, as you may remember from your maths classes, is an equation we can solve to find the stationary point(s) of our function. The difference is that $S$ is not a function that takes some value of an independent variable (e.g. $x$) and returns a value, instead $S$ is a *functional*. Which means that it takes a function as its input and returns a value. In this case, the function it takes as its input are the $q_i(t)$ that specify our path. This is what distinguishes calculus of variations from regular calculus. In calculus of variations, we are working with functionals, whereas in regular calculus we are working with functions. 

Say we add functions $\delta q_i$ to $q_i$, thereby adjusting the path our object takes. Crucially, however, the end times and end points of our path will be unchanged. The new, adjusted action we will denote by $S'$ and it is:
\begin{align*}
S' &= \int_{t_0}^{t_f} L(q_i+\delta q_i, \dot{q}_i+\delta \dot{q}_i, t) dt
\end{align*}
We will expand $L$ using the first couple of terms of its Taylor series around $a=q_i$ and $a=\dot{q}_i$ (we have multiple $a$s as our function depends on multiple variables):
\begin{align*}
S' &\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \delta q_i + \dfrac{\partial L}{\partial \dot{q}_i}\delta \dot{q}_i dt.
\end{align*}
Next, we will apply integration by parts to the $\dfrac{\partial L}{\partial \dot{q}_i}\delta \dot{q}_i$ term:
\begin{align*}
S' &\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \delta q_i dt + \left[\sum_{i=1}^n\dfrac{\partial L}{\partial \dot{q}_i}\delta q_i\right]_{t_0}^{t_f} - \int_{t_0}^{t_f} \sum_{i=1}^n\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\delta q_i dt.
\end{align*}
At the end times of our path, $\delta q_i$ should be zero, as the endpoints are fixed. Hence the term $\left[\dfrac{\partial L}{\partial \dot{q}_i}\delta q_i\right]_{t_0}^{t_f}$ should evaluate to zero. Our expression then simplifies to:
\begin{align*}
S' &\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \delta q_i dt - \int_{t_0}^{t_f} \sum_{i=1}^n\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\delta q_i dt \\
&\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\delta q_i \left(\dfrac{\partial L}{\partial q_i} - \dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\right) dt.
\end{align*}
The first of these terms is present in $S$, too. Hence:
\begin{align*}
S' - S &\approx \int_{t_0}^{t_f}\sum_{i=1}^n\delta q_i \left(\dfrac{\partial L}{\partial q_i} - \dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\right) dt.
\end{align*}
We can denote this first-order approximation to how much $S$ changes as $\delta S$. It is a first-order approximation as this is the highest power of $\delta q_i$ and $\delta \dot{q}_i$ it depends on. First order also describes the highest order derivatives of $L$ with respect to $q_i$ and $\dot{q}_i$ that our approximation depends upon. Hence:
\begin{align}
\delta S &= \int_{t_0}^{t_f} \sum_{i=1}^n\delta q_i \left(\dfrac{\partial L}{\partial q_i} - \dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\right) dt. \label{deltaSELE}
\end{align}
Hence $\delta S$ will equal zero irrespective of $\delta q_i$ if and only if Equations \eqref{ELE} hold. Hence from setting $\delta S=0$ we get the Euler-Lagrange equations. 
