## Derivation of Equations \eqref{ELE}
~~~
<figure id="stationary-action" style="float: right;">
<img src="/fidello/Stationary%20action.svg" width="400px"/>
<figcaption>Some of the paths a system could take.</figcaption>
</figure>
~~~
The derivation of Equations \eqref{ELE} requires knowledge of a field of mathematics known as the calculus of variations as well as a physics principle known as the stationary action principle. If you are interested in learning more, I would suggest you check out Veritasium's video on this topic ([https://www.youtube.com/watch?v=Q10_srZ-pbs](https://www.youtube.com/watch?v=Q10_srZ-pbs)). 

There is a quantity known as the action, which is denoted as $S$, that is defined as:
\begin{align}
    S &= \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) dt. \label{act}
\end{align}

One interesting property of the action is that nature consistently sets it to its stationary point, and this is regardless of context. Regardless of whether we are examining the motion of planets with Newton's second law, or studying the behaviour of black holes with general relativity, or checking out the behaviour of subatomic particles, there is always an action that is at its stationary point. Traditionally, this stationary point is typically thought of as a minimum (hence why this is also known as the principle of **least** action), but strictly speaking it is not always a minimum. When you use the calculus of variations to find the stationary point of Equation \eqref{act}, Equations \eqref{ELE} pop out. 

I will also provide you my own derivation of Equations \eqref{ELE} using the principle of stationary action. Of all the paths that an object could explore, the one that it actually explores is the one that sets the action to its stationary point. Figure \figref{stationary-action} provides some example paths that we examine when applying the principle of stationary action. In the terminology of the calculus of variations, we can express the fact the action is at its stationary point as:
\begin{align}
\delta S &= 0.\label{deltaS}
\end{align}
This is analogous to saying, in regular calculus, that our function's first derivative is zero. 

Say we add functions $\delta q_i$ to $q_i$, thereby adjusting the path our object takes. Crucially, however, the end times and end points of our path will be unchanged. The new, adjusted action we will denote by $S'$ and it is:
\begin{align*}
S' &= \int_{t_0}^{t_f} L(q_i+\delta q_i, \dot{q}_i+\delta \dot{q}_i, t) dt \\
&\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \delta q_i + \dfrac{\partial L}{\partial \dot{q}_i}\delta \dot{q}_i dt.
\end{align*}

Next, we will apply integration by parts to the $\dfrac{\partial L}{\partial \dot{q}_i}\delta \dot{q}_i$ term:
\begin{align*}
S' &\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \delta q_i dt + \left[\sum_{i=1}^n\dfrac{\partial L}{\partial \dot{q}_i}\delta q_i\right]_{t_0}^{t_f} - \int_{t_0}^{t_f} \sum_{i=1}^n\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\delta q_i dt.
\end{align*}
At the end times of our path, $\delta q_i$ should be zero, as the endpoints are fixed. Hence the term $\left[\dfrac{\partial L}{\partial \dot{q}_i}\delta q_i\right]_{t_0}^{t_f}$ should evaluate to zero. Our expression then simplifies to:
\begin{align*}
S' &\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \delta q_i dt + \left[\sum_{i=1}^n\dfrac{\partial L}{\partial \dot{q}_i}\delta q_i\right]_{t_0}^{t_f} - \int_{t_0}^{t_f} \sum_{i=1}^n\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\delta q_i dt \\
&\approx \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) + \sum_{i=1}^n\delta q_i \left(\dfrac{\partial L}{\partial q_i} - \dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\right) dt.
\end{align*}
The first of these terms is present in $S$, too. Hence:
\begin{align*}
S' - S &\approx \int_{t_0}^{t_f}\sum_{i=1}^n\delta q_i \left(\dfrac{\partial L}{\partial q_i} - \dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\right) dt.
\end{align*}
We can denote this first-order approximation to how much $S$ changes as $\delta S$. Hence:
\begin{align}
\delta S &= \int_{t_0}^{t_f} \sum_{i=1}^n\delta q_i \left(\dfrac{\partial L}{\partial q_i} - \dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i}\right) dt. \label{deltaSELE}
\end{align}
Equation \eqref{deltaSELE} will be zero irrespective of $\delta q_i$ if and only if Equations \eqref{ELE} hold.
