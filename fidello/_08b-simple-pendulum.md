## Example: simple pendulum
To see the Hamilton-Jacobi equation in action, let us apply it to the simple pendulum. Replacing $p_{\theta}$ with $\dfrac{\partial S}{\partial \theta}$:

\begin{align*}
H &= \dfrac{\left(\dfrac{\partial S}{\partial \theta}\right)^2}{2ml^2}+mgl\sin{\theta}. 
\end{align*}

Substituting into Equation \eqref{HJE} yields:

\begin{align*}
\dfrac{\left(\dfrac{\partial S}{\partial \theta}\right)^2}{2ml^2}+mgl\sin{\theta} + \dfrac{\partial S}{\partial t} &= 0.
\end{align*}

Substituting Equation \eqref{HJE} yields:

\begin{align*}
\dfrac{\left(\dfrac{dW}{d\theta}\right)^2}{2ml^2}+mgl\sin{\theta} -E &= 0.
\end{align*}

Multiplying both sides by $2ml^2$:

\begin{align*}
\left(\dfrac{dW}{d\theta}\right)^2 + 2m^2gl^3\sin{\theta}-2ml^2E &=0.
\end{align*}

Rearranging yields:
\begin{align*}
\left(\dfrac{dW}{d\theta}\right)^2 &= 2ml^2E - 2m^2gl^3\sin{\theta}.
\end{align*}

Square-rooting both sides yields:

\begin{align*}
\dfrac{dW}{d\theta} &= \pm\sqrt{2ml^2E-2m^2gl^3\sin{\theta}} \\
&= \pm \sqrt{2ml^2(E-mgl\sin{\theta})}.
\end{align*}

Integrating yields:

\begin{align*}
W(\theta) &= \pm\sqrt{2m}l \int \sqrt{E-mgl\sin{\theta}} d\theta.
\end{align*}

Next, we assume $\dfrac{\partial S}{\partial E}$ will equal a constant, say $\beta$. Then we differentiate Equation \eqref{SWE}:
\begin{align*}
\dfrac{\partial W}{\partial E} - t &= \beta.
\end{align*}

Rearranging yields:

\begin{align*}
t+\beta &= \dfrac{\partial W}{\partial E}.
\end{align*}

Hence:
\begin{align*}
t+\beta &= \pm \sqrt{2m}l\int \dfrac{1}{2\sqrt{E-mgl\sin{\theta}}} d\theta \\
&= \pm \int \sqrt{\dfrac{ml^2}{2(E-mgl\sin{\theta})}} d\theta.
\end{align*}
The integrand &mdash; given it is integrated with respect to $\theta$ to obtain $t$ (plus some constant) &mdash; must be $\dfrac{dt}{d\theta}$ which is one over $\dfrac{d\theta}{dt}$ or $\dot{\theta}$. Therefore:

\begin{align*}
\dfrac{dt}{d\theta} &= \pm \sqrt{\dfrac{ml^2}{2(E-mgl\sin{\theta})}}.
\end{align*}

Hence $\dot{\theta}$ is:

\begin{align}
\dot{\theta} &= \pm \sqrt{\dfrac{2(E-mgl\sin{\theta})}{ml^2}}.\label{dtheta}
\end{align}

Squaring Equation \eqref{dtheta} yields:
\begin{align}
\dot{\theta}^2 &= \dfrac{2(E-mgl\sin{\theta})}{ml^2}.\label{dtheta2}
\end{align}
Differentiating the left-hand side of Equation \eqref{dtheta2} with respect to $\theta$:
\begin{align}
\dfrac{d\dot{\theta}^2}{d\theta} &= 2\dot{\theta}\dfrac{d\dot{\theta}}{d\theta} \nonumber \\
&= 2\dfrac{d\dot{\theta}}{d\theta}\dfrac{d\theta}{dt} \nonumber\\
&= 2\dfrac{d\dot{\theta}}{dt} \nonumber\\
&= 2\ddot{\theta}.\label{2ddottheta}
\end{align}
Differentiating the right-hand side of Equation \eqref{dtheta2} and equating to the right-hand side of Equation \eqref{2ddottheta} yields:
\begin{align}
2\ddot{\theta} &= -\dfrac{-2mgl\cos{\theta}}{ml^2} \nonumber\\
&= \dfrac{-2g\cos{\theta}}{l}.\nonumber
\end{align}
Dividing both sides by $2$ yields:
\begin{align}
\ddot{\theta} &= -\dfrac{g}{l}\cos{\theta}.\nonumber
\end{align}
Which we can see matches Equation \eqref{simppen}.