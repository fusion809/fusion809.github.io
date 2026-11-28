#### Runge-Kutta-Fehlberg 4/5th-order method
An arguably better method for numerically solving ordinary differential equations is the Runge-Kutta-Fehlberg fourth/fifth-order method (RKF45), as it adapts our step size $h$ to ensure that the error in our solution never exceeds a certain threshold. In this method, we use:

\begin{align*}
k_1 &= hf(x_i, y_i) \\
k_2 &= hf\left(x_i + \dfrac{h}{4}, y_i + \dfrac{k_1}{4}\right) \\
k_3 &= hf\left(x_i + \dfrac{3h}{8}, y_i + \dfrac{3k_1}{32} + \dfrac{9k_2}{32}\right) \\
k_4 &= hf\left(x_i + \dfrac{12h}{13}, y_i + \dfrac{1932k_1}{2197} - \dfrac{7200k_2}{2197} + \dfrac{7296k_3}{2197}\right) \\
k_5 &= hf\left(x_i + h, y_i + \dfrac{439k_1}{216} - 8k_2 + \dfrac{3680k_3}{513} - \dfrac{845k_4}{4104}\right) \\
k_6 &= hf\left(x_i + \dfrac{h}{2}, y_i - \dfrac{8k_1}{27} + 2k_2 - \dfrac{3544k_3}{2565} + \dfrac{1859k_4}{4104} - \dfrac{11k_5}{40}\right) \\
z_{i+1} &= y_i + \dfrac{16k_1}{135} + \dfrac{6656k_3}{12825} + \dfrac{28561k_4}{56430} - \dfrac{9k_5}{50} + \dfrac{2k_6}{55} \\
y_{i+1} &= y_i + \dfrac{25k_1}{216} + \dfrac{1408k_3}{2565} + \dfrac{2197k_4}{4104} - \dfrac{k_5}{5} \\
\mathrm{TE} &= |y_{i+1} - z_{i+1}| \\
s &= \left(\dfrac{\epsilon}{2\mathrm{TE}}\right)^{1/5}.
\end{align*}
Here $\epsilon$ is the maximum absolute error in our approximation we are willing to accept. $\mathrm{TE}$ is our approximation of the truncation error as our estimate for $y_{i+1}$ is fourth-order and $z_{i+1}$ is the fifth-order estimate. If $\mathrm{TE} \geq \epsilon$, then we rescale $h$ using $h=h\cdot s$ and repeat the approximation again. If not, then we accept our approximation to the next $y$ value and rescale $h$ using $h=h\cdot s$. This method is very efficient and is the numerical technique I use for integrating ordinary differential equations in all of my JavaScript-based ordinary differential equations solver webpages. 

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
How would we numerically integrate this system? We could apply RKF45 to determine $y_{i+1}$ (where $y_0$ and $y'_0$ will be our initial conditions) using this scheme:
\begin{align*}
k_1 &= hu_i\\
l_1 &= -hy_i \\
k_2 &= h\left(u_i+\dfrac{l_1}{4}\right) \\
l_2 &= -h\left(y_i+\dfrac{k_1}{4}\right) \\
k_3 &= h\left(u_i+\dfrac{3l_1}{32}+\dfrac{9l_2}{32}\right) \\
l_3 &= -h\left(y_i+\dfrac{3k_1}{32}+\dfrac{9k_2}{32}\right) \\
k_4 &= h\left(u_i+\dfrac{1932l_1}{2197}-\dfrac{7200l_2}{2197}+\dfrac{7296l_3}{2197}\right) \\
l_4 &= -h\left(y_i+\dfrac{1932k_1}{2197}-\dfrac{7200k_2}{2197}+\dfrac{7296k_3}{2197}\right) \\
k_5 &= h\left(u_i+\dfrac{439l_1}{216}-8l_2+\dfrac{3680l_3}{513}-\dfrac{845l_4}{4104}\right) \\
l_5 &= -h\left(y_i+\dfrac{439k_1}{216}-8k_2+\dfrac{3680k_3}{513}-\dfrac{845k_4}{4104}\right) \\
k_6 &= h\left(u_i-\dfrac{8l_1}{27}+2l_2-\dfrac{3544l_3}{2565}+\dfrac{1859l_4}{4104}-\dfrac{11l_5}{40}\right) \\
l_6 &= -h\left(y_i-\dfrac{8k_1}{27}+2k_2-\dfrac{3544k_3}{2565}+\dfrac{1859k_4}{4104}-\dfrac{11k_5}{40}\right) \\
u_{4,i+1} &= u_{i} + \dfrac{25l_1}{216} + \dfrac{1408l_3}{2565} + \dfrac{2197l_4}{4104} - \dfrac{l_5}{5}\\
u_{5,i+1} &= u_{i} + \dfrac{16l_1}{135} + \dfrac{6656l_3}{12825} + \dfrac{28561l_4}{56430} - \dfrac{9l_5}{50} + \dfrac{2l_6}{55}\\
y_{4,i+1} &= y_{i} + \dfrac{25k_1}{216} + \dfrac{1408k_3}{2565} + \dfrac{2197k_4}{4104} - \dfrac{k_5}{5}\\
y_{5,i+1} &= y_{i} + \dfrac{16k_1}{135} + \dfrac{6656k_3}{12825} + \dfrac{28561k_4}{56430} - \dfrac{9k_5}{50} + \dfrac{2k_6}{55}\\
\mathrm{TE} &= \max_{f \in \{u, y\}} \Big|f_{5,i+1}-f_{4,i+1}\Big|\\
s &= \left(\dfrac{\epsilon}{2\mathrm{TE}}\right)^{1/5}.
\end{align*}
I should clarify that we are maximizing over the error estimates in our two functions $u$ and $y$ to calculate $\mathrm{TE}$. So whichever error estimate is bigger out of those for $u$ and $y$, we choose for $\mathrm{TE}$. If $\mathrm{TE} < \epsilon$, we continue onto the next step with our new $y_i$ being the old $y_{5,i+1}$ and the new $u_i$ being the old $u_{5,i+1}$ and obtain our new $h$ by multiplying it by $s$. If $\mathrm{TE} \geq \epsilon$, we multiply $h$ by $s$ and repeat the above calculations again.

So we will use the initial conditions $y(0)=0$ and $y'(0)=1$. This gives the solution $y(x)=\sin{x}$ and integrate our ordinary differential equation on the interval $x\in [0, 2\pi]$.

\codelisting{RKF45coup}{Runge-Kutta-Fehlberg 4/5th-order method solution to $y''=-y$ with $y(0)=0$ and $y'(0)=1$}{./code/RKF45coup}{fidello/code/RKF45coup.jl}
\output{./code/RKF45coup}

In Code Listing \coderef{RKF45coup}, we used $\epsilon=10^{-8}$ as our error tolerance. You may notice that our maximum errors in the output table do indeed exceed this. This is possible as our truncation error estimate is solely for that step, it is not the cumulative error for all steps up until the point we are in our calculation. 