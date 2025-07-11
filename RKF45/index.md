@def title="Runge-Kutta-Fehlberg 4/5th-order method"

The **Runge-Kutta-Fehlberg 4th-order method with 5th-order error checking** (RKF45) is a numerical integration technique for ordinary differential equations (ODEs) that utilizes an adaptive step size to achieve a prescribed error tolerance. It is used by each page in this website that solves an ODE system. In it, we assume our ODE system has the form

\begin{align*}
\dfrac{d\mathbf{x}}{dt} &= f(\mathbf{x}, t). \\
\end{align*}

To numerically integrate this system, we discretize our domain of $t$ values with an adaptive step size of $h$ and at each step $i$, we utilize the integration scheme

\begin{align*}
\mathbf{k_1} &= hf(\mathbf{x}_i, t_i) \\
\mathbf{k_2} &= hf\left(\mathbf{x}_i + \dfrac{\mathbf{k_1}}{4}, t_i + \dfrac{h}{4}\right) \\
\mathbf{k_3} &= hf\left(\mathbf{x}_i + \dfrac{3\mathbf{k_1}}{32} + \dfrac{9\mathbf{k_2}}{32}, t_i + \dfrac{3h}{8}\right) \\
\mathbf{k_4} &= hf\left(\mathbf{x}_i + \dfrac{1932\mathbf{k_1}}{2197} - \dfrac{7200\mathbf{k_2}}{2197} + \dfrac{7296\mathbf{k_3}}{2197}, t_i + \dfrac{12h}{13}\right) \\
\mathbf{k_5} &= hf\left(\mathbf{x}_i + \dfrac{439\mathbf{k_1}}{216} - 8\mathbf{k_2} + \dfrac{3680\mathbf{k_3}}{513} - \dfrac{845\mathbf{k_4}}{4104}, t_i + h\right) \\
\mathbf{k_6} &= hf\left(\mathbf{x}_i - \dfrac{8\mathbf{k_1}}{27}  + 2\mathbf{k_2} - \dfrac{3544\mathbf{k_3}}{2565} + \dfrac{1859\mathbf{k_4}}{4104} - \dfrac{11\mathbf{k_4}}{40}, t_i + \dfrac{h}{2}\right) \\
\mathbf{X}_{1, i} &= \mathbf{x}_i + \dfrac{25\mathbf{k_1}}{216} + \dfrac{1408\mathbf{k_3}}{2565} + \dfrac{2197\mathbf{k_4}}{4104} - \dfrac{\mathbf{k_5}}{5} \\
\mathbf{X}_{2, i} &= \mathbf{x}_i + \dfrac{16\mathbf{k_1}}{135} + \dfrac{6656\mathbf{k_3}}{12825} + \dfrac{28561\mathbf{k_4}}{56430} - \dfrac{9\mathbf{k_5}}{50}+\dfrac{2\mathbf{k_6}}{55} \\
\mathbf{R}_i &= \dfrac{|\mathbf{X}_{1,i} - \mathbf{X}_{2,i}|}{h} \\
\mathbf{S}_i &= \left[\dfrac{\epsilon}{2\mathbf{R}_i}\right]^{1/4}\\
r &= \max{\mathbf{R}_i} \\
s &= \min{\mathbf{S}_i} \\
h &= sh.
\end{align*}

If $r \leq \epsilon $, we let $\mathbf{X}_{1,i}$ be our value of $\mathbf{x}_{i+1}$ and proceed to the next step. Otherwise, we repeat the calculation with our updated step size. 

In the $\mathbf{R}_i$ calculation, $|\mathbf{v}|$ denotes a vector whose elements are the absolute value of the elements of $\mathbf{v}$, not its vector magnitude. Here $r$ is the maximum error estimate we have for step $i$.

