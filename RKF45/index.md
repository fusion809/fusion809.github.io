@def title="Runge-Kutta-Fehlberg 4/5th-order method"

The **Runge-Kutta-Fehlberg 4th-order method with 5th-order error checking** (RKF45) is a numerical integration technique for [ordinary differential equations](https://en.wikipedia.org/wiki/Ordinary_differential_equation) (ODEs) that utilizes an adaptive step size to achieve a prescribed error tolerance, $\epsilon$. It is used by each page in this website that solves an ODE system. In it, we assume our ODE system has the form

\begin{align*}
\dfrac{d\mathbf{x}}{dt} &= f(\mathbf{x}, t). \\
\end{align*}

To numerically integrate this system, we discretize our domain of $t$ values. We will not know how many points we discretize our domain into until after we have finished the integration, due to our adaptive step size. At each step $i$, we will use the step size $h_i$ and we will utilize the integration scheme

\begin{align*}
\mathbf{k}_{1,i} &= h_if(\mathbf{x}_i, t_i) \\
\mathbf{k}_{2,i} &= h_if\left(\mathbf{x}_i + \dfrac{\mathbf{k}_{1,i}}{4}, t_i + \dfrac{h}{4}\right) \\
\mathbf{k}_{3,i} &= h_if\left(\mathbf{x}_i + \dfrac{3\mathbf{k}_{1,i}}{32} + \dfrac{9\mathbf{k}_{2,i}}{32}, t_i + \dfrac{3h}{8}\right) \\
\mathbf{k}_{4,i} &= h_if\left(\mathbf{x}_i + \dfrac{1932\mathbf{k}_{1,i}}{2197} - \dfrac{7200\mathbf{k}_{2,i}}{2197} + \dfrac{7296\mathbf{k}_{3,i}}{2197}, t_i + \dfrac{12h}{13}\right) \\
\mathbf{k}_{5,i} &= h_if\left(\mathbf{x}_i + \dfrac{439\mathbf{k}_{1,i}}{216} - 8\mathbf{k}_{2,i} + \dfrac{3680\mathbf{k}_{3,i}}{513} - \dfrac{845\mathbf{k}_{4,i}}{4104}, t_i + h\right) \\
\mathbf{k}_{6,i} &= h_if\left(\mathbf{x}_i - \dfrac{8\mathbf{k}_{1,i}}{27}  + 2\mathbf{k}_{2,i} - \dfrac{3544\mathbf{k}_{3,i}}{2565} + \dfrac{1859\mathbf{k}_{4,i}}{4104} - \dfrac{11\mathbf{k}_{4,i}}{40}, t_i + \dfrac{h}{2}\right) \\
\mathbf{X}_{1, i} &= \mathbf{x}_i + \dfrac{25\mathbf{k}_{1,i}}{216} + \dfrac{1408\mathbf{k}_{3,i}}{2565} + \dfrac{2197\mathbf{k}_{4,i}}{4104} - \dfrac{\mathbf{k}_{5,i}}{5} \\
\mathbf{X}_{2, i} &= \mathbf{x}_i + \dfrac{16\mathbf{k}_{1,i}}{135} + \dfrac{6656\mathbf{k}_{3,i}}{12825} + \dfrac{28561\mathbf{k}_{4,i}}{56430} - \dfrac{9\mathbf{k}_{5,i}}{50}+\dfrac{2\mathbf{k}_{6,i}}{55} \\
\mathbf{R}_i &= \dfrac{|\mathbf{X}_{1,i} - \mathbf{X}_{2,i}|}{h} \\
\mathbf{S}_i &= \left[\dfrac{\epsilon}{2\mathbf{R}_i}\right]^{1/4}\\
r_i &= \max{\mathbf{R}_i} \\
s_i &= \min{\mathbf{S}_i} \\
h_i &= s_ih_i.
\end{align*}

If $r \leq \epsilon $, we let $\mathbf{X}_{1,i}$ be our value of $\mathbf{x}_{i+1}$, $h_{i+1}$ be $h_i$, increment $i$ by 1, and proceed to the next step. Otherwise, we repeat the calculation with our updated step size. In the $\mathbf{R}_i$ calculation, $|\mathbf{v}|$ denotes a vector whose elements are the absolute value of the elements of $\mathbf{v}$, not its vector magnitude. Here $r$ is the maximum error estimate we have for step $i$.

$h$ can become unmanageably small if the problem one is solving is too numerically unstable. This can be seen in the [double elastic pendulum solver](/doubleElasticPendulum/) webpage if you reduce $\epsilon$ to too small a value or increase $t_f$ too much. 
