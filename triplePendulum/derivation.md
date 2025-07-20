@def hassim=false;
@def title="Derivation of the equations of motion of the triple pendulum"
@def mintoclevel=1
@def maxtoclevel=5

In this article, the equations of motion of the triple pendulum will be derived via the Euler-Lagrange equations with dissipation. See [this article](/triplePendulum/) for a solver of this system of equations. 

\tableofcontents

## Positions, velocities and generalized basis vectors
~~~
<figure>
    <img src="/triplePendulum/Triple pendulum.svg" width="500px"></img>
    <figcaption style="font-weight: bold; font-size: 18px;">Figure 1: Diagram of the triple pendulum.</b></figcaption>
</figure>
~~~

### Pendulum bob 1
\begin{align*}
& x_{1b} &= l_1 \cos{\theta_1}, & y_{1b} &= l_1 \sin{\theta_1}, & \dot{x}_{1b} &= -l_1 \dot{\theta}_1 \sin{\theta_1}, & \dot{y}_{1b} &= l_1 \dot{\theta}_1 \cos{\theta_1}\\
& \therefore v_{1b} &= l_1 \dot{\theta}_1 \\
& \vec{v}_{1b} &= l_1 \dot{\theta}_1\begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{1b, \theta_1} &= l_1 \begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{1b, \theta_2} &= \hat{e}_{1b, \theta_3} &= \vec{0} \\
\end{align*}
### Pendulum rod 1
\begin{align*}
& x_{1r} &= \dfrac{l_1 \cos{\theta_1}}{2}, & y_{1r} &= \dfrac{l_1 \sin{\theta_1}}{2}, & \dot{x}_{1r} &= -\dfrac{l_1 \dot{\theta}_1 \cos{\theta_1}}{2}, & \dot{y}_{1r} &= \dfrac{l_1\dot{\theta}_1 \cos{\theta_1}}{2} \\
& \therefore v_{1r} &= \dfrac{l_1 \dot{\theta}_1}{2} \\
& \vec{v}_{1r} &= \dfrac{l_1 \dot{\theta}_1}{2}\begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{1r, \theta_1} &= \dfrac{l_1}{2} \begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{1r, \theta_2} &= \hat{e}_{1r, \theta_3} &= \vec{0}
\end{align*}
### Pendulum bob 2
Defining $\Delta_{ij} = \theta_i - \theta_j$, we get
\begin{align*}
& x_{2b} &= l_1\cos{\theta_1} + l_2 \cos{\theta_2}, & y_{2b} &= l_1 \sin{\theta_1} + l_2 \sin{\theta_2} & \dot{x}_{2b} &= -l_1 \dot{\theta_1}\sin{\theta_1} - l_2\dot{\theta_2}\sin{\theta_2}, &\dot{y}_{2b} &= l_1 \dot{\theta_1} \cos{\theta_1} + l_2 \dot{\theta_2}\cos{\theta_2} \\
& \therefore v_{2b} &= \sqrt{l_1^2 \dot{\theta}_1^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + l_2^2 \dot{\theta}_2^2} \\
& \vec{v}_{2b} &= \begin{bmatrix}
-l_1 \dot{\theta_1}\sin{\theta_1} - l_2\dot{\theta_2}\sin{\theta_2} \\
l_1 \dot{\theta_1} \cos{\theta_1} + l_2 \dot{\theta_2}\cos{\theta_2}
\end{bmatrix}, \hat{e}_{2b, \theta_1} &= l_1 \begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{2b, \theta_2} &= l_2 \begin{bmatrix}
-\sin{\theta_2} \\
\cos{\theta_2}
\end{bmatrix}, & \hat{e}_{2b, \theta_3} &= \vec{0}
\end{align*}
### Pendulum rod 2
\begin{align*}
& x_{2r} &= l_1\cos{\theta_1} + \dfrac{l_2 \cos{\theta_2}}{2}, & y_{2r} &= l_1 \sin{\theta_1} + \dfrac{l_2 \sin{\theta_2}}{2} & \dot{x}_{2r} &= -l_1 \dot{\theta_1}\sin{\theta_1} - \dfrac{l_2\dot{\theta_2}\sin{\theta_2}}{2}, &\dot{y}_{2r} &= l_1 \dot{\theta_1} \cos{\theta_1} + \dfrac{l_2 \dot{\theta_2}\cos{\theta_2}}{2} \\
& \therefore v_{2r} &= \sqrt{l_1^2 \dot{\theta}_1^2 + l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + \dfrac{l_2^2 \dot{\theta}_2^2}{4}} \\
& \vec{v}_{2r} &= \begin{bmatrix}
-l_1 \dot{\theta_1}\sin{\theta_1} - \dfrac{l_2\dot{\theta_2}\sin{\theta_2}}{2} \\
l_1 \dot{\theta_1} \cos{\theta_1} + \dfrac{l_2 \dot{\theta_2}\cos{\theta_2}}{2}
\end{bmatrix}, \hat{e}_{2r, \theta_1} &= l_1 \begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{2r, \theta_2} &= \dfrac{l_2}{2} \begin{bmatrix}
-\sin{\theta_2} \\
\cos{\theta_2}
\end{bmatrix}, & \hat{e}_{2r, \theta_3} &= \vec{0}
\end{align*}
### Pendulum bob 3
\begin{align*}
& x_{3b} &= l_1\cos{\theta_1} + l_2 \cos{\theta_2} + l_3 \cos{\theta_3}, & y_{3b} &= l_1 \sin{\theta_1} + l_2\sin{\theta_2} + l_3 \sin{\theta_3} & \dot{x}_{3b} &= -l_1 \dot{\theta_1}\sin{\theta_1} - l_2 \dot{\theta}_2\sin{\theta_2} - l_3\dot{\theta_3}\sin{\theta_3}, &\dot{y}_{3b} &= l_1 \dot{\theta_1} \cos{\theta_1} + l_2 \dot{\theta}_2 \cos{\theta_2} + l_3 \dot{\theta_3}\cos{\theta_3}\\
& \therefore v_{3b} &= \sqrt{l_1^2 \dot{\theta}_1^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + l_2^2 \dot{\theta}_2^2 + 2l_2 l_3 \dot{\theta}_2 \dot{\theta}_3 \cos{\Delta_{32}} + 2l_1l_3 \dot{\theta}_1\dot{\theta}_3 \cos{\Delta_{31}} + l_3^2 \dot{\theta}_3^2}\\
& \vec{v}_{3b} &= \begin{bmatrix}
-l_1 \dot{\theta_1}\sin{\theta_1} - l_2\dot{\theta_2}\sin{\theta_2} - l_3\dot{\theta_3}\sin{\theta_3} \\
l_1 \dot{\theta_1} \cos{\theta_1} + l_2 \dot{\theta_2}\cos{\theta_2} + l_3 \dot{\theta_3}\cos{\theta_3}
\end{bmatrix}, \hat{e}_{3b, \theta_1} &= l_1 \begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{3b, \theta_2} &= l_2 \begin{bmatrix}
-\sin{\theta_2} \\
\cos{\theta_2}
\end{bmatrix}, & \hat{e}_{3b, \theta_3} &= l_3 \begin{bmatrix}
-\sin{\theta_3}
\cos{\theta_3}
\end{bmatrix}
\end{align*}
### Pendulum rod 3
\begin{align*}
& x_{3r} &= l_1\cos{\theta_1} + l_2 \cos{\theta_2} + \dfrac{l_3 \cos{\theta_3}}{2}, & y_{3r} &= l_1 \sin{\theta_1} + l_2\sin{\theta_2} + \dfrac{l_3 \sin{\theta_3}}{2} & \dot{x}_{3r} &= -l_1 \dot{\theta_1}\sin{\theta_1} - l_2 \dot{\theta}_2\sin{\theta_2} - \dfrac{l_3\dot{\theta_3}\sin{\theta_3}}{2}, &\dot{y}_{3r} &= l_1 \dot{\theta_1} \cos{\theta_1} + l_2 \dot{\theta}_2 \cos{\theta_2} + \dfrac{l_3 \dot{\theta_3}\cos{\theta_3}}{2}\\
& \therefore v_{3r} &= \sqrt{l_1^2 \dot{\theta}_1^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + l_2^2 \dot{\theta}_2^2 + l_1l_3 \dot{\theta}_1\dot{\theta}_3\cos{\Delta_{31}} + l_2l_3\dot{\theta}_2 \dot{\theta}_3 \cos{\Delta_{32}} + \dfrac{l_3^2 \dot{\theta}_3^2}{4}}\\
& \vec{v}_{3r} &= \begin{bmatrix}
-l_1 \dot{\theta_1}\sin{\theta_1} - \dfrac{l_2\dot{\theta_2}\sin{\theta_2}}{2}  - \dfrac{l_3\dot{\theta_3}\sin{\theta_3}}{2}\\
l_1 \dot{\theta_1} \cos{\theta_1} + \dfrac{l_2 \dot{\theta_2}\cos{\theta_2}}{2} + \dfrac{l_3 \dot{\theta_3}\cos{\theta_3}}{2}
\end{bmatrix}, \hat{e}_{3r, \theta_1} &= l_1 \begin{bmatrix}
-\sin{\theta_1} \\
\cos{\theta_1}
\end{bmatrix}, & \hat{e}_{3r, \theta_2} &= l_2 \begin{bmatrix}
-\sin{\theta_2} \\
\cos{\theta_2}
\end{bmatrix}, & \hat{e}_{3r, \theta_3} &= \dfrac{l_3}{2} \begin{bmatrix}
-\sin{\theta_3} \\
\cos{\theta_3}
\end{bmatrix}
\end{align*}

## Kinetic energy
\begin{align*}
T &= \dfrac{m_{1b}}{2} v_{1b}^2 + \dfrac{m_{1r}}{2} v_{1r}^2 + \dfrac{m_{1r}I_{\mathrm{cm},1r} \omega_{1r}^2}{2} + \dfrac{m_{2b}}{2} v_{2b}^2 + \dfrac{m_{2r}}{2} v_{2r}^2 + \dfrac{m_{2r}I_{\mathrm{cm},2r} \omega_{2r}^2}{2} + \dfrac{m_{3b}}{2} v_{3b}^2 + \dfrac{m_{3r}}{2} v_{3r}^2 + \dfrac{m_{3r}I_{\mathrm{cm},3r} \omega_{3r}^2}{2} \\
&= \dfrac{m_{1b}l_1^2 \dot{\theta}_1^2}{2} + \dfrac{m_{1r}l_1^2 \dot{\theta}_1^2}{8} + \dfrac{m_{1r}l_1^2 \dot{\theta}_1^2}{24} + \dfrac{m_{2b}}{2}(l_1^2 \dot{\theta}_1^2 + 2l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \cos{\Delta_{21}} + l_2^2 \dot{\theta}_2^2) + \dfrac{m_{2r}}{2}\left(l_1^2 \dot{\theta}_1^2 + l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \cos{\Delta_{21}} + \dfrac{l_2^2 \dot{\theta}_2^2}{4}\right)  + \dfrac{m_{2r}l_2^2 \dot{\theta}_2^2}{24} + \dfrac{m_{3b}}{2} \left(l_1^2 \dot{\theta}_1^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + l_2^2 \dot{\theta}_2^2 + 2l_2 l_3 \dot{\theta}_2 \dot{\theta}_3 \cos{\Delta_{32}} + 2l_1l_3 \dot{\theta}_1\dot{\theta}_3 \cos{\Delta_{31}} + l_3^2 \dot{\theta}_3^2\right)+ \dfrac{m_{3r}}{2}\left(l_1^2 \dot{\theta}_1^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + l_2^2 \dot{\theta}_2^2 + l_1l_3 \dot{\theta}_1\dot{\theta}_3\cos{\Delta_{31}} + l_2l_3\dot{\theta}_2 \dot{\theta}_3 \cos{\Delta_{32}} + \dfrac{l_3^2 \dot{\theta}_3^2}{4}\right)  + \dfrac{m_{3r}l_3^2 \dot{\theta}_3^2}{24} \\
&= \dfrac{1}{2}\left(m_{1b}+\dfrac{m_{1r}}{3} + m_{2b} + m_{2r} + m_{3b} + m_{3r}\right)l_1^2 \dot{\theta}_1^2 + \dfrac{1}{2}\left(m_{2b} + \dfrac{m_{2r}}{3} + m_{3b} + m_{3r}\right)l_2^2 \dot{\theta}_2^2 + \dfrac{1}{2}\left(m_{3b} + \dfrac{m_{3r}}{3}\right)l_3^2 \dot{\theta}_3^2 + \left(m_{2b} + \dfrac{m_{2r}}{2} + m_{3b} + m_{3r}\right)l_1l_2\dot{\theta}_1\dot{\theta}_2\cos{\Delta_{21}} + \left(m_{3b} + \dfrac{m_{3r}}{2}\right)l_3(l_2\dot{\theta}_2\dot{\theta}_3\cos{\Delta_{32}}+l_1\dot{\theta}_1\dot{\theta}_3\cos{\Delta_{31}})
\end{align*}

Defining $M_1 = m_{1b}+\dfrac{m_{1r}}{3} + m_{2b} + m_{2r} + m_{3b} + m_{3r}$, $M_2 = m_{2b} + \dfrac{m_{2r}}{3} + m_{3b} + m_{3r}$, $M_3 = m_{3b} + \dfrac{m_{3r}}{3}$, $\mu_1 = m_{1b}+\dfrac{m_{1r}}{2} + m_{2b} + m_{2r} + m_{3b} + m_{3r}$, $\mu_2 = m_{2b} + \dfrac{m_{2r}}{2} + m_{3b} + m_{3r}$ and $\mu_3 = m_{3b} + \dfrac{m_{3r}}{2}$.

\begin{align*}
T &= \dfrac{M_1 l_1^2 \dot{\theta}_1^2}{2} + \dfrac{M_2 l_2^2 \dot{\theta}_2^2}{2} + \dfrac{M_3 l_3^2 \dot{\theta}_3^2}{2} + \mu_2 l_1l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + \mu_3 l_3\dot{\theta}_3(l_2\dot{\theta}_2\cos{\Delta_{32}}+l_1\dot{\theta}_1\cos{\Delta_{31}}).
\end{align*}

## Potential energy
\begin{align*}
V &= m_{1b} gy_{1b} + m_{1r}gy_{1r} + m_{2b}gy_{2b} + m_{2r}gy_{2r} + m_{3b}gy_{3b} + m_{3r}gy_{3r} \\
&= m_{1b}gl_1 \sin{\theta_1} + \dfrac{m_{1r}gl_1\sin{\theta_1}}{2} + m_{2b} g(l_1\sin{\theta_1} + l_2\sin{\theta_2}) + m_{2r}g\left(l_1\sin{\theta_1}+\dfrac{l_2\sin{\theta_2}}{2}\right) + m_{3b} g(l_1\sin{\theta_1} + l_2\sin{\theta_2}+l_3\sin{\theta_3}) + m_{3r}g\left(l_1\sin{\theta_1} + l_2\sin{\theta_2} +\dfrac{l_3\sin{\theta_3}}{2}\right) \\
&= \mu_1 gl_1 \sin{\theta_1} + \mu_2 gl_2\sin{\theta_2} + \mu_3 gl_3\sin{\theta_3}.
\end{align*}

## Lagrangian
\begin{align*}
\mathcal{L} &= T - V\\
&= \dfrac{M_1 l_1^2 \dot{\theta}_1^2}{2} + \dfrac{M_2 l_2^2 \dot{\theta}_2^2}{2} + \dfrac{M_3 l_3^2 \dot{\theta}_3^2}{2} + \mu_2 l_1l_2 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}} + \mu_3 l_3\dot{\theta}_3(l_2\dot{\theta}_2\cos{\Delta_{32}}+l_1\dot{\theta}_1\cos{\Delta_{31}}) - \mu_1 gl_1 \sin{\theta_1} - \mu_2 gl_2\sin{\theta_2} - \mu_3 gl_3\sin{\theta_3} \\
&= \dfrac{M_1 l_1^2 \dot{\theta}_1^2}{2} + \dfrac{M_2 l_2^2 \dot{\theta}_2^2}{2} + \dfrac{M_3 l_3^2 \dot{\theta}_3^2}{2} + \mu_2 l_2(l_1 \dot{\theta}_1 \dot{\theta}_2 \cos{\Delta_{21}}-g\sin{\theta_2}) + \mu_3 l_3(\dot{\theta}_3(l_2\dot{\theta}_2\cos{\Delta_{32}}+l_1\dot{\theta}_1\cos{\Delta_{31}})-g\sin{\theta_3}) - \mu_1 gl_1 \sin{\theta_1}.
\end{align*}

## Generalized dissipative force
### $\theta_1$
\begin{align*}
Q_{\theta_1} &= -(b_{1b}+c_{1b}|v_{1b}|)\vec{v}_{1b} \cdot \hat{e}_{1b, \theta_1}-(b_{1r}+c_{1r}|v_{1r}|)\vec{v}_{1r} \cdot \hat{e}_{1r, \theta_1} -(b_{2b}+c_{2b}|v_{2b}|)\vec{v}_{2b} \cdot \hat{e}_{2b, \theta_1}-(b_{2r}+c_{2r}|v_{2r}|)\vec{v}_{2r} \cdot \hat{e}_{2r, \theta_1} -(b_{3b}+c_{3b}|v_{3b}|)\vec{v}_{3b} \cdot \hat{e}_{3b, \theta_1}-(b_{3r}+c_{3r}|v_{3r}|)\vec{v}_{3r} \cdot \hat{e}_{3r, \theta_1}.
\end{align*}

We will not substitute our values of $v_{2b}$ to $v_{3r}$ as they will only complicate our equation
\begin{align*}
Q_{\theta_1} &=-(b_{1b}+c_{1b}l_1|\dot{\theta}_1|)l_1^2 \dot{\theta}_1-\left(b_{1r}+c_{1r}\dfrac{l_1|\dot{\theta}_1|}{2}\right)\dfrac{l_1^2 \dot{\theta}_1}{4} -(b_{2b}+c_{2b}|v_{2b}|)(l_1^2\dot{\theta}_1 + l_1l_2\dot{\theta}_2\cos{\Delta_{21}})-(b_{2r}+c_{2r}|v_{2r}|)(l_1^2\dot{\theta}_1 + \dfrac{l_1l_2\dot{\theta}_2\cos{\Delta_{21}}}{2}) -(b_{3b}+c_{3b}|v_{3b}|)(l_1^2\dot{\theta}_1 + l_1l_2\dot{\theta}_2\cos{\Delta_{21}}+l_1l_3\dot{\theta}_3\cos{\Delta_{31}})-(b_{3r}+c_{3r}|v_{3r}|)(l_1^2\dot{\theta}_1 + l_1l_2\dot{\theta}_2\cos{\Delta_{21}}+\dfrac{l_1l_3\dot{\theta}_3\cos{\Delta_{31}}}{2}).
\end{align*}

### $\theta_2$
The generalized dissipative force for $\theta_2$ is (pendulum 1 terms are ignored because their generalized basis vectors are zero)

\begin{align*}
Q_{\theta_2} &= -(b_{2b}+c_{2b}|v_{2b}|)\vec{v}_{2b} \cdot \hat{e}_{2b, \theta_2} - (b_{2r}+c_{2r}|v_{2r}|) \vec{v}_{2r} \cdot \hat{e}_{2r, \theta_2} -(b_{3b}+c_{3b}|v_{3b}|)\vec{v}_{3b} \cdot \hat{e}_{3b, \theta_2} - (b_{3r}+c_{3r}|v_{3r}|) \vec{v}_{3r} \cdot \hat{e}_{3r, \theta_2}\\
&= -(b_{2b}+c_{2b}|v_{2b}|)(l_2^2 \dot{\theta}_2 + l_1l_2 \dot{\theta}_1 \cos{\Delta_{21}}) - (b_{2r}+c_{2r}|v_{2r}|) (\dfrac{l_2^2 \dot{\theta}_2}{4} + \dfrac{l_1l_2 \dot{\theta}_1 \cos{\Delta_{21}}}{2}) -(b_{3b}+c_{3b}|v_{3b}|)(l_2^2\dot{\theta}_2 + l_1l_2 \dot{\theta}_1\cos{\Delta_{21}} + l_2l_3\dot{\theta}_3 \cos{\Delta_{32}}) - (b_{3r}+c_{3r}|v_{3r}|) (l_2^2\dot{\theta}_2 + l_1l_2 \dot{\theta}_1\cos{\Delta_{21}} + \dfrac{l_2l_3\dot{\theta}_3 \cos{\Delta_{32}}}{2}).
\end{align*}

### $\theta_3$
\begin{align*}
Q_{\theta_3} &= -(b_{3b}+c_{3b}|v_{3b}|)\vec{v}_{3b} \cdot \hat{e}_{3b, \theta_3} - (b_{3r}+c_{3r}|v_{3r}|) \vec{v}_{3r} \cdot \hat{e}_{3r, \theta_3}\\
&= -(b_{3b}+c_{3b}|v_{3b}|)(l_3^2\dot{\theta}_3 + l_1l_3 \dot{\theta}_1\cos{\Delta_{31}} + l_2l_3\dot{\theta}_2 \cos{\Delta_{32}}) - (b_{3r}+c_{3r}|v_{3r}|) \left(\dfrac{l_3^2\dot{\theta}_3}{4} + \dfrac{l_1l_3 \dot{\theta}_1\cos{\Delta_{31}}}{2} + \dfrac{l_2l_3\dot{\theta}_3 \cos{\Delta_{32}}}{2}\right).
\end{align*}

## Left-hand side of the Euler-Lagrange equations
### $\theta_1$
\begin{align*}
p_{\theta_1} &= \dfrac{\partial \mathcal{L}}{\partial \dot{\theta}_1} \\
&= M_1 l_1^2 \dot{\theta}_1 + \mu_2 l_1l_2 \dot{\theta}_2\cos{\Delta_{21}} + \mu_3 l_1 l_3\dot{\theta}_3\cos{\Delta_{31}} \\
\dot{p}_{\theta_1} &=  M_1 l_1^2 \ddot{\theta}_1 + \mu_2 l_1l_2 (\ddot{\theta}_2\cos{\Delta_{21}} - \dot{\theta}_2(\dot{\theta}_2-\theta_1)\sin{\Delta_{21}}) + \mu_3 l_1 l_3(\ddot{\theta}_3\cos{\Delta_{31}} - \dot{\theta}_3(\dot{\theta}_3-\dot{\theta_1})\sin{\Delta_{31}}) \\
F_{\theta_1} &= \dfrac{\partial \mathcal{L}}{\partial \theta_1} \\
&= -\mu_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2\dfrac{\partial \Delta_{21}}{\partial \theta_1}\sin{\Delta_{21}} - \mu_3l_1 l_3\dot{\theta}_1\dot{\theta}_3\dfrac{\partial \Delta_{31}}{\partial \theta_1}\sin{\Delta_{31}} - \mu_1 gl_1\cos{\theta_1} \\
&= \mu_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2\sin{\Delta_{21}} + \mu_3l_1 l_3\dot{\theta}_1\dot{\theta}_3\sin{\Delta_{31}} - \mu_1 gl_1\cos{\theta_1}
\end{align*}

\begin{align*}
\delta'_{\theta_1} \mathcal{L} &= \dot{p}_{\theta_1} - F_{\theta_1} \\
&= M_1 l_1^2 \ddot{\theta}_1 + \mu_2 l_1l_2 (\ddot{\theta}_2\cos{\Delta_{21}} - \dot{\theta}_2(\dot{\theta}_2-\theta_1)\sin{\Delta_{21}}) + \mu_3 l_1 l_3(\ddot{\theta}_3\cos{\Delta_{31}} - \dot{\theta}_3(\dot{\theta}_3-\dot{\theta_1})\sin{\Delta_{31}}) - \mu_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2\sin{\Delta_{21}} + \mu_3l_1 l_3\dot{\theta}_1\dot{\theta}_3\sin{\Delta_{31}} + \mu_1 gl_1\cos{\theta_1}\\
&= M_1 l_1^2 \ddot{\theta}_1 + \mu_2 l_1l_2 (\ddot{\theta}_2\cos{\Delta_{21}} - [\dot{\theta}_2(\dot{\theta}_2-\theta_1)+\dot{\theta}_1\dot{\theta}_2]\sin{\Delta_{21}}) + \mu_3 l_1 l_3(\ddot{\theta}_3\cos{\Delta_{31}} - [\dot{\theta}_3(\dot{\theta}_3-\dot{\theta_1})+\dot{\theta}_1\dot{\theta}_3]\sin{\Delta_{31}}) + \mu_1 gl_1\cos{\theta_1} \\
&= M_1 l_1^2 \ddot{\theta}_1 + \mu_2 l_1l_2 (\ddot{\theta}_2\cos{\Delta_{21}} - \dot{\theta}_2(\dot{\theta}_2-\theta_1)\sin{\Delta_{21}}) + \mu_3 l_1 l_3(\ddot{\theta}_3\cos{\Delta_{31}} - \dot{\theta}_3(\dot{\theta}_3-\dot{\theta_1})\sin{\Delta_{31}}) - \mu_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2\sin{\Delta_{21}} + \mu_3l_1 l_3\dot{\theta}_1\dot{\theta}_3\sin{\Delta_{31}} + \mu_1 gl_1\cos{\theta_1}\\
&= M_1 l_1^2 \ddot{\theta}_1 + \mu_2 l_1l_2 (\ddot{\theta}_2\cos{\Delta_{21}} - \dot{\theta}_2^2\sin{\Delta_{21}}) + \mu_3 l_1 l_3(\ddot{\theta}_3\cos{\Delta_{31}} - \dot{\theta}_3^2\sin{\Delta_{31}}) + \mu_1 gl_1\cos{\theta_1}.
\end{align*}

### $\theta_2$
\begin{align*}
p_{\theta_2} &= \dfrac{\partial \mathcal{L}}{\partial \dot{\theta}_2} \\
&= M_2 l_2^2 \dot{\theta}_2 + \mu_2 l_1l_2 \dot{\theta}_1\cos{\Delta_{21}} + \mu_3 l_2l_3 \dot{\theta}_3 \cos{\Delta_{32}}\\
\dot{p}_{\theta_2} &= M_2 l_2^2 \ddot{\theta}_2 + \mu_2 l_1l_2 (\ddot{\theta}_1\cos{\Delta_{21}} - \dot{\theta}_1 (\dot{\theta}_2-\dot{\theta}_1)\sin{\Delta_{21}})+ \mu_3 l_2l_3 (\ddot{\theta}_3 \cos{\Delta_{32}} -\dot{\theta}_3 (\dot{\theta}_3-\dot{\theta}_2)\sin{\Delta_{32}})\\
F_{\theta_2} &= -\mu_2 l_2(l_1\dot{\theta}_1\dot{\theta}_2 \dfrac{\partial \Delta_{21}}{\partial \theta_2}\sin{\Delta_{21}}+g\cos{\theta_2}) - \mu_3 l_2l_3\dot{\theta}_2\dot{\theta}_3\dfrac{\partial \Delta_{32}}{\partial \theta_2}\sin{\Delta_{32}}\\
&= -\mu_2 l_2(l_1\dot{\theta}_1\dot{\theta}_2 \sin{\Delta_{21}}+g\cos{\theta_2}) + \mu_3 l_2l_3\dot{\theta}_2\dot{\theta}_3\sin{\Delta_{32}}\\
\delta'_{\theta_2} \mathcal{L} &= M_2 l_2^2 \ddot{\theta}_2 + \mu_2 l_1l_2 (\ddot{\theta}_1\cos{\Delta_{21}} - \dot{\theta}_1 (\dot{\theta}_2-\dot{\theta}_1)\sin{\Delta_{21}})+ \mu_3 l_2l_3 (\ddot{\theta}_3 \cos{\Delta_{32}} -\dot{\theta}_3 (\dot{\theta}_3-\dot{\theta}_2)\sin{\Delta_{32}}) + \mu_2 l_2(l_1\dot{\theta}_1\dot{\theta}_2 \sin{\Delta_{21}}+g\cos{\theta_2}) - \mu_3 l_2l_3\dot{\theta}_2\dot{\theta}_3\sin{\Delta_{32}} \\
&= M_2 l_2^2 \ddot{\theta}_2 + \mu_2 l_1l_2 (\ddot{\theta}_1\cos{\Delta_{21}} + \dot{\theta}_1^2\sin{\Delta_{21}})+ \mu_3 l_2l_3 (\ddot{\theta}_3 \cos{\Delta_{32}} -\dot{\theta}_3^2\sin{\Delta_{32}}) +\mu_2 l_2g\cos{\theta_2} \\
&= M_2 l_2^2 \ddot{\theta}_2 + \mu_2 l_2 (l_1(\ddot{\theta}_1\cos{\Delta_{21}} + \dot{\theta}_1^2\sin{\Delta_{21}})+g\cos{\theta_2})+ \mu_3 l_2l_3 (\ddot{\theta}_3 \cos{\Delta_{32}} -\dot{\theta}_3^2\sin{\Delta_{32}}).
\end{align*}

### $\theta_3$
\begin{align*}
p_{\theta_3} &= \dfrac{\partial \mathcal{L}}{\partial \dot{\theta}_3} \\
&= M_3 l_3^2 \dot{\theta}_3 + \mu_3 l_3(l_2 \dot{\theta}_2\cos{\Delta_{32}}+l_1\dot{\theta}_1\cos{\Delta_{31}}) \\
\dot{p}_{\theta_3} &= M_3 l_3^2 \ddot{\theta}_3 + \mu_3 l_3(l_2 (\ddot{\theta}_2\cos{\Delta_{32}} - \dot{\theta}_2 (\dot{\theta}_3-\dot{\theta}_2)\sin{\Delta_{32}})+l_1(\ddot{\theta}_1\cos{\Delta_{31}}-\dot{\theta}_1(\dot{\theta}_3-\dot{\theta}_1)\sin{\Delta_{31}})) \\
F_{\theta_3} &= \dfrac{\partial \mathcal{L}}{\partial \theta_3} \\
&= -\mu_3 l_3\left[\dot{\theta}_3 (l_2\dot{\theta}_2 \dfrac{\partial \Delta_{32}}{\partial \theta_3}\sin{\Delta_{32}}+l_1\dot{\theta}_1\dfrac{\partial \Delta_{31}}{\partial \theta_3}\sin{\Delta_{31}}) + g\cos{\theta_3}\right] \\
&= -\mu_3 l_3\left[\dot{\theta}_3 (l_2\dot{\theta}_2 \sin{\Delta_{32}}+l_1\dot{\theta}_1\sin{\Delta_{31}}) + g\cos{\theta_3}\right] \\
\delta'_{\theta_3} \mathcal{L} &= M_3l_3^2 \ddot{\theta}_3 + \mu_3 l_3(l_2 (\ddot{\theta}_2\cos{\Delta_{32}} - \dot{\theta}_2 (\dot{\theta}_3-\dot{\theta}_2)\sin{\Delta_{32}})+l_1(\ddot{\theta}_1\cos{\Delta_{31}}-\dot{\theta}_1(\dot{\theta}_3-\dot{\theta}_1)\sin{\Delta_{31}})) + \mu_3 l_3\left[\dot{\theta}_3 (l_2\dot{\theta}_2 \sin{\Delta_{32}}+l_1\dot{\theta}_1\sin{\Delta_{31}}) + g\cos{\theta_3}\right]\\
&= M_3l_3^2 \ddot{\theta}_3 + \mu_3 l_3\left[l_2 (\ddot{\theta}_2\cos{\Delta_{32}} + \dot{\theta}_2^2\sin{\Delta_{32}})+l_1(\ddot{\theta}_1\cos{\Delta_{31}}+\dot{\theta}_1^2\sin{\Delta_{31}})+g\cos{\theta_3}\right].
\end{align*}

## Final system
Hence given our equations of motion are $\delta'_{\theta_i}\mathcal{L} = Q_{\theta_i}$, we could write them in matrix form as (given how long $Q_{\theta_i}$ is, we will not expand on it)

\begin{align*}
\begin{bmatrix}
M_1 l_1^2 & \mu_2 l_1l_2 \cos{\Delta_{21}} & \mu_3 l_1l_3 \cos{\Delta_{31}} \\
\mu_2 l_1l_2 \cos{\Delta_{21}} & M_2 l_2^2 & \mu_3 l_2l_3 \cos{\Delta_{32}} \\
\mu_3 l_1l_3 \cos{\Delta_{31}} & \mu_3 l_2 l_3 \cos{\Delta_{32}} & M_3 l_3^2 
\end{bmatrix} \begin{bmatrix}
\ddot{\theta}_1 \\
\ddot{\theta}_2 \\
\ddot{\theta}_3
\end{bmatrix} &= \begin{bmatrix}
Q_{\theta_1} + \mu_2 l_1l_2\dot{\theta}_2^2 \sin{\Delta_{21}} + \mu_3l_1l_3 \dot{\theta}_3^2 \sin{\Delta_{31}} - \mu_1gl_1\cos{\theta_1}\\
Q_{\theta_2} - \mu_2 l_2 (l_1\dot{\theta}_1^2\sin{\Delta_{21}}+g\cos{\theta_2}) + \mu_3 l_2l_3\dot{\theta}_3^2 \sin{\Delta_{32}}\\
Q_{\theta_3} - \mu_3 l_2l_3 \dot{\theta}_2^2 \cos{\Delta_{32}} - \mu_3 l_1 l_3 \dot{\theta}_1^2 \sin{\Delta_{31}} - \mu_3 gl_3 \cos{\theta}_3
\end{bmatrix}.
\end{align*}