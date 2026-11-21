### Double pendulum
\figure{doublependiag}{right}{600px}{/fidello/Double%20pendulum.svg}{Double pendulum diagram.}
I thought it would be instructive for me to show you how to derive the equations of motion for the simplest sort of double pendulum. One where there is no air resistance, no friction, and the rods are both massless. Then the coordinates of the first bob are:
\begin{align*}
x_1 &= l_1 \cos{\theta_1} \\
y_1 &= l_1 \sin{\theta_1}.
\end{align*}
As for the simple pendulum, $v_1^2 = l_1^2 \dot{\theta}_1^2$. As for the second bob, they would be:
\begin{align*}
x_2 &= x_1 + l_2 \cos{\theta_2} \\
&= l_1 \cos{\theta_1} + l_2\cos{\theta_2} \\
y_2 &= y_1 + l_2 \sin{\theta_2} \\
&= l_1 \sin{\theta_1} + l_2\sin{\theta_2}.
\end{align*}
Velocity components are:
\begin{align*}
\dot{x}_2 &= -l_1 \dot{\theta}_1 \sin{\theta_1} - l_2 \dot{\theta}_2 \sin{\theta_2} \\
\dot{y}_2 &= l_1 \dot{\theta}_1 \cos{\theta_1} + l_2 \dot{\theta}_2 \cos{\theta_2} \\
v_2^2 &= \dot{x}_2^2+\dot{y}_2^2 \\
&= l_1^2 \dot{\theta}_1^2 \sin^2{\theta_1} + l_2^2 \dot{\theta}_2^2\sin^2{\theta_2} + 2l_1 l_2\dot{\theta}_1\dot{\theta}_2\sin{\theta_1}\sin{\theta_2} + l_1^2 \dot{\theta}_1^2\cos^2{\theta_1} + l_2^2\dot{\theta}_2^2\cos^2{\theta_2} + 2l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \cos{\theta_1}\cos{\theta_2} \\
&= l_1^2 \dot{\theta}_1^2 + l_2^2\dot{\theta}_2^2 + 2l_1l_2 \dot{\theta}_1\dot{\theta}_2 \cos{(\theta_2-\theta_1)}.
\end{align*}
Hence the kinetic energy is:
\begin{align*}
T &= \dfrac{m_1l_1^2 \dot{\theta}_1^2}{2} + \dfrac{m_2 \left(l_1^2 \dot{\theta}_1^2 + l_2^2\dot{\theta}_2^2 + 2l_1l_2 \dot{\theta}_1\dot{\theta}_2 \cos{(\theta_2-\theta_1)}\right)}{2}.
\end{align*}
As for the potential energy:
\begin{align*}
V &= m_1 gy_1 + m_2 gy_2 \\
&= m_1 gl_1 \sin{\theta_1} + m_2 g(l_1\sin{\theta_1} + l_2\sin{\theta_2}).
\end{align*}
Hence the Lagrangian is:
\begin{align*}
L &= T -V \\
&= \dfrac{m_1l_1^2 \dot{\theta}_1^2}{2} + \dfrac{m_2 \left(l_1^2 \dot{\theta}_1^2 + l_2^2\dot{\theta}_2^2 + 2l_1l_2 \dot{\theta}_1\dot{\theta}_2 \cos{(\theta_2-\theta_1)}\right)}{2} - m_1 gl_1 \sin{\theta_1} - m_2 g(l_1\sin{\theta_1} + l_2\sin{\theta_2}) \\
&= \dfrac{l_1^2 (m_1 + m_2)}{2} \dot{\theta}_1^2 + \dfrac{m_2 l_2^2 \dot{\theta}_2^2}{2} + m_2 l_1l_2\dot{\theta}_1\dot{\theta}_2 \cos{(\theta_2-\theta_1)} - (m_1+m_2)gl_1\sin{\theta_1} - m_2 gl_2 \sin{\theta_2}.
\end{align*}
Hence the canonical momentums are:
\begin{align*}
p_{\theta_1} &= \dfrac{\partial L}{\partial \dot{\theta}_1} \\
&= (m_1+m_2) l_1^2 \dot{\theta}_1 + m_2 l_1 l_2\dot{\theta}_2 \cos{(\theta_2-\theta_1)} \\
p_{\theta_2} &= \dfrac{\partial L}{\partial \dot{\theta}_2} \\
&= m_2 l_2^2 \dot{\theta}_2 + m_2 l_1 l_2\dot{\theta}_1 \cos{(\theta_2-\theta_1)} \\
\dot{p}_{\theta_1} &= (m_1+m_2) l_1^2 \ddot{\theta}_1 + m_2 l_1 l_2\ddot{\theta}_2 \cos{(\theta_2-\theta_1)} - m_2 l_1l_2 \dot{\theta}_2(\dot{\theta}_2-\dot{\theta}_1)\sin{(\theta_2-\theta_1)}\\
\dot{p}_{\theta_2} &= m_2 l_2^2 \ddot{\theta}_2 + m_2 l_1 l_2\ddot{\theta}_1 \cos{(\theta_2-\theta_1)} - m_2 l_1l_2 \dot{\theta}_1(\dot{\theta}_2-\dot{\theta}_1)\sin{(\theta_2-\theta_1)}.
\end{align*}
As for $\dfrac{\partial L}{\partial \theta_i}$:
\begin{align*}
\dfrac{\partial L}{\partial \theta_1} &= m_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{(\theta_2-\theta_1)} - (m_1+m_2)gl_1\cos{\theta_1}\\
\dfrac{\partial L}{\partial \theta_2} &= -m_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{(\theta_2-\theta_1)} - m_2 gl_2\cos{\theta_2}.
\end{align*}
Hence the Euler-Lagrange equations are:
\begin{align*}
(m_1+m_2) l_1^2 \ddot{\theta}_1 + m_2 l_1 l_2\ddot{\theta}_2 \cos{(\theta_2-\theta_1)} - m_2 l_1l_2 \dot{\theta}_2(\dot{\theta}_2-\dot{\theta}_1)\sin{(\theta_2-\theta_1)} &= m_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{(\theta_2-\theta_1)} - (m_1+m_2)gl_1\cos{\theta_1} \\
m_2 l_2^2 \ddot{\theta}_2 + m_2 l_1 l_2\ddot{\theta}_1 \cos{(\theta_2-\theta_1)} - m_2 l_1l_2 \dot{\theta}_1(\dot{\theta}_2-\dot{\theta}_1)\sin{(\theta_2-\theta_1)} &= -m_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{(\theta_2-\theta_1)} - m_2 gl_2\cos{\theta_2}.
\end{align*}
Let us get the second derivatives onto the left with everything else on the right-hand side of the equations:
\begin{align*}
(m_1+m_2) l_1^2 \ddot{\theta}_1 + m_2 l_1 l_2\ddot{\theta}_2 \cos{(\theta_2-\theta_1)} &= m_2 l_1l_2 \dot{\theta}_2(\dot{\theta}_2-\dot{\theta}_1)\sin{(\theta_2-\theta_1)} + m_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{(\theta_2-\theta_1)} - (m_1+m_2)gl_1\cos{\theta_1} \\
&= m_2 l_1l_2 \dot{\theta}_2^2\sin{(\theta_2-\theta_1)} - (m_1+m_2)gl_1\cos{\theta_1} \\
m_2 l_2^2 \ddot{\theta}_2 + m_2 l_1 l_2\ddot{\theta}_1 \cos{(\theta_2-\theta_1)} &= m_2 l_1l_2 \dot{\theta}_1(\dot{\theta}_2-\dot{\theta}_1)\sin{(\theta_2-\theta_1)} -m_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2 \sin{(\theta_2-\theta_1)} - m_2 gl_2\cos{\theta_2} \\
&= -m_2 l_1l_2 \dot{\theta}_1^2\sin{(\theta_2-\theta_1)} - m_2 gl_2\cos{\theta_2}.
\end{align*}
Ideally, we want one equation for $\ddot{\theta}_1$ (and $\ddot{\theta}_1$ alone) and another $\ddot{\theta}_2$. So let us put our system into matrix form:
\begin{align*}
\begin{bmatrix}
(m_1+m_2)l_1^2 & m_2 l_1 l_2 \cos{(\theta_2-\theta_1)} \\
m_2 l_1l_2 \cos{(\theta_2-\theta_1)} & m_2 l_2^2
\end{bmatrix} \begin{bmatrix}
\ddot{\theta_1} \\
\ddot{\theta_2}
\end{bmatrix} &= \begin{bmatrix}
m_2 l_1l_2 \dot{\theta}_2^2\sin{(\theta_2-\theta_1)} - (m_1+m_2)gl_1\cos{\theta_1} \\
-m_2 l_1l_2 \dot{\theta}_1^2\sin{(\theta_2-\theta_1)} - m_2 gl_2\cos{\theta_2}
\end{bmatrix}.
\end{align*}
Inverting the square matrix on the left, which we will call $A$, we get:
\begin{align*}
A^{-1} &= \dfrac{1}{(m_1+m_2)l_1^2 \cdot m_2 l_2^2 - (m_2 l_1 l_2 \cos{(\theta_2-\theta_1)})^2} \begin{bmatrix}
m_2 l_2^2 & -m_2 l_1 l_2 \cos{(\theta_2-\theta_1)}\\
-m_2 l_1l_2 \cos{(\theta_2-\theta_1)} & (m_1+m_2)l_1^2 
\end{bmatrix} \\
&= \dfrac{1}{m_1m_2l_1^2 l_2^2 + m_2^2 l_1^2 l_2^2 - m_2^2 l_1^2l_2^2 \cos^2{(\theta_2-\theta_1)}} \begin{bmatrix}
m_2 l_2^2 & -m_2 l_1 l_2 \cos{(\theta_2-\theta_1)}\\
-m_2 l_1l_2 \cos{(\theta_2-\theta_1)} & (m_1+m_2)l_1^2 
\end{bmatrix} \\
&= \dfrac{1}{m_1m_2l_1^2 l_2^2 + m_2^2 l_1^2 l_2^2\sin^2{(\theta_2-\theta_1)}} \begin{bmatrix}
m_2 l_2^2 & -m_2 l_1 l_2 \cos{(\theta_2-\theta_1)}\\
-m_2 l_1l_2 \cos{(\theta_2-\theta_1)} & (m_1+m_2)l_1^2 
\end{bmatrix} \\
&= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
m_2 l_2^2 & -m_2 l_1 l_2 \cos{(\theta_2-\theta_1)}\\
-m_2 l_1l_2 \cos{(\theta_2-\theta_1)} & (m_1+m_2)l_1^2 
\end{bmatrix}.
\end{align*}
Hence the solution to this system of equations is:
\begin{align*}
\begin{bmatrix}
\ddot{\theta_1} \\
\ddot{\theta_2}
\end{bmatrix} &= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
m_2 l_2^2 & -m_2 l_1 l_2 \cos{(\theta_2-\theta_1)}\\
-m_2 l_1l_2 \cos{(\theta_2-\theta_1)} & (m_1+m_2)l_1^2 
\end{bmatrix} \begin{bmatrix}
m_2 l_1l_2 \dot{\theta}_2^2\sin{(\theta_2-\theta_1)} - (m_1+m_2)gl_1\cos{\theta_1} \\
-m_2 l_1l_2 \dot{\theta}_1^2\sin{(\theta_2-\theta_1)} - m_2 gl_2\cos{\theta_2}
\end{bmatrix} \\
&= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
m_2^2 l_1 l_2^3 \dot{\theta}_2^2 \sin{(\theta_2-\theta_1)} - (m_1+m_2)m_2 gl_1 l_2^2 \cos{\theta_1} + m_2^2 l_1^2 l_2^2 \dot{\theta}_1^2 \cos{(\theta_2-\theta_1)}\sin{(\theta_2-\theta_1)} + m_2^2 gl_1 l_2^2 \cos{(\theta_2-\theta_1)}\cos{\theta_2} \\
-m_2^2 l_1^2 l_2^2 \dot{\theta}_2^2 \cos{(\theta_2-\theta_1)}\sin{(\theta_2-\theta_1)} + (m_1+m_2)gl_1^2 l_2\cos{(\theta_2-\theta_1)}\cos{\theta_1} - (m_1+m_2)m_2l_1^3 l_2 \dot{\theta}_1^2 \sin{(\theta_2-\theta_1)} - (m_1+m_2)m_2 gl_1^2 l_2 \cos{\theta_2}
\end{bmatrix} \\
&= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
\sin{(\theta_2-\theta_1)} \left[m_2^2 l_1 l_2^3 \dot{\theta}_2^2 + m_2^2 l_1^2 l_2^2 \dot{\theta}_1^2 \cos{(\theta_2-\theta_1)}\right]- m_1 m_2 gl_1 l_2^2 \cos{\theta_1} + m_2^2 gl_1 l_2^2 (\cos{(\theta_2-\theta_1)}\cos{\theta_2}-\cos{\theta_1})  \\
-m_2 l_1^2 l_2\sin{(\theta_2-\theta_1)}\left[m_2 l_2 \dot{\theta}_2^2 \cos{(\theta_2-\theta_1)}+(m_1+m_2)l_1 \dot{\theta}_1^2\right] + (m_1+m_2)gl_1^2 l_2\left[\cos{(\theta_2-\theta_1)}\cos{\theta_1}  - \cos{\theta_2}\right]
\end{bmatrix} \\
&= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
\sin{(\theta_2-\theta_1)} \left[m_2^2 l_1 l_2^3 \dot{\theta}_2^2 + m_2^2 l_1^2 l_2^2 \dot{\theta}_1^2 \cos{(\theta_2-\theta_1)}\right]- m_1 m_2 gl_1 l_2^2 \cos{\theta_1} + m_2^2 gl_1 l_2^2 (\cos^2{\theta_2}\cos{\theta_1}+\sin{\theta_2}\sin{\theta_1}\cos{\theta_2}-\cos{\theta_1})  \\
-m_2 l_1^2 l_2\sin{(\theta_2-\theta_1)}\left[m_2 l_2 \dot{\theta}_2^2 \cos{(\theta_2-\theta_1)}+(m_1+m_2)l_1 \dot{\theta}_1^2\right] + (m_1+m_2)gl_1^2 l_2\left[\cos{\theta_2}\cos^2{\theta_1}+\sin{\theta_1}\sin{\theta_2}\cos{\theta_1}  - \cos{\theta_2}\right]
\end{bmatrix} \\
&= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
\sin{(\theta_2-\theta_1)} \left[m_2^2 l_1 l_2^3 \dot{\theta}_2^2 + m_2^2 l_1^2 l_2^2 \dot{\theta}_1^2 \cos{(\theta_2-\theta_1)}\right]- m_1 m_2 gl_1 l_2^2 \cos{\theta_1} + m_2^2 gl_1 l_2^2 \left[-\cos{\theta_1}\sin^2{\theta_2}  + \sin{\theta_1}\sin{\theta_2}\cos{\theta_2}\right]  \\
-m_2 l_1^2 l_2\sin{(\theta_2-\theta_1)}\left[m_2 l_2 \dot{\theta}_2^2 \cos{(\theta_2-\theta_1)}+(m_1+m_2)l_1 \dot{\theta}_1^2\right] + (m_1+m_2)gl_1^2 l_2\left[-\cos{\theta_2}\sin^2{\theta_1}  + \sin{\theta_1}\sin{\theta_2}\cos{\theta_1}\right]
\end{bmatrix} \\
&= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
\sin{(\theta_2-\theta_1)} \left[m_2^2 l_1 l_2^3 \dot{\theta}_2^2 + m_2^2 l_1^2 l_2^2 \dot{\theta}_1^2 \cos{(\theta_2-\theta_1)}\right]- m_1 m_2 gl_1 l_2^2 \cos{\theta_1} + m_2^2 gl_1 l_2^2 \sin{\theta_2}\left[-\cos{\theta_1}\sin{\theta_2}  + \sin{\theta_1}\cos{\theta_2}\right]  \\
-m_2 l_1^2 l_2\sin{(\theta_2-\theta_1)}\left[m_2 l_2 \dot{\theta}_2^2 \cos{(\theta_2-\theta_1)}+(m_1+m_2)l_1 \dot{\theta}_1^2\right] + (m_1+m_2)gl_1^2 l_2 \sin{\theta_1}\left[-\cos{\theta_2}\sin{\theta_1}  + \sin{\theta_2}\cos{\theta_1}\right]
\end{bmatrix} \\
&= \dfrac{1}{(m_1+m_2\sin^2{(\theta_2-\theta_1)})m_2l_1^2 l_2^2} \begin{bmatrix}
\sin{(\theta_2-\theta_1)} \left[m_2^2 l_1 l_2^3 \dot{\theta}_2^2 + m_2^2 l_1^2 l_2^2 \dot{\theta}_1^2 \cos{(\theta_2-\theta_1)}\right]- m_1 m_2 gl_1 l_2^2 \cos{\theta_1} - m_2^2 gl_1 l_2^2 \sin{\theta_2}\sin{(\theta_2-\theta_1)}  \\
-m_2 l_1^2 l_2\sin{(\theta_2-\theta_1)}\left[m_2 l_2 \dot{\theta}_2^2 \cos{(\theta_2-\theta_1)}+(m_1+m_2)l_1 \dot{\theta}_1^2\right] + (m_1+m_2)gl_1^2 l_2 \sin{\theta_1}\sin{(\theta_2-\theta_1)}
\end{bmatrix}.
\end{align*}