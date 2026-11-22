### PDEs: Separation of variables
I thought I should also cover one technique for solving partial differential equations, too. Say we had the heat equation, Equation \eqref{heat}, we would assume that $u(x,t)$ can be written as $u(x,t)=X(x)T(t)$. Then Equation \eqref{heat} becomes:
\begin{align*}
X(x)\dot{T}(t) &= \alpha X''(x)T(t).
\end{align*}
Dividing both sides by $u(x,t)$ yields:
\begin{align*}
\dfrac{\dot{T}(t)}{T(t)} &= \alpha \dfrac{X''(x)}{X(x)}.
\end{align*}
The left-hand side depends only on $t$ and the right-hand side only depends on $x$. Hence both sides must equal a constant, which we will call $-\lambda$. The importance of the negative sign will become apparent soon. Hence we can write a pair of ordinary differential equations we need to solve:
\begin{align*}
\dfrac{\dot{T}(t)}{T(t)} &= -\lambda \\
\dot{T}(t) &= -\lambda T(t) \\
\alpha \dfrac{X''(x)}{X(x)} &= -\lambda \\
\alpha X''(x) &= -\lambda X(x) \\
X''(x) &= \dfrac{-\lambda}{\alpha} X(x).
\end{align*}
By assuming exponential solutions, we get:
\begin{align*}
T(t) &= C_1 e^{-\lambda t} \\
X(x) &= C_2 e^{i\sqrt{\dfrac{\lambda}{\alpha}}x} + C_3 e^{-i\sqrt{\dfrac{\lambda}{\alpha}}x}.
\end{align*}
$C_1$, $C_2$ and $C_3$ are constants that depend on our boundary and initial conditions. Since $e^{i\theta}\cos{\theta}+i\sin{\theta}$:
\begin{align*}
X(x) &= C_2 \cos{\sqrt{\dfrac{\lambda}{\alpha}}x} + C_3 \sin{\sqrt{\dfrac{\lambda}{\alpha}}x}.
\end{align*}
The $C_2$ and $C_3$ we used here will be different to those used when we presented $X(x)$ in exponential form. Then we apply our boundary and initial conditions, which are often of the form:
\begin{align*}
u(x,0)&=f(x) \\
u(0,t)&=0 \\
u(L,t)&=0.
\end{align*}
Where $x\in [0,L]$ is our spatial domain and $f(x)$ is some specified function. In this event, we impose our spatial boundary conditions on $X(x)$:
\begin{align*}
X(0) &= 0 \\
\implies C_2 \cos{0} + C_3\sin{0} &= 0 \\
C_2 &= 0 \\
X(L) &= 0 \\
C_3\sin{\sqrt{\dfrac{\lambda}{\alpha}}L} &= 0.
\end{align*}
This, in turn, implies that
\begin{align*}
\sqrt{\dfrac{\lambda}{\alpha}}L &= n\pi,\hspace{0.2cm}n\in \mathbb{Z}.
\end{align*}
Squaring both sides yields:
\begin{align*}
\dfrac{\lambda}{\alpha} L^2 &= n^2 \pi^2 \\
\lambda &= \dfrac{\alpha n^2 \pi^2}{L^2}.
\end{align*}
Hence we do not have a just a singular solution for $\lambda$ and we should really call $\lambda_n$ to reflect this reality. $X(x)$, $T(t)$ and $u(x,t)$ all depend on $\lambda_n$ and hence should get a $n$ subscript.
\begin{align*}
\lambda_n &= \dfrac{\alpha n^2 \pi^2}{L^2} \\
X_n(x) &= C_{3n} \sin{\sqrt{\dfrac{n^2 \pi^2}{L^2} x}} \\
&= C_{3n} \sin{\dfrac{|n| \pi}{L} x} \\
T_n(t) &= C_{1n}e^{-\dfrac{\alpha n^2 \pi^2}{L^2} t}.
\end{align*}
Why have I given the constants $n$ indices? Because they, too, could vary between terms with different $\lambda$. Why the $|n|$ well we square rooted $n^2$, so it will be the absolute value of $n$. We will call $C_{1n}C_{3n}$, $a_n$ when we write $u_n(x,t)$:
Hence:
\begin{align*}
u_n(x,t) &= a_n \sin{\dfrac{|n| \pi}{L} x} e^{-\dfrac{\alpha n^2 \pi^2}{L^2} t}.
\end{align*}
Keep in mind that $u_{-n}(x,t)=u_n(x,t)$ and that $u_0=0$ (as $\sin{0}=0$) we obtain the following series as our solution:
\begin{align*}
u(x,t) &= \sum_{n=1}^{\infty} a_n \sin{\dfrac{n \pi}{\sqrt{\alpha} L} x} e^{-\dfrac{n^2 \pi^2}{L^2} t}.
\end{align*}
Imposing that $u(x,0)=f(x)$ yields:
\begin{align*}
f(x) &= \sum_{n=0}^{\infty} a_n \sin{\dfrac{n \pi}{\sqrt{\alpha} L} x}.
\end{align*}
How can we determine $a_n$ from this? Well, it turns out that if we multiply this by $\sin{\dfrac{m \pi}{L} x}$ and integrate over $x\in [0,L]$ we get:
\begin{align*}
\int_0^L f(x) \sin{\dfrac{m \pi}{\sqrt{\alpha} L} x} dx &= \sum_{n=0}^{\infty} a_n \int_0^L\sin{\dfrac{n \pi}{L} x}\sin{\dfrac{m \pi}{L} x} dx.
\end{align*}
Interestingly:
\begin{align*}
\sin{A}\sin{B} &= \dfrac{1}{2}(\cos{(A-B)-\cos{A+B}}).
\end{align*}
This yields:
\begin{align*}
\sin{\dfrac{n \pi}{L} x}\sin{\dfrac{m \pi}{L} x} &= \dfrac{1}{2}\left(\cos{\dfrac{m-n}{\pi}{L}x} - \cos{\dfrac{m+n}{\pi}{L}x}\right).
\end{align*}
Hence our integral becomes:
\begin{align*}
\int_0^L f(x) \sin{\dfrac{m \pi}{\sqrt{\alpha} L} x} dx &= \dfrac{1}{2}\sum_{n=0}^{\infty} a_n \int_0^L \left(\cos{\dfrac{(m-n)\pi}{L}x} - \cos{\dfrac{(m+n)\pi}{L}x}\right) dx \\
&= \dfrac{1}{2}\sum_{n=0}^{\infty} a_n \left(\left[\dfrac{L}{(m-n)\pi}\sin{\dfrac{(m-n)\pi}{L}x}\right]_0^L - \left[\dfrac{L}{(m+n)\pi}\sin{\dfrac{m+n}{\pi}{L}x}\right]_0^L\right).
\end{align*}
Typically we will get zero from these integrals. Evaluating at $x=0$ will always yield $0$ for both integrals. $x=L$ will also yield $0$ if $m\neq n$. If $m=n$ though, we end up dividing by 0 for the first integral, which means we cannot simply declare it zero. Although the second antiderivative, when evaluated at $x=L$ yields $0$. So we will need to take the limit as $m$ approaches $n$ to determine the value of the first integral at its upper bound:
\begin{align*}
\lim_{m\rightarrow n} \dfrac{L}{(m-n)\pi}\sin{\dfrac{(m-n)\pi}{L}L} &= \lim_{m\rightarrow n} \dfrac{L}{(m-n)\pi}\sin{(m-n)\pi}.
\end{align*}
Let $x=(m-n)\pi$:
\begin{align*}
\lim_{m\rightarrow n} \dfrac{L}{(m-n)\pi}\sin{(m-n)\pi} &= \lim_{x\rightarrow 0} \dfrac{L}{x} \sin{x}.
\end{align*}
We can use the Maclaurin series to approximate $\sin{x}$ as $x\rightarrow 0$:
\begin{align*}
\lim_{x\rightarrow 0} \dfrac{L}{x} \sin{x} &= \lim_{x\rightarrow 0} \dfrac{L}{x} x \\
&= \lim_{x\rightarrow 0} L \\
&= L.
\end{align*}
Hence:
\begin{align*}
\left[\dfrac{L}{(m-n)\pi}\sin{\dfrac{(m-n)\pi}{L}x}\right]_0^L &= L\delta_{mn} \\
\left[\dfrac{L}{(m+n)\pi}\sin{\dfrac{m+n}{\pi}{L}x}\right]_0^L &= 0.
\end{align*}
Hence:
\begin{align*}
\int_0^L f(x) \sin{\dfrac{m \pi}{\sqrt{\alpha} L} x} dx &= \sum_{n=0}^{\infty} a_n L \delta_{mn} \\
&= La_m.
\end{align*}
Here we have used the fact that $\delta_{mn}$ is zero if $m\neq n$ and one if $m=n$ and hence the only non-zero term corresponds to $m=n$. Hence we can calculate $a_n$ as:
\begin{align*}
a_n &= \dfrac{1}{L} \int_0^L f(x) \sin{\dfrac{n \pi}{\sqrt{\alpha} L} x} dx.
\end{align*}