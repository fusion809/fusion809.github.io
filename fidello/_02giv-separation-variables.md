### PDEs: Separation of variables
I thought I should also cover one technique for solving partial differential equations, too. Say we had the heat equation, Equation \eqref{heat}, we would assume that $u(x,t)$ can be written as $u(x,t)=X(x)T(t)$. Then Equation \eqref{heat} becomes:
\begin{align*}
X(x)\dot{T}(t) &= \alpha X''(x)T(t).
\end{align*}
Dividing both sides by $u(x,t)$ yields:
\begin{align*}
\dfrac{\dot{T}(t)}{T(t)} &= \alpha \dfrac{X''(x)}{X(x)}.
\end{align*}
The left-hand side depends only on $t$ and the right-hand side only depends on $x$. Hence both sides must equal a constant, which we will call $-\lambda$. The importance of the negative sign will become apparent soon. Hence we can write a pair of ordinary differential equations we need to solve. For $T(t)$, the equation is:
\begin{align*}
\dfrac{\dot{T}(t)}{T(t)} &= -\lambda \\
\dot{T}(t) &= -\lambda T(t)
\end{align*}
As for $X(x)$, its equation is:
\begin{align*}
\alpha \dfrac{X''(x)}{X(x)} &= -\lambda \\
\alpha X''(x) &= -\lambda X(x) \\
X''(x) &= -\dfrac{\lambda}{\alpha} X(x).
\end{align*}
These are both linear, homogeneous differential equations with constant coefficients. Hence, our solutions will involve exponential functions. Solving the polynomial equations we get from this knowledge yields:
\begin{align}
T(t) &= C_1 e^{-\lambda t} \label{Tt}\\
X(x) &= C_2 e^{i\sqrt{\dfrac{\lambda}{\alpha}}x} + C_3 e^{-i\sqrt{\dfrac{\lambda}{\alpha}}x}.\label{XxP}
\end{align}
$C_1$, $C_2$ and $C_3$ are constants that depend on our boundary and initial conditions. Since $e^{i\theta}=\cos{\theta}+i\sin{\theta}$:
\begin{align}
X(x) &= C_2 \cos{\sqrt{\dfrac{\lambda}{\alpha}}x} + C_3 \sin{\sqrt{\dfrac{\lambda}{\alpha}}x}.\label{Xx}
\end{align}
So this is why we had opted to give $\lambda$ a negative sign earlier. As real exponentials will not satisfy our boundary conditions of $X(0)=X(L)=0$ without $C_2$ and $C_3$ both being zero. The $C_2$ and $C_3$ we used in Equation \eqref{Xx} will be different to those in Equation \eqref{XxP}. Then we apply our boundary and initial conditions, which are often of the form:
\begin{align*}
u(x,0)&=f(x) \\
u(0,t)&=0 \\
u(L,t)&=0.
\end{align*}
Where $x\in [0,L]$ is our spatial domain and $f(x)$ is some specified function. We will impose our spatial boundary conditions on $X(x)$, $X(0)=0$ implies that:
\begin{align*}
C_2 \cos{0} + C_3\sin{0} &= 0.
\end{align*}
$\sin{0}=0$, $\cos{0}=1$ and hence:
\begin{align*}
C_2 &= 0.
\end{align*}
As for our other boundary condition:
\begin{align*}
X(L) &= 0 \\
C_3\sin{\sqrt{\dfrac{\lambda}{\alpha}}L} &= 0.
\end{align*}
This, in turn, implies that:
\begin{align*}
\sqrt{\dfrac{\lambda}{\alpha}}L &= n\pi,\hspace{0.2cm}n\in \mathbb{Z}.
\end{align*}
Squaring both sides yields:
\begin{align*}
\dfrac{\lambda}{\alpha} L^2 &= n^2 \pi^2
\end{align*}
Multiplying by $\dfrac{\alpha}{L^2}$ yields:
\begin{align*}
\lambda &= \dfrac{\alpha n^2 \pi^2}{L^2}.
\end{align*}
Hence we do not have a just a singular solution for $\lambda$ and we should really call it $\lambda_n$ to reflect this reality. $X(x)$, $T(t)$ and $u(x,t)$ all depend on $\lambda_n$ and hence should get a $n$ subscript.
\begin{align*}
\lambda_n &= \dfrac{\alpha n^2 \pi^2}{L^2} \\
X_n(x) &= C_{3n} \sin{\sqrt{\dfrac{n^2 \pi^2}{L^2} x}} \\
&= C_{3n} \sin{\dfrac{|n| \pi}{L} x} \\
T_n(t) &= C_{1n}e^{-\dfrac{\alpha n^2 \pi^2}{L^2} t}.
\end{align*}
Why have I given the constants $n$ indices? Because they, too, could vary between terms with different $\lambda$. Why the $|n|$? Well, we square rooted $n^2$, so it will be the absolute value of $n$. We will call $C_{1n}C_{3n}$, $a_n$ when we write $u_n(x,t)$:
Hence:
\begin{align*}
u_n(x,t) &= a_n \sin{\dfrac{|n| \pi}{L} x} e^{-\dfrac{\alpha n^2 \pi^2}{L^2} t}.
\end{align*}
Keep in mind that $u_{-n}(x,t)=u_n(x,t)$ and that $u_0=0$ (as $\sin{0}=0$), hence we only need to sum over $n \in \mathbb{Z}^{+}$. Hence, we obtain the following series as our solution:
\begin{align*}
u(x,t) &= \sum_{n=1}^{\infty} a_n \sin{\dfrac{n \pi}{L} x} e^{-\dfrac{n^2 \pi^2}{L^2} t}.
\end{align*}
Imposing that $u(x,0)=f(x)$ yields:
\begin{align}
f(x) &= \sum_{n=0}^{\infty} a_n \sin{\dfrac{n \pi}{L} x}.\label{fxFourier}
\end{align}
How can we determine $a_n$ from this? Well, it turns out that if we multiply this by $\sin{\dfrac{m \pi}{L} x}$ and integrate over $x\in [0,L]$ we get:
\begin{align*}
\int_0^L f(x) \sin{\dfrac{m \pi}{L} x} dx &= \sum_{n=0}^{\infty} a_n \int_0^L\sin{\dfrac{n \pi}{L} x}\sin{\dfrac{m \pi}{L} x} dx.
\end{align*}
Interestingly:
\begin{align*}
\sin{A}\sin{B} &= \dfrac{1}{2}(\cos{(A-B)-\cos{(A+B)}}).
\end{align*}
This yields:
\begin{align*}
\sin{\left(\dfrac{n \pi}{L} x\right)}\sin{\left(\dfrac{m \pi}{L} x\right)} &= \dfrac{1}{2}\left(\cos{\left(\dfrac{m-n}{\pi}{L}x\right)} - \cos{\left(\dfrac{m+n}{\pi}{L}x\right)}\right).
\end{align*}
Hence our integral becomes:
\begin{align*}
\int_0^L f(x) \sin{\dfrac{m \pi}{L} x} dx &= \dfrac{1}{2}\sum_{n=0}^{\infty} a_n \int_0^L \left(\cos{\left(\dfrac{(m-n)\pi}{L}x\right)} - \cos{\left(\dfrac{(m+n)\pi}{L}x\right)}\right) dx \\
&= \dfrac{1}{2}\sum_{n=0}^{\infty} a_n \left(\left[\dfrac{L}{(m-n)\pi}\sin{\left(\dfrac{(m-n)\pi}{L}x\right)}\right]_0^L - \left[\dfrac{L}{(m+n)\pi}\sin{\left(\dfrac{(m+n)\pi}{L}x\right)}\right]_0^L\right) \\
&= \dfrac{1}{2}\sum_{n=0}^{\infty} a_n \left(\dfrac{L}{(m-n)\pi}\sin{((m-n)\pi)}-\dfrac{L}{(m-n)\pi}\sin{0} - \left[\dfrac{L}{(m+n)\pi}\sin{((m+n)\pi)}-\dfrac{L}{(m+n)\pi}\sin{0}\right]\right).
\end{align*}
$m+n$ and $m-n$ is always an integer so the sines should be zero. $\sin{0}=0$, too. Hence, we will get $0$ from this if $m\neq n$. If $m=n$ though, we end up dividing by 0 for the first term, which means we cannot simply declare it zero. So we will need to take the limit as $m$ approaches $n$ to determine the value of the first integral at its upper bound:
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
\int_0^L f(x) \sin{\dfrac{m \pi}{L} x} dx &= \dfrac{1}{2}\sum_{n=0}^{\infty} a_n L \delta_{mn} \\
&= \dfrac{La_m}{2}.
\end{align*}
Here we have used the fact that $\delta_{mn}$ is zero if $m\neq n$ and one if $m=n$ and hence the only non-zero term corresponds to $m=n$. Hence we can calculate $a_n$ as:
\begin{align*}
a_n &= \dfrac{2}{L} \int_0^L f(x) \sin{\dfrac{n \pi}{L} x} dx.
\end{align*}
Interestingly, Equation \eqref{fxFourier} is called the Fourier series decomposition of $f(x)$ in terms of sines only. Fourier series are a way of rewriting periodic functions in terms of sine and cosine waves. Specifically in general, if we wish to write $f(x)$ as a Fourier series over the domain $x \in [-L,L]$ we would write (keep in mind the series coefficients $a_n$ and $b_n$ are not the same as those we have been using so far):
\begin{align*}
f(x) &= a_0 + \sum_{n=1}^{\infty} a_n \cos{\dfrac{n\pi}{L}x} + b_n \sin{\dfrac{n\pi}{L}x}. 
\end{align*}
Where:
\begin{align*}
a_0 &= \dfrac{1}{2L} \int_{-L}^{L} f(x) dx \\
a_n &= \dfrac{1}{L} \int_{-L}^{L} f(x) \cos{\left(\dfrac{2\pi nx}{L}\right)} dx \\
b_n &= \dfrac{1}{L} \int_{-L}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx.
\end{align*}
We can see our $b_n$ here matches the $a_n$ from the $u_n(x,t)$, except that $a_n$ is calculated by integrating only over $x\in [0,L]$. So what is up with this? Well, interestingly we can only represent $f(x)$ solely with sine terms if and only if $f(x)$ is an odd function, just like $\sin{\left(\dfrac{n\pi x}{L}\right)}$ itself. Meaning that $f(-x)=-f(x)$. Hence we can be split our $b_n$ integral as:
\begin{align*}
b_n &= \dfrac{1}{L} \int_{-L}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx \\
&= \dfrac{1}{L} \left[\int_{-L}^{0} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx + \int_{0}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx \right]
\end{align*}
In our first integral, we will replace $x$ with $-x$ and hence $dx\rightarrow-dx$ and our domain becomes $x\in [L, 0]$ instead of $x\in [-L, 0]$. 
\begin{align*}
b_n &= \dfrac{1}{L} \left[\int_{L}^{0} f(-x) \sin{\left(-\dfrac{2\pi nx}{L}\right)} -dx + \int_{0}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx \right]\\
&= \dfrac{1}{L} \left[\int_{L}^{0} -f(x) \cdot -\sin{\left(\dfrac{2\pi nx}{L}\right)} \cdot -dx + \int_{0}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx \right]\\
&= \dfrac{1}{L} \left[-\int_{L}^{0} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx + \int_{0}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx \right]\\
&= \dfrac{1}{L} \left[\int_{0}^L f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx + \int_{0}^{L} f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx \right]\\
&= \dfrac{2}{L} \int_{0}^L f(x) \sin{\left(\dfrac{2\pi nx}{L}\right)} dx.
\end{align*}
So we see that our $a_n$ from $u_n(x,t)$ are in fact equivalent to $b_n$, at least for odd functions. 

I should also mention that while Fourier series are technically intended for periodic functions, our function actually does not need to repeat itself within the domain we are calculating it for. Further, how the function behaves outside that domain is irrelevant to the calculation of Fourier series. So the Fourier series can, in fact, match our function perfectly within our domain despite not actually being periodic. Although, this lack of periodicity will cause the Fourier series and the original function to not match outside the domain in which we calculated the series. 

Interestingly, our differential equation for $X_n$ is an example of a type of ODE problem known as a Sturm-Liouville problem. Sturm-Liouville problems, in general, are of the form:[^1]
\begin{align}
\dfrac{d}{dx}\left[p(x) \dfrac{dy}{dx}\right] + q(x)y &= -\lambda r(x) y \label{SLP}
\end{align}
Where $p>0$, $r\geq 0$ and $p$, $q$ and $r$ are well-defined and continuous functions throughout the domain in which we are solving this problem (we will denote this domain as $[a,b]$). There are also boundary conditions imposed as part of a Sturm-Liouville problem. $\lambda$ can take on a set of discrete values, which are called eigenvalues. Each eigenvalue corresponds to a different solution $y$ which is called the corresponding eigenfunction. And yes, as you have probably guessed, the "eigen" prefix is German in origin. I should say that eigenvalues also arise in the study of matrices. Namely, if you have the matrix equation:
\begin{align*}
\mathbf{A} \mathbf{x} &= \lambda \mathbf{x}
\end{align*}
then $\lambda$ is said to be an eigenvalue of $\mathbf{A}$ and $\mathbf{x}$ is an eigenvector of $\mathbf{A}$. 

So in the example we solved for the heat equation, $X_n$ were the eigenfunctions and $\lambda_n$ were the eigenvalues. An interesting quality of eigenfunctions is that they are orthogonal. That means that:
\begin{align*}
\langle y_m | y_n \rangle &= \int_a^b y_m y_n r(x) dx \\
&= c_m \delta_{mn}.
\end{align*}
Where $c_m$ is some constant and $\langle y_m | y_n \rangle$ is called the inner product between $y_m$ and $y_n$. Interestingly, the dot product between two vectors is also referred to as the inner product. Just like the function inner product, the dot product also yields zero if two vectors are orthogonal. There are a lot of parallels between Sturm-Liouville theory and the study of matrices and vectors. 

This orthogonality quality of eigenfunctions is quite useful to us when we are solving linear partial differential equations using separation of variables and at least one of the equations we get from this is a Sturm-Liouville problem. As it often makes finding the coefficients of the series expansion of our solution easier. Another partial differential equation whose solution gives rise to Sturm-Liouville problems is the Schrödinger equation.