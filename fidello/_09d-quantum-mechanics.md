## Quantum mechanics
We have seen in this article that there are multiple different formulations of classical mechanics. Similarly, there are multiple different formulations of quantum mechanics. 

### Path-integral formulation
One formulation of quantum mechanics is the path-integral formulation developed by Richard Feynmann. In it, we expect our particle or wave to explore all possible paths and we integrate (or sum) over all those paths in order to compute the probability amplitude of our particle. The exploration of all possible paths is a consequence of the randomness of quantum mechanics. Specifically, we compute the integral:

\begin{align}
K(\mathbf{x}_f, t_f; \mathbf{x}_0, t_0) &= \int_{\mathbf{x}(t_0)=\mathbf{x}_0}^{\mathbf{x}(t_f)=\mathbf{x}_f} e^{\frac{iS(\mathbf{x})}{\hbar}} \mathcal{D}\mathbf{x}.
\end{align}

Where $K$ is the propagator, which specifies the probability amplitude for travelling from $\mathbf{x}_0$ at time $t_0$ to $\mathbf{x}_f$ at time $t_f$. $S$ is the action of the particle, $\hbar$ is the reduced Planck's constant, $i$ is the imaginary number and $\mathcal{D}\mathbf{x}$ refers to the fact we are integrating over all possible paths. The probability amplitude is a quantity that, if you take the absolute value of, and square this quantity you get the probability of the particle travelling from $\mathbf{x}_i$ at time $t_i$ to $\mathbf{x}_f$ at time $t_f$. 

One major strength of the path-integral formulation of quantum mechanics is that it provides an explanation for how quantum mechanics gives rise to classical mechanics on larger scales. 

#### Complex number geometry
So I will need to explain some geometry of complex numbers with you before I can show you how this formulation of quantum mechanics gives rise to classical mechanics on larger scales. $i$, as you know is equal to $\sqrt{-1}$. Complex numbers are numbers that involve both a real and an imaginary part. That is, they are of the form $a+bi$ where $a$ and $b$ are both real numbers and $i$ is the imaginary number. $a$ is called the real part of the number and $b$ is called the imaginary part of the number. Importantly, $a$ and $b$ can be zero, consequently all real numbers are technically complex numbers. We can display complex numbers in the complex plane (also known as the Argand plane) as vectors. Specifically, the $y$ component of this vector is the imaginary part of our number and the $x$ component is the real part of the number. 

~~~
<figure id="complex">
<img src="Complex plane.svg" width="300px" float="right">
<figcaption>Number 2+3i displayed in the complex plane.</figcaption>
</figure>
~~~

Figure \figref{complex} shows the number $2+3i$ displayed in the complex plane. The length of the vector, $l$, can be calculated using the Pythagorean theorem as $l=\sqrt{2^2+3^2}=\sqrt{13}$. $l$ can also be called the modulus of our complex number. As for $\theta$, the angle from the positive $x$-axis measure counterclockwise it can be found using the definition of the tangent function as $\theta =\tan^{-1}{\dfrac{3}{2}}$. Note also that $3=l\sin{\theta}$ and $2=l\cos{\theta}$ by the definition of these trigonometric functions. 

Interestingly, these two qualities, the length of the complex number's vector and its angle, measured counterclockwise from the positive $x$-axis, can be used to uniquely define it. To see how, I must first share with you an interesting quality of the exponential function. Namely, that:
\begin{align}
e^{iz} &= \cos{z} + i\sin{z}. 
\end{align}
Why is this the case? Well, I can explain it to you, but it does require knowledge of a field of mathematics known as Taylor and Maclaurin series. See, if we have a function $f(x)$, it is usually possible to write an infinite sum (or series) that equals $f(x)$ using $f(x)$'s derivatives. Namely:

\begin{align}
f(x) &= \sum_{n=0}^{\infty} \dfrac{d^n f}{dx^n}\Big|_{x=a} \dfrac{(x-a)^n}{n!}. \label{Taylor}
\end{align}

Where we are using $\dfrac{d^n f}{dx^n}\Big|_{x=a}$ to denote the $n$th derivative of $f(x)$ evaluated at $x=a$ and $n!$ denotes the factorial of $n$. As for $a$ it is an arbitrary number we get to choose. Equation \eqref{Taylor} is known as the Taylor series of $f(x)$. If $a=0$, then we get the Maclaurin series:
\begin{align}
f(x) &= \sum_{n=0}^{\infty} \dfrac{d^n f}{dx^n}\Big|_{x=0} \dfrac{x^n}{n!}. \label{Maclaurin}
\end{align}
Now, I know this probably seems pretty abstract to you right now. But let us try to figure out the Maclaurin series for some example functions for which the $n$th derivative is fairly easy to calculate. Say, for instance, we calculate the Maclaurin series for $e^x$. Why $e^x$? Because the $n$th derivative of $e^x$ is just $e^x$, and evaluating this at $x=0$ yields 1. Hence:

\begin{align}
e^x &= \sum_{n=0}^{\infty} \dfrac{x^n}{n!}
\end{align}
is the Maclaurin series of $e^x$. Okay, let us try our hand at another function's Maclaurin series. How about $\cos{x}$? We know the first derivative of $\cos{x}$ is $-\sin{x}$. Hence the second derivative is $-\cos{x}$. The third derivative would be $\sin{x}$. Fourth derivative would be $\cos{x}$. Fifth derivative would be $-\sin{x}$ and six derivative would be $-\cos{x}$. Well, I have noticed a pattern here! Seems to me like $(-1)^{\lfloor\frac{n+1}{2}\rfloor}$ (where $\lfloor ... \rfloor$ rounds down decimals, e.g. it would take 3.6 and return 3) will recreate the sign of this function. If $n$ is even, the trigonometric function will be $\cos{x}$ and if it is odd it will be $\sin{x}$. Evaluating these functions at $x=0$, we see that if $n$ is even we will get $1$ from the trigonometric function and 0 if $n$ is odd. Hence:

\begin{align}
\dfrac{d^nf}{dx^n}\Big|_{x=0} &= \begin{cases}
(-1)^{\lfloor\frac{n+1}{2}\rfloor} & n\hspace{0.1cm}\mathrm{even} \\
0 & n \hspace{0.1cm}\mathrm{odd}.
\end{cases}
\end{align}
So our sum will essentially just be over even numbers. One way we can denote a sum over even numbers is by replacing $n$ with $2n$ and summing over all numbers. As $2$ times any integer is always an even number. In this even our $(-1)^{\lfloor\frac{n+1}{2}\rfloor}$ becomes:

\begin{align}
(-1)^{\lfloor\frac{2n+1}{2}\rfloor} &=(-1)^{\lfloor\frac{2n+1}{2}\rfloor} \nonumber \\
&= (-1)^{\lfloor{n+\frac{1}{2}\rfloor}} \nonumber \\
&= (-1)^{n}. \label{neg1fl}
\end{align}
We arrived at Equation \eqref{neg1fl} simply because $n$ is an integer so adding a half to it should make no difference when we are rounding down. Consequently, our Maclaurin series becomes:

\begin{align}
\cos{x} &= \sum_{n=0}^{\infty} \dfrac{(-1)^n x^{2n}}{(2n)!}
\end{align}

So here is some homework for you: try constructing the Maclaurin series for $\sin{x}$ you should arrive at:
\begin{align}
\sin{x} &= \sum_{n=0}^{\infty} \dfrac{(-1)^nx^{2n+1}}{(2n+1)!}.
\end{align}

Now this brings to me why I brought up the Maclaurin and Taylor series in the first place. As replacing $x$ with $iz$ in the Maclaurin series for $e^x$ yields:
\begin{align}
e^{iz} &= \sum_{n=0}^{\infty} \dfrac{(iz)^n}{n!}.
\end{align}
You may be thinking, so what? Well, powers of $i$ have interesting properties. Namely, $i^0=1$, $i^1=i$, $i^2=-1$, $i^3=-i$ and $i^4=1$. Generally, $i^{4n}$ (where $n$ is an integer) is equal to $1$, $i^{4n+1}$ is $i$, $i^{4n+2}$ is $-1$ and $i^{4n+3}$ is $-i$. Focusing on even powers, we find that $i^{2n}=(-1)^n$. Whereas odd powers yields $i^{2n+1}=(-1)^n i$. Consequently:
\begin{align}
e^{iz} &= \sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n}}{(2n)!} + \dfrac{(-1)^n iz^{2n+1}}{(2n+1)!} \\
&= \sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n}}{(2n)!} + \sum_{n=0}^{\infty} \dfrac{(-1)^n iz^{2n+1}}{(2n+1)!} \\
&= \sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n}}{(2n)!} + i\sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n+1}}{(2n+1)!}.
\end{align}
The first of these sums we can recognize as $\cos{z}$ whereas the second is $\sin{z}$. Hence:
\begin{align}
e^{iz} &= \cos{z} + i\sin{z}.
\end{align}
An interesting thing to note that $e^{iz}$ is that its modulus is one. How do I know this? Well, the modulus is defined as just the length of the complex number's vector in the complex plane, right? So for $e^{iz}$ that would be:
\begin{align*}
l &= \sqrt{\cos{z}^2+\sin^2{z}}
&= \sqrt{1} \\
&= 1.
\end{align*}
Here I used the trigonometric identity that $\cos^2{z}+\sin^2{z}=1$ for all values of $z$. Okay, so getting back to how $l$ and $\theta$ can uniquely define a complex number. Well, interestingly you can write any complex number in the form $le^{i\theta}$. Why? Well, recall that $e^{i\theta}$ is just $\cos{\theta}+i\sin{\theta}$. Looking at Figure \figref{complex} we can clearly see that the $x$ component of our vector is just $l\cos{\theta}$, right? And we can also see that the $y$ component is just $l\sin{\theta}$ and if we expand our $le^{i\theta}$ we get:
\begin{align}
le^{i\theta} &= l\cos{\theta} + i\sin{\theta}.
\end{align}
Hence $le^{i\theta} = x+iy$. For instance, for the complex number $2+3i$ we can rewrite it as:
\begin{align}
2+3i &= \sqrt{13}e^{i\tan^{-1}{\frac{3}{2}}}.
\end{align}

#### Relation to classical mechanics
So, with this background information in mind, we now know that $e^{\frac{iS}{\hbar}}$ is like a 2D unit vector as it points in particular direction in the $xy$ plane and has unit length. Specifically, it points in the direction $\dfrac{S}{\hbar}$ radians counterclockwise from the positive $x$ axis. Each integer multiple of $2\pi$, which corresponds to a complete revolution around the Cartesian plane, that is present in $\dfrac{S}{\hbar}$ makes no difference to its overall direction. What matters is what is left after subtracting these integer multiples of $2\pi$. 

Interestingly, when we integrate over all possible paths on classical mechanics scales we find that the $e^{\frac{iS}{\hbar}}$, for most paths, cancel each other out. To see why, we must first start with the fact that $\hbar$ is a very tiny number &mdash; it is approximately $1.055\times 10^{-34}J\cdot s$. Consequently, $\dfrac{S}{\hbar}$ is a massive number when we are analysing classical mechanical systems. As typically in these scenarios $S$ is $\gt0.01J\cdot s$. Consequently, even a relatively small difference in $S$, can lead to a majorly different direction for $e^{\frac{iS}{\hbar}}$.

This leads to the actions of nearby paths often being in opposite directions and cancelling each other when we calculate the path integral. The one exception is the path corresponding to the stationary action. Remember the equation $\delta S=0$ that we solve to find the path of our classical mechanical system? It says that a tiny change to the stationary action path should make no change (to first order) to the action and consequently, when we integrate $e^{\frac{iS}{\hbar}}$ for that path and its neighbours, the stationary action path term is the only one that survives. Consequently, in the classical limit our object only explores the path that extremizes the action, and as we have seen an extremized action leads to the equations of motion in classical mechanics.

### Operator formalism
Another connection between what we have covered and quantum mechanics is in the operator formalism of quantum mechanics. In it, we represent the state of our physical system using a function called the wavefunction which is denoted as $\psi$, which can also thought of as a vector, as it is a complex number. We can calculate the wavefunction using Schrödinger's equation:

\begin{align}
i\hbar \dfrac{\partial \psi}{\partial t} &= \left(-\dfrac{\hbar^2}{2m}\Delta + V(\mathbf{x})\right)\psi.\label{Schrod}
\end{align}

$V(\mathbf{x})$ is the potential energy and $\Delta$ is the Laplacian which is essentially just $\dfrac{\partial^2}{\partial x^2} + \dfrac{\partial^2}{\partial y^2}+\dfrac{\partial^2}{\partial z^2}$. The right-hand side of Equation \eqref{Schrod} can also be represented by $\hat{H}\psi$ where $\hat{H}$ is the Hamiltonian operator and calculates, from the wavefunction, the total energy of our system.
