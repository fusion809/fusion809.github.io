## Complex numbers
Complex numbers are numbers that involve both a real and an imaginary part. That is, they are of the form $a+bi$ where $a$ and $b$ are both real numbers and $i$ is the imaginary number. $a$ is called the real part of the number and $b$ is called the imaginary part of the number. Importantly, $a$ and $b$ can be zero, consequently all real numbers are technically complex numbers. We can display complex numbers in the complex plane (also known as the Argand plane) as vectors. Specifically, the $y$ component of this vector is the imaginary part of our number and the $x$ component is the real part of the number. 

~~~
<figure id="complex">
<img src="Complex plane.svg" width="300px" float="right">
<figcaption>Number 2+3i displayed in the complex plane.</figcaption>
</figure>
~~~

Figure \figref{complex} shows the number $2+3i$ displayed in the complex plane. The length of the vector, $l$, can be calculated using the Pythagorean theorem as $l=\sqrt{2^2+3^2}=\sqrt{13}$. $l$ can also be called the modulus of our complex number. As for $\theta$, the angle from the positive $x$-axis measure counterclockwise it can be found using the definition of the tangent function as $\theta =\tan^{-1}{\dfrac{3}{2}}$. Note also that $3=l\sin{\theta}$ and $2=l\cos{\theta}$ by the definition of these trigonometric functions. 

Interestingly, these two qualities, the length of the complex number's vector and its angle, measured counterclockwise from the positive $x$-axis, can be used to uniquely define it. To see how, I must first share with you an interesting quality of the exponential function. Namely, that:
\begin{align*}
e^{iz} &= \cos{z} + i\sin{z}. 
\end{align*}
Why is this the case? Replacing $x$ with $iz$ in Equation \eqref{ex} yields:
\begin{align*}
e^{iz} &= \sum_{n=0}^{\infty} \dfrac{(iz)^n}{n!}.
\end{align*}
You may be thinking, so what? Well, powers of $i$ have interesting properties. Namely, $i^0=1$, $i^1=i$, $i^2=-1$, $i^3=-i$ and $i^4=1$. Generally, $i^{4n}$ (where $n$ is an integer) is equal to $1$, $i^{4n+1}$ is $i$, $i^{4n+2}$ is $-1$ and $i^{4n+3}$ is $-i$. Focusing on even powers, we find that $i^{2n}=(-1)^n$ (where makes sense as $i^2=-1$ and $i^{2n}=(i^2)^n$). Whereas odd powers yields $i^{2n+1}=(-1)^n i$. Consequently:
\begin{align*}
e^{iz} &= \sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n}}{(2n)!} + \dfrac{(-1)^n iz^{2n+1}}{(2n+1)!} \\
&= \sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n}}{(2n)!} + \sum_{n=0}^{\infty} \dfrac{(-1)^n iz^{2n+1}}{(2n+1)!} \\
&= \sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n}}{(2n)!} + i\sum_{n=0}^{\infty} \dfrac{(-1)^n z^{2n+1}}{(2n+1)!}.
\end{align*}
The first of these sums we can recognize as $\cos{z}$ whereas the second is $\sin{z}$. Hence:
\begin{align*}
e^{iz} &= \cos{z} + i\sin{z}.
\end{align*}
An interesting thing to note that $e^{iz}$ is that its modulus is one. How do I know this? Well, the modulus is defined as just the length of the complex number's vector in the complex plane, right? So for $e^{iz}$ that would be:
\begin{align*}
l &= \sqrt{\cos{z}^2+\sin^2{z}}\\
&= \sqrt{1} \\
&= 1.
\end{align*}
Here I used the trigonometric identity:
\begin{align}
\cos^2{z}+\sin^2{z}=1 \hspace{0.1cm} \forall z \in \mathbb{R}.
\end{align}
Where $\mathbb{R}$ is the set of real numbers. Okay, so getting back to how $l$ and $\theta$ can uniquely define a complex number. Well, interestingly you can write any complex number in the form $le^{i\theta}$. Why? Well, recall that $e^{i\theta}$ is just $\cos{\theta}+i\sin{\theta}$. Looking at Figure \figref{complex} we can clearly see that the $x$ component of our vector is just $l\cos{\theta}$, right? And we can also see that the $y$ component is just $l\sin{\theta}$ and if we expand our $le^{i\theta}$ we get:
\begin{align*}
le^{i\theta} &= l\cos{\theta} + i\sin{\theta}.
\end{align*}
Hence $le^{i\theta} = x+iy$. For instance, for the complex number $2+3i$ we can rewrite it as:
\begin{align*}
2+3i &= \sqrt{13}e^{i\tan^{-1}{\frac{3}{2}}}.
\end{align*}
### Complex analysis
I have definitely not covered everything within the field of complex numbers. There's a field of maths known as complex analysis where we combine calculus and complex numbers. In it, we learn, among other things, how to integrate along curves (or contours) in the complex plane and that there is an extension of Taylor series to the complex plane called the Laurent series. In it, we write:
\begin{align}
f(z) &= \sum_{n=-\infty}^{\infty} a_n (z-z_0)^n
\end{align}
where:
\begin{align}
a_n &= \dfrac{1}{2\pi i} \oint_C \dfrac{f(z)}{(z-z_0)^{n+1}} dz.
\end{align}
This is a contour integral over a curve $C$ that includes the point $z=z_0$. $z_0$ is analogous to the $a$ in our Taylor series as it is arbitrary. There's a formula called Cauchy's differentiation formula that says that:
\begin{align}
f^{(n)}(z_0) &= \dfrac{n!}{2\pi i} \oint_C \dfrac{f(z)}{(z-z_0)^{n+1}} dz, \hspace{0.1cm}n \in \{0, 1, 2, ...\}.
\end{align}
Where $f^{(n)}(z_0)$ is the $n$th derivative of $f(z)$ evaluated at $z=z_0$. Which helps us see how the Laurent series is connected to the Taylor series as for non-negative $n$, $a_n = \dfrac{f^{(n)}(z_0)}{n!}$. As for negative $n$, $a_n=0$ if $f(z)$ is defined everywhere in the complex plane. Why do I mention this? Well, I personally consider complex analysis one of the most beautiful fields of mathematics and I thought I would mention it in case you are ever in a position to take a complex analysis class while at university. 