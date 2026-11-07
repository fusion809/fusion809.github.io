## Maclaurin and Taylor series
If we have a function $f(x)$, it is usually possible to write an infinite sum (or series) that equals $f(x)$ using $f(x)$'s derivatives. Namely:

\begin{align}
f(x) &= \sum_{n=0}^{\infty} \dfrac{d^n f}{dx^n}\Big|_{x=a} \dfrac{(x-a)^n}{n!} \label{Taylor} \\
&= f(a) + \sum_{n=1}^{\infty}\dfrac{d^n f}{dx^n}\Big|_{x=a} \dfrac{(x-a)^n}{n!} \label{Taylorfirstterm}. \\
\end{align}

Where we are using $\dfrac{d^n f}{dx^n}\Big|_{x=a}$ to denote the $n$th derivative of $f(x)$ evaluated at $x=a$ and $n!$ denotes the factorial of $n$. If you are unfamiliar, well $n!=n\times (n-1) \times (n-2) \times(n-3) ... \times 1$. For instance, $5!=5\times 4 \times 3 \times 2 \times 1=120$. As for $a$ it is an arbitrary number we get to choose. Equation \eqref{Taylor} is known as the Taylor series of $f(x)$. 

If $a=0$, then we get the Maclaurin series:
\begin{align}
f(x) &= \sum_{n=0}^{\infty} \dfrac{d^n f}{dx^n}\Big|_{x=0} \dfrac{x^n}{n!} \label{Maclaurin} \\
&= f(0) + \sum_{n=1}^{\infty} \dfrac{d^n}{dx^n}\Big|_{x=0} \dfrac{x^n}{n!}. \label{Maclaurinfirstterm}
\end{align}

Why do these equations hold? Well, try taking the $n$th derivative of Equations \eqref{Taylor} and \eqref{Maclaurin} and substituting $x=a$ and $x=0$, respectively. You'll find that you get $\dfrac{d^n f}{dx^n}\Big|_{x=a}$ and $\dfrac{d^n f}{dx^n}\Big|_{x=0}$, respectively. Consequently, our series matches our original function in its $n$th derivative for $n\in \mathbb{Z}^{+}$. Likewise, if you substitute in $x=a$ and $x=0$, respectively, into these series, you should get $f(a)$ and $f(0)$, respectively. So we see that this is like a polynomial of infinite order that matches our original function for its value and derivatives value at one point. It turns out that this means it contains enough information to match our original function at all other points, too.  

Now, I know this probably seems pretty abstract to you right now. But let us try to figure out the Maclaurin series for some example functions for which the $n$th derivative is fairly easy to calculate. Say, for instance, we calculate the Maclaurin series for $e^x$. Why $e^x$? Because the $n$th derivative of $e^x$ is just $e^x$, and evaluating this at $x=0$ yields 1. Hence:

\begin{align}
e^x &= \sum_{n=0}^{\infty} \dfrac{x^n}{n!} \label{ex}
\end{align}
is the Maclaurin series of $e^x$. Okay, let us try our hand at another function's Maclaurin series. How about $\cos{x}$? We know the first derivative of $\cos{x}$ is $-\sin{x}$. Hence the second derivative is $-\cos{x}$. The third derivative would be $\sin{x}$. Fourth derivative would be $\cos{x}$. Fifth derivative would be $-\sin{x}$ and sixth derivative would be $-\cos{x}$. Well, I have noticed a pattern here! Seems to me like $(-1)^{\lfloor\frac{n+1}{2}\rfloor}$ (where $\lfloor ... \rfloor$ rounds down decimals, e.g. it would take 3.6 and return 3) will recreate the sign of this function. If $n$ is even, the trigonometric function will be $\cos{x}$ and if it is odd it will be $\sin{x}$. Evaluating these functions at $x=0$, we see that if $n$ is even we will get $1$ from the trigonometric function and 0 if $n$ is odd. Hence:

\begin{align*}
\dfrac{d^nf}{dx^n}\Big|_{x=0} &= \begin{cases}
(-1)^{\lfloor\frac{n+1}{2}\rfloor} & n\hspace{0.1cm}\mathrm{even} \\
0 & n \hspace{0.1cm}\mathrm{odd}.
\end{cases}
\end{align*}
So our sum will essentially just be over even numbers. One way we can denote a sum over even numbers is by replacing $n$ with $2n$ and summing over all numbers. As $2$ times any integer is always an even number. In this event, our $(-1)^{\lfloor\frac{n+1}{2}\rfloor}$ becomes:

\begin{align}
(-1)^{\lfloor\frac{2n+1}{2}\rfloor} &=(-1)^{\lfloor\frac{2n+1}{2}\rfloor} \nonumber \\
&= (-1)^{\lfloor n+\frac{1}{2}\rfloor} \nonumber \\
&= (-1)^{n}. \label{neg1fl}
\end{align}
We arrived at Equation \eqref{neg1fl} simply because $n$ is an integer so adding a half to it should make no difference when we are rounding down. Consequently, our Maclaurin series becomes:

\begin{align*}
\cos{x} &= \sum_{n=0}^{\infty} \dfrac{(-1)^n x^{2n}}{(2n)!}
\end{align*}

So here is some homework for you: try constructing the Maclaurin series for $\sin{x}$. You should arrive at:
\begin{align*}
\sin{x} &= \sum_{n=0}^{\infty} \dfrac{(-1)^nx^{2n+1}}{(2n+1)!}.
\end{align*}

### Classical limit of special relativity's kinetic energy
Maclaurin and Taylor series have broad applications in physics. One application they have is in helping us connect more general and robust laws of physics to older approximations. For instance, according to Einstein the kinetic energy of a particle is given by:
\begin{align*}
T(v) &= mc^2 \dfrac{1}{\sqrt{1-\dfrac{v^2}{c^2}}}-mc^2 \\
&= mc^2\left(\dfrac{1}{\sqrt{1-\dfrac{v^2}{c^2}}}-1\right).
\end{align*}
Using the Maclaurin series, we are going to show that this is well approximated by $T=\dfrac{mv^2}{2}$ when $v << c$. To make things simpler for ourselves, we will consider the Maclaurin series of the function:

\begin{align*}
f(x) &= \dfrac{1}{\sqrt{1-x}} \\
&= (1-x)^{-1/2}.
\end{align*}

(This is as opposed to directly calculating the Maclaurin series of our $T(v)$. If we directly applied the Maclaurin series formula to $T(v)$ the derivatives we would need to calculate would get a bit unruly.) In this case:
\begin{align*}
f'(x) &= -\dfrac{1}{2}(1-x)^{-3/2} \cdot -1 \\
&= \dfrac{(1-x)^{-3/2}}{2} \\
f''(x) &= -\dfrac{3(1-x)^{-5/2}}{4} \cdot -1 \\
&= \dfrac{3(1-x)^{-5/2}}{4} \\
f^{(3)}(x) &= -\dfrac{15(1-x)^{-7/2}}{8} \cdot -1 \\
&= \dfrac{15(1-x)^{-7/2}}{8}
\end{align*}
Where $f^{(3)}(x)$ denotes the 3rd derivative of $f$ with respect to $x$. Note how for the first derivative the coefficient is just $\dfrac{1}{2}$, for the second derivative it is $\dfrac{3}{4}$ and for the third it is $\dfrac{15}{8}$. Seems to me that the denominator is just $2^n$ where $n$ is the order of the derivative and the numerator of the coefficient is $(2n-1) \times (2n-3) \times (2n-5)\times ... \times 1$. Hence:
\begin{align*}
\dfrac{d^nf}{dx^n} &= \dfrac{\left(\displaystyle \prod_{i=1}^n 2i-1\right)(1-x)^{-(2n+1)/2}}{2^{n}} \hspace{0.1cm}\forall n\in \mathbb{Z}^{+}.
\end{align*}
$\displaystyle \prod_{i=1}^n 2i-1$ is a way of writing a multi-term product that is somewhat analogous to the summation notation $\sum_{i=1}^n 2i-1$. In this sum, we add up 2i-1 for $i=1$, $i=2$, $i=3$, and so forth, all the way up until $i=n$. Whereas in the product $\displaystyle \prod_{i=1}^n 2i-1$ we would multiply 2i-1 for $i=1$, $i=2$, $i=3$, and so forth up until $i=n$. Likewise $\mathbb{Z}^{+}$ denotes the set of all positive integers (so $\{1, 2, 3, ...\}$). For instance if $n=3$:
\begin{align*}
\displaystyle \prod_{i=1}^3 2i-1 &= (2(1)-1) \times (2(2)-1) \times (2(3)-1) \\
&= 1 \times 3 \times 5 \\
&= 15.
\end{align*}

Substituting in $x=0$ into $f(x)$ and our derivatives yields:
\begin{align*}
f(0) &= 1 \\
f'(0) &= \dfrac{1}{2} \\
\dfrac{d^nf}{dx^n}\Big|_{x=0} &= \dfrac{\displaystyle \prod_{i=1}^n 2i-1}{2^{n}}.
\end{align*}

In this case, $x=\dfrac{v^2}{c^2}$ and hence:
\begin{align*}
T(v) &= mc^2 \left(1 + \sum_{n=1}^{\infty} \dfrac{\left(\displaystyle \prod_{i=1}^n 2i-1\right)\left(\dfrac{v^2}{c^2}\right)^n}{2^{n}n!}-1\right)\\
&= mc^2 \sum_{n=1}^{\infty} \dfrac{\left(\displaystyle \prod_{i=1}^n 2i-1\right)\left(\dfrac{v^2}{c^2}\right)^n}{2^{n}n!}.
\end{align*}
Hence the first few terms of our Maclaurin series are:
\begin{align*}
T(v) &= mc^2 \left(\dfrac{v^2}{2c^2}+\dfrac{3v^4}{8c^4}+\dfrac{15v^6}{48c^6}+\dots\right) \\
&= \dfrac{mv^2}{2} + \dfrac{3mv^4}{8c^2} + \dfrac{15mv^6}{48c^4} + \dots.
\end{align*}
So what happens when $v$ is much smaller than $c$? Well, we would expect every term with a non-zero power of $c$ in the denominator to be negligible in its contribution to $T$ and this leaves only the first term. Hence we obtain that for $v << c$:
\begin{align*}
T(v) &\approx \dfrac{mv^2}{2}.
\end{align*}
