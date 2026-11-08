## Maclaurin and Taylor series
If we have a function $f(x)$, it is usually possible to write an infinite sum (or series) that equals $f(x)$ using $f(x)$'s derivatives.Namely:

\begin{align}
f(x) &= \sum_{n=0}^{\infty} \dfrac{d^n f}{dx^n}\Big|_{x=a} \dfrac{(x-a)^n}{n!} \label{Taylor} \\
&= f(a) + \sum_{n=1}^{\infty}\dfrac{d^n f}{dx^n}\Big|_{x=a} \dfrac{(x-a)^n}{n!} \label{Taylorfirstterm}. \\
\end{align}

Where we are using $\dfrac{d^n f}{dx^n}\Big|_{x=a}$ to denote the $n$th derivative of $f(x)$ evaluated at $x=a$ and $n!$ denotes the factorial of $n$. If you are unfamiliar, well $n!=n\times (n-1) \times (n-2) \times(n-3) ... \times 1$. For instance, $5!=5\times 4 \times 3 \times 2 \times 1=120$. As for $a$ it is an arbitrary number we get to choose. Equation \eqref{Taylor} is known as the Taylor series of $f(x)$. 

If $a=0$, then we get the Maclaurin series:
\begin{align}
f(x) &= \sum_{n=0}^{\infty} \dfrac{d^n f}{dx^n}\Big|_{x=0} \dfrac{x^n}{n!} \label{Maclaurin} \\
&= f(0) + \sum_{n=1}^{\infty} \dfrac{d^n f}{dx^n}\Big|_{x=0} \dfrac{x^n}{n!}. \label{Maclaurinfirstterm}
\end{align}

Why do these equations hold? Well, try taking the $n$th derivative of Equations \eqref{Taylor} and \eqref{Maclaurin} and substituting $x=a$ and $x=0$, respectively. You'll find that you get $\dfrac{d^n f}{dx^n}\Big|_{x=a}$ and $\dfrac{d^n f}{dx^n}\Big|_{x=0}$, respectively. Consequently, our series matches our original function in its $n$th derivative for $n\in \mathbb{Z}^{+}$. Likewise, if you substitute in $x=a$ and $x=0$ into the Taylor and Maclaurin series, respectively, you should get $f(a)$ and $f(0)$, respectively. So we see that this is like a polynomial of infinite order that matches our original function for its value and derivatives value at one point. It turns out that this means it contains enough information to match our original function at other points, too.

There are exceptions where Taylor series do not accurately compute the function they pertain to. For instance, points where the function is undefined can mess things up. As an example, $\ln{(x+1)}$'s Maclaurin series will only approximate $f(x)$ well for $-1 < x < 1$ as at $x=-1$ the function is undefined. Beyond these points it will fail to converge (meaning it will not equal some finite number as terms for higher $n$ will not get rapidly smaller as $n$ gets bigger).

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

As another example, we will find the Maclaurin series for $\ln{(x+1)}$ so we can see why it only converges well for $|x| < -1$ where $|x|$ is the absolute value of $x$. The first derivative is $\dfrac{1}{x+1}$. Second derivative is $-\dfrac{1}{(x+1)^2}$. Third derivative is $\dfrac{2}{(x+1)^3}$. Fourth derivative is $-\dfrac{6}{(x+1)^4}$. So it seems to me that the pattern is that the numerator is $(-1)^{n-1} (n-1)!$ and the denominator is $(x+1)^n$. Hence:
\begin{align}
\dfrac{d^n}{dx^n} \ln{(x+1)} &= \dfrac{(-1)^{n-1}(n-1)!}{(x+1)^n}, \hspace{0.1cm} n \geq 1.
\end{align}
And:
\begin{align}
\dfrac{d^n}{dx^n} \ln{(x+1)} \Big|_{x=0} &= (-1)^{n-1}(n-1)!, \hspace{0.1cm} n \geq 1.
\end{align}
Hence our Maclaurin series is:
\begin{align}
\ln{(x+1)} &= \ln{1} + \sum_{n=1}^{\infty} \dfrac{(-1)^{n-1}(n-1)! x^n}{n!} \\
&= \sum_{n=1}^{\infty} \dfrac{(-1)^{n-1}x^n}{n}.
\end{align}
If $|x|>1$, we would expect the power of $x$ in our numerator to potentially grow faster than $n$ and hence our terms could become progressively larger as $n\rightarrow \infty$. As for if $x=\pm 1$, well our numerator would be smaller than our denominator, but our denominator would grow rather slowly, which means we would have to calculate a lot of terms to figure out what value it ends up equalling, assuming it does converge. If $x=1$, we get:
\begin{align}
\ln{2} &= \sum_{n=1}^{\infty} \dfrac{(-1)^{n-1}}{n}.
\end{align}
Which is actually true. As this sum does, indeed, equal $\ln{2}$. But because the denominator grows very slowly, we would need a lot of terms in the series in order to get a good approximation for what it equals. If $x=-1$, then we get:
\begin{align}
\ln{0} &= \sum_{n=1}^{\infty} \dfrac{(-1)^{n-1} (-1)^{n}}{n} \\
&= -\sum_{n=1}^{\infty} \dfrac{1}{n}.
\end{align}
This sum does not converge, hence it actually equals $\infty$. Hence why $|x| < 1$ is a requirement for this Maclaurin series to converge, or at least converge at a reasonable pace, although it does also converge, albeit slowly, for $x=1$. 

Maclaurin and Taylor series have broad applications in physics. One application they have is in helping us connect more general and robust laws of physics to older approximations. For more on this, see the section on [the classical limit of special relativity's kinetic energy](#classical_limit_of_special_relativitys_kinetic_energy).