## Maclaurin and Taylor series
If we have a function $f(x)$, it is usually possible to write an infinite sum (or series) that equals $f(x)$ using $f(x)$'s derivatives. Namely:

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