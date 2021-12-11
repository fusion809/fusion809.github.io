@def title="Approximating square roots without a calculator"
@def tags = ["maths", "numerical methods"]

A popular technique for approximating the square root of a number is to use the tangent line approximation. Say we wish to find $a_s=\sqrt{a}$, and we know that the square root of a neighbouring number $b$ is $b_s$, then we would use the approximation: $\sqrt{a} \approx \sqrt{b} + \dfrac{a-b}{2\sqrt{b}}$. This approximation, while often satisfactory, can be substantially off if $b$ and $a$ are relatively quite different. In this case, you can use Newton's method to refine the approximation. To introduce this technique, let us define a function $f(x) = x^2 - a$, hence $f'(x) = 2x$. Newton's method then gives us this scheme for approximating $a_s$:

\begin{eqnarray}
x_{n+1} &=& x_n - \dfrac{f(x_n)}{f'(x_n)} \\\\
&=& x_n - \dfrac{x_n^2 - a}{2x_n} \\\\
&=& x_n - \dfrac{x_n}{2} + \dfrac{a}{2x_n} \\\\
&=& \dfrac{x_n}{2} + \dfrac{a}{2x_n}.
\end{eqnarray}

Where $x_{n}$ is our $n$th Newton's method approximation to $a_s$.

# Example 1, finding $\sqrt{2}$ to five decimal places
A classic example where the tangent line approximation is pretty poor is when we are trying to find $\sqrt{2}$. The tangent line approximation gives:

\begin{eqnarray}
    \sqrt{2} &\approx& \sqrt{1} + \dfrac{2-1}{2\sqrt{1}} \\\\
    &=& 1 + \dfrac{1}{2} \\\\
    &=& \dfrac{3}{2} \\\\
    &=& 1.5.
\end{eqnarray}

When of course $\sqrt{2}$ is to five decimal places $1.41421$, and the above approximation is only accurate to one significant figure (or not even that if you round it to the nearest integer). We can improve on this approximation using Newton's method. Our first iteration is:

\begin{eqnarray}
    x_1 &=& \dfrac{x_0}{2} + \dfrac{2}{2x_0} \\\\
    &=& \dfrac{\dfrac{3}{2}}{2} + \dfrac{2}{2\times \dfrac{3}{2}} \\\\
    &=& \dfrac{3}{4} + \dfrac{2}{3} \\\\
    &=& \dfrac{3\times 3}{4\times 3} + \dfrac{2\times 4}{3\times 4} \\\\
    &=& \dfrac{9}{12} + \dfrac{8}{12} \\\\
    &=& \dfrac{17}{12} \\\\
    &\approx& 1.417.
\end{eqnarray}

Now this is clearly not accurate to five decimal places, so another iteration of Newton's method is required. 

\begin{eqnarray}
    x_2 &=& \dfrac{x_1}{2} + \dfrac{2}{2x_1} \\\\
    &=& \dfrac{\dfrac{17}{12}}{2} + \dfrac{1}{\dfrac{17}{12}} \\\\
    &=& \dfrac{17}{24} + \dfrac{12}{17} \\\\
    &=& \dfrac{17 \times 17}{24 \times 17} + \dfrac{12 \times 24}{17 \times 24} \\\\
    &=& \dfrac{289}{408} + \dfrac{288}{408} \\\\
    &=& \dfrac{577}{408} \\
    &\approx& 1.4142157.
\end{eqnarray}

This is accurate to four decimal places, the fifth is not as we would round it up to 2, not 1. Running Newton's method again yields:

\begin{eqnarray}
    x_3 &=& \dfrac{x_2}{2} + \dfrac{2}{2x_2} \\\\
    &=& \dfrac{\dfrac{577}{408}}{2} + \dfrac{1}{\dfrac{577}{408}} \\\\
    &=& \dfrac{577}{816} + \dfrac{408}{577} \\\\
    &=& \dfrac{577 \times 577}{816 \times 577} + \dfrac{408 \times 816}{577 \times 816} \\\\
    &=& \dfrac{332929}{470832} + \dfrac{332928}{470832} \\\\
    &=& \dfrac{665857}{470832} \\
    &\approx& 1.41421356237469.
\end{eqnarray}

Which is an accurate approximation of $\sqrt{2}$ to 11 decimal places. 

# Example 2, approximating $\sqrt{76}$ to five decimal places
The tangent line approximation yields:

\begin{eqnarray}
    \sqrt{76} &\approx& \sqrt{81} - \dfrac{5}{2\sqrt{81}} \\\\
    &=& 9 - \dfrac{5}{2\times 9} \\\\
    &=& 9 - \dfrac{5}{18} \\\\
    &=& \dfrac{9\times 18}{18} - \dfrac{5}{18} \\\\
    &=& \dfrac{162}{18} - \dfrac{5}{18} \\\\
    &=& \dfrac{157}{18} \\\\
    &\approx& 8.722.
\end{eqnarray}

The actual value of $\sqrt{76}$ is $8.717797887081348$, so this approximation will not suffice. Applying Newton's method yields:

\begin{eqnarray}
    x_1 &=& \dfrac{x_0}{2} + \dfrac{76}{2x_0} \\\\
    &=& \dfrac{157}{18 \times 2} + \dfrac{76}{2 \times \dfrac{157}{18}} \\\\
    &=& \dfrac{157}{36} + \dfrac{38}{\dfrac{157}{18}} \\\\
    &=& \dfrac{157}{36} + \dfrac{38 \times 18}{157} \\\\
    &=& \dfrac{157 \times 157}{36 \times 157} + \dfrac{684 \times 36}{157 \times 36} \\\\
    &=& \dfrac{24649}{5652} + \dfrac{24624}{5652} \\\\
    &=& \dfrac{49273}{5652} \\\\
    &\approx& 8.717799.
\end{eqnarray}

Which is accurate to 5 decimal places. 