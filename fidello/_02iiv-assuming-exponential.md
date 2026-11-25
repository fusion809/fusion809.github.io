### Assuming an exponential solution
Another method we could use for this sort of ODE is only applicable to ODEs with constant coefficients. That is, ODEs of the form:
\begin{align*}
\sum_{i=0}^{n} a_i \dfrac{d^i y}{dx^i} &= 0.
\end{align*}
Where $a_i$ are constant numbers. In this method, we assume that:
\begin{align*}
y &= \sum_{j=1}^{n} C_j e^{r_j x}.
\end{align*}
Where $C_j$ is a constant that depends on our initial conditions and $r_j$ depends on our equation. Susbstituting into our ODE yields:
\begin{align*}
\sum_{i=0}^{n} \sum_{j=1}^{n} a_i C_j r_j^i e^{r_j x} &= 0 \\
\sum_{j=1}^{n} C_j e^{r_j x}\sum_{i=0}^{n} a_i r_j^i &= 0.
\end{align*}
Hence we have to solve:
\begin{align*}
\sum_{i=0}^{n} a_i r_j^i &= 0.  
\end{align*}
for $j\in \{1,2,3,...,n\}$. This is really just a polynomial equation of degree $n$ that we need to solve and $r_j$ is just the $j$th solution to this equation. For $y'=y$, for instance $a_0=-1$ and $a_1=1$ and hence we need to solve:
\begin{align*}
-r_j^0 + r_j &= 0
\implies r_j &= 1.
\end{align*}
Hence our solution is:
\begin{align*}
y(x) &= \sum_{j=1}^1 C_j e^{x} \\
&= C_1 e^x.
\end{align*}
And with $y(0)=1$ we obtain $C_1=1$. Hence $y(x)=e^x$. To see how it works, let us try to solve $y'''+6y''+11y'+6y=0$. Hence our polynomial equation is:
\begin{align*}
r_j^3 + 6r_j^2 + 11r_j + 6 = 0.
\end{align*}
Substituting $r_j=-1$ into the left-hand side of this equation yields $0$. Hence $r_j+1$ is a factor. Using synthetic division (polynomial long division will work equally well) we get:
\begin{align*}
\begin{array}{cc}
    \begin{array}{r} \\ -1 \\ \end{array}
    &
    \begin{array}{|rrrr} \ 
        1 & 6 & 11 & 6 \\
          & -1 & -5  & -6    \\
        \hline 
        1 & 5 & 6 & 0
    \end{array}
\end{array}.
\end{align*}
Hence our polynomial equation becomes:
\begin{align*}
(r_j+1)(r_j^2+5r_j+6) &= 0.
\end{align*}
3 and 2 multiply to give 6 and add to give 5, hence the quadratic factor can be factorized as $(r_j+2)(r_j+3)$:
\begin{align*}
(r_j+1)(r_j+2)(r_j+3) &= 0.
\end{align*}
Applying the null factor law and rearranging yields:
\begin{align*}
r_1 &= -1 \\
r_2 &= -2 \\
r_3 &= -3.
\end{align*}

Hence $y(x) = C_1 e^{-x} + C_2 e^{-2x} + C_3 e^{-3x}$ where $C_1$, $C_2$ and $C_3$ depend on our initial conditions. 