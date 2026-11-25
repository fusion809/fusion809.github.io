### Inversion
Most square matrices like $\mathbf{A}$ can be inverted. This is the matrix equivalent to the reciprocal. So if you want to solve the equation: $3x=5$ you would multiply both sides by $\dfrac{1}{3}$ because this is the reciprocal of $3$, the coefficient of $x$. Whereas to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ we could multiply both sides by $\mathbf{A}^{-1}$, the inverse of $\mathbf{A}$ to get:
\begin{align*}
\mathbf{A}^{-1}\mathbf{A} \mathbf{x} &= \mathbf{A}^{-1}\mathbf{b} \\
\mathbf{I}\mathbf{x} &= \mathbf{A}^{-1}\mathbf{b} \\
\mathbf{x} &= \mathbf{A}^{-1}\mathbf{b}.
\end{align*}
What is $\mathbf{I}$? It is the so called identity matrix. It is essentially the square matrix equivalent to $1$. As if you multiply it by any matrix with the required dimensions to multiply with it, multiplication by it should not change the matrix at all. Its elements are ones down the main diagonal and zero in all other elements. For instance, the $2\times 2$ identity matrix is:
\begin{align*}
\mathbf{I}_2 &= \begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}.
\end{align*}
I should also mention that some square matrices have no inverse. This is usually the case when one or more rows can be written as linear combinations of the other rows. For instance:
\begin{align*}
\mathbf{A} &= \begin{bmatrix}
1 & 2 & 3 \\
2 & 4 & 6 \\
4 & 5 & 6
\end{bmatrix}
\end{align*}
is singular because row 2 is row 1 times two. 

It is also possible for a matrix to be invertible, but for inverting the matrix to not be a great idea for another reason. Beyond invertibility, when dealing with a square matrix we also often care about its *condition* number. A matrix with a large, yet finite condition number may be technically invertible but even small errors (such as due to floating point arithmetic errors) in the elements of that matrix can massively change the contents of its inverse. Whereas a matrix with a condition number close to one should not suffer this problem at all. I should also mention that singular matrices have an infinite condition number.