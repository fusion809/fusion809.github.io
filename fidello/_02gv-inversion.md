### Inversion
Square matrices like $\mathbf{A}$ can be inverted. Which is the matrix equivalent to the reciprocal. So if you want to solve the equation: $3x=5$ you would multiply both sides by $\dfrac{1}{3}$ because this is the reciprocal of $3$, the coefficient of $x$. Whereas to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ we could multiply both sides by $\mathbf{A}^{-1}$, the inverse of $\mathbf{A}$ to get:
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