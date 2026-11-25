### Multiplication
If you multiply two matrices, you multiply the elements along the corresponding row of the first matrix by the elements of the corresponding column of the second matrix and sum them up. For instance, say we multiplied $\mathbf{A}$ by this matrix:
\begin{align*}
\mathbf{B} &= \begin{bmatrix}
5 & 6 \\
7 & 8
\end{bmatrix}
\end{align*}
Then:
\begin{align*}
\mathbf{A} \mathbf{B} &= \begin{bmatrix}
1 \times 5 + 2 \times 7 & 1 \times 6 + 2 \times 8 \\
3 \times 5 + 4 \times 7 & 3 \times 6 + 4 \times 8
\end{bmatrix} \\
&= \begin{bmatrix}
19 & 22 \\
43 & 50
\end{bmatrix}.
\end{align*}
As you can see the first row first column element of the product matrix was obtained by multiplying the first row of $\mathbf{A}$ by the first column of $\mathbf{B}$. Likewise, the second row first column element of the product matrix was obtained by multiplying the second row of $\mathbf{A}$ by the first column of $\mathbf{B}$. 

Due to how we multiply matrices, the matrix product is only defined if we are multiplying by a $m\times n$ matrix by a $n\times p$ matrix. The number of columns of the first matrix must match the number of rows of the second. And the result of this multiplication will be a $m\times p$ matrix. Further, because of this definition of matrix multiplication, the operation is non-commutative. In other words, the order in which we multiply matrices matters. So that if we left-multiply one matrix, say $\mathbf{A}$, by another, say $\mathbf{B}$, this is not the equivalent to right-multiplying $\mathbf{A}$ by $\mathbf{B}$. In other words, in general:
\begin{align*}
\mathbf{A} \mathbf{B} \neq \mathbf{B}\mathbf{A}.
\end{align*}
Although there are some matrices $\mathbf{A}$ and $\mathbf{B}$ that are exceptions to this rule.
