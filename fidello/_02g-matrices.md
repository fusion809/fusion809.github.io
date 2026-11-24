## Matrices
Matrices have two potential dimensions &mdash; rows and columns. As such the dimensions or size of a matrix is described as being $m\times n$ where $m$ is the number of rows and $n$ is the number of columns. If a matrix has an equal number of rows and columns, it is said to be a square matrix. For instance:
\begin{align*}
\mathbf{A} &= \begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}.
\end{align*}
If we want to refer to a specific element of this matrix, we can use its row number and column number. For instance: $A_{1,1} = 1$ and $A_{2,2}=4$. 

Vectors are a type of matrix; specifically, they are one-dimensional matrices. Either they only have one row (in which case they are a row vector) or one column (in which case they are a column vector). 

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
Due to how we multiply matrices, the matrix product is only defined if we are multiplying by a $m\times n$ matrix by a $n\times p$ matrix. The number of columns of the first matrix must match the number of rows of the second. And the result of this multiplication will be a $m\times p$ matrix. Further, because of this definition of matrix multiplication, the operation is non-commutative. Which means that if we left-multiply one matrix, say $\mathbf{A}$, by another, say $\mathbf{B}$, this is not the equivalent to right-multiply $\mathbf{A}$ by $\mathbf{B}$. In other words, in general:
\begin{align*}
\mathbf{A} \mathbf{B} \neq \mathbf{B}\mathbf{A}.
\end{align*}
Although there are some matrices $\mathbf{A}$ and $\mathbf{B}$ that are exceptions to this rule.

We can write a series of linear simultaneous equations using matrices. Specifically, say we wanted to express:
\begin{align*}
x + 2y &= 4 \\
3x + 4y &= 9
\end{align*}
in matrix form. Then we could rewrite it as:
\begin{align*}
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix} \begin{bmatrix}
x \\
y
\end{bmatrix} &= \begin{bmatrix}
4 \\
9
\end{bmatrix}.
\end{align*}
If we denote:

\begin{align*}
\mathbf{x} &= \begin{bmatrix}
x \\
y
\end{bmatrix} \\
\mathbf{b} &= \begin{bmatrix}
4 \\
9
\end{bmatrix}.
\end{align*}
Then our equation can be written as:
\begin{align}
\mathbf{A}\mathbf{x} &= \mathbf{b}.\label{Axb}
\end{align}
And yes, $\mathbf{x}$ and $\mathbf{b}$ are both column vectors. 

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
So how do we invert a matrix? We use something analogous to the elimination method of solving simultaneous equations. Specifically, a technique called Gaussian elimination. We can also adapt this technique to solve matrix equations like Equation \eqref{Axb}. To invert $\mathbf{A}$ we first augment it with $\mathbf{I}$:
\begin{align*}
[\mathbf{A}|\mathbf{I}] &= \left[
\begin{array}{cc|cc}
1 & 2 & 1 & 0\\
3 & 4 & 0 & 1
\end{array}
\right].
\end{align*}
Our aim in inverting $\mathbf{A}$ is to turn the left of this augmented matrix into the identity matrix. At this point we can read $\mathbf{A}^{-1}$ off the right side of the augmented matrix.

Let us multiply row 1 by 2:
\begin{align*}
\left[
\begin{array}{cc|cc}
2 & 4 & 2 & 0\\
3 & 4 & 0 & 1
\end{array}
\right].
\end{align*}
Next we will subtract row 1 from row 2 and replace row 1 with it:
\begin{align*}
\left[
\begin{array}{cc|cc}
1 & 0 & -2 & 1\\
3 & 4 & 0 & 1
\end{array}
\right].
\end{align*}
Next we will subtract three times row 1 from row 2:
\begin{align*}
\left[
\begin{array}{cc|cc}
1 & 0 & -2 & 1\\
0 & 4 & 6 & -2
\end{array}
\right].
\end{align*}
Next we will divide row 2 by four:
\begin{align*}
\left[
\begin{array}{cc|cc}
1 & 0 & -2 & 1\\
0 & 1 & \dfrac{3}{2} & -\dfrac{1}{2}
\end{array}
\right].
\end{align*}
Hence:
\begin{align*}
\mathbf{A}^{-1} &= \begin{bmatrix}
-2 & 1 \\
\dfrac{3}{2} & -\dfrac{1}{2}
\end{bmatrix}.
\end{align*}

If we want to solve Equation \eqref{Axb}, we use a similar procedure. Except we put $\mathbf{b}$ in place of $\mathbf{I}$ when we form the augmented matrix. The aim is still to convert the $\mathbf{A}$ to $\mathbf{I}$:
\begin{align*}
\left[
\begin{array}{cc|c}
1 & 2 & 4\\
3 & 4 & 9
\end{array}
\right].
\end{align*}
We will do the same row operations. First we will multiply row 1 by 2 and replace row 1 with it. We will denote this operation as $2R_1 \rightarrow R_1'$:
\begin{align*}
\left[
\begin{array}{cc|c}
2 & 4 & 8\\
3 & 4 & 9
\end{array}
\right].
\end{align*}
$R_2 - R_1 \rightarrow R_1'$:
\begin{align*}
\left[
\begin{array}{cc|c}
1 & 0 & 1\\
3 & 4 & 9
\end{array}
\right].
\end{align*}
$R_2 - 3R_1 \rightarrow R_2'$:
\begin{align*}
\left[
\begin{array}{cc|c}
1 & 0 & 1\\
0 & 4 & 6
\end{array}
\right].
\end{align*}
$\dfrac{1}{4}R_2 \rightarrow R_2'$:
\begin{align*}
\left[
\begin{array}{cc|c}
1 & 0 & 1\\
0 & 1 & \dfrac{3}{2}
\end{array}
\right].
\end{align*}
Hence:
\begin{align*}
\mathbf{x} &= \dfrac{bmatrix}
1 \\
\dfrac{3}{2}
\end{bmatrix}.
\end{align*}