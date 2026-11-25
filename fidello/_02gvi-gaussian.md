### Gaussian elimination
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
\mathbf{x} &= \begin{bmatrix}
1 \\
\dfrac{3}{2}
\end{bmatrix}.
\end{align*}
