### Addition/subtraction
Addition and subtraction between matrices is element-wise. So:
\begin{align*}
\mathbf{C} &= \mathbf{A}-\mathbf{B}.
\end{align*}
Means that $\mathbf{C}$ has elements that are equal to the corresponding element of $\mathbf{A}$ minus the corresponding element in $\mathbf{B}$. For instance:
\begin{align*}
\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix} - \begin{bmatrix}
0 & 1 \\
5 & 6
\end{bmatrix} &= \begin{bmatrix}
1-0 & 2-1 \\
3-5 & 4-6
\end{bmatrix} \\
&= \begin{bmatrix}
1 & 1 \\
-2 & -2
\end{bmatrix}.
\end{align*}