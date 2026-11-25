### Simultaneous equations
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