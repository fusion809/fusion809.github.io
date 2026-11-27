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
And yes, $\mathbf{x}$ and $\mathbf{b}$ are both column vectors. Typically, in these equations $\mathbf{A}$ is a square matrix. If it is a square matrix, we will get exact solutions to our equation. 

If it has more rows than columns, then we have an overdetermined system and our solution for $\mathbf{x}$ will be approximate. Overdetermined systems arise when we are fitting lines of best fit to bivariate data. In that example, you will notice that some of the scatter plot data will not be exactly on the line of best fit. This is precisely because of the approximate nature of our solution. Underdetermined systems have the opposite issue, as we have more columns than rows. This gives us freedom to choose one or more elements of our solution vector $\mathbf{x}$. 