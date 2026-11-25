### Determinant
There's one last matrix concept I should probably cover with you; it is called the determinant. It is a single number value you can obtain from the elements of a square matrix. For a matrix of size $2\times 2$ it is defined as:
\begin{align*}
\det\begin{bmatrix}
a & b \\
c & d
\end{bmatrix} &= \begin{vmatrix}
a & b \\
c & d
\end{vmatrix} \\
&= ad-bc. 
\end{align*}
It gets more complex with bigger matrices. For instance, for a $3\times 3$ matrix it is:
\begin{align*}
\det\begin{bmatrix}
A_{1,1} & A_{1,2} & A_{1,3} \\
A_{2,1} & A_{2,2} & A_{2,3} \\
A_{3,1} & A_{3,2} & A_{3,3}
\end{bmatrix} &= \begin{vmatrix}
A_{1,1} & A_{1,2} & A_{1,3} \\
A_{2,1} & A_{2,2} & A_{2,3} \\
A_{3,1} & A_{3,2} & A_{3,3}
\end{vmatrix} \\
&= A_{1,1} \begin{vmatrix}
A_{2,2} & A_{2,3} \\
A_{3,2} & A_{3,3}
\end{vmatrix} - A_{1,2}\begin{vmatrix}
A_{2,1} & A_{2,3} \\
A_{3,1} & A_{3,3}
\end{vmatrix} + A_{1,3}\begin{vmatrix}
A_{2,1} & A_{2,2} \\
A_{3,1} & A_{3,2}
\end{vmatrix} \\
&= A_{1,1}(A_{2,2}A_{3,3}-A_{2,3}A_{3,2}) - A_{1,2}(A_{2,1}A_{3,3}-A{2,3}A_{3,1}) + A_{1,3}(A_{2,1}A_{3,2}-A{2,2}A_{3,1}).
\end{align*}