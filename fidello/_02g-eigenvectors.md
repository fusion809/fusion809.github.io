## Eigenvectors and eigenvalues of matrices
If you have the matrix equation:
\begin{align*}
\mathbf{A} \mathbf{x} &= \lambda \mathbf{x}
\end{align*}
then $\lambda$ is said to be an eigenvalue of $\mathbf{A}$ and $\mathbf{x}$ is an eigenvector of $\mathbf{A}$. A square matrix of dimensions $n\times n$ will have $n$ distinct eigenvalues and eigenvectors. Each eigenvector will be of dimensions $n\times 1$. Interestingly, each eigenvector of a matrix is orthogonal to the others. In other words, they are perpendicular to each other &mdash; the angle between them is $90^{\circ}$. This also means that their dot product will be zero.

We solve for $\lambda$ by solving the polynomial equation that results from:
\begin{align*}
\det{\left(\mathbf{A}-\lambda \mathbf{I}\right)} &= 0.
\end{align*}
Where $\det{\left(\mathbf{A}-\lambda \mathbf{I}\right)}$ denotes the determinant of matrix $\mathbf{A}-\lambda \mathbf{I}$ and $\mathbf{I}$ is the identity matrix. Then we figure out $\mathbf{x}$ by solving:
\begin{align*}
(\mathbf{A}-\lambda\mathbf{I})\mathbf{x} = \mathbf{0}.
\end{align*}
When we solve for the eigenvectors, we will find that one unknown is left and this means that the scaling of the vector is arbitrary. Often we choose the scaling to make it so that $\mathbf{x} \cdot \mathbf{x} = 1$. When this condition is satisfied we say that our eigenvectors are normalized. Say, for example, we have the matrix:
\begin{align*}
\mathbf{A} &= \begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}.
\end{align*}
Hence the determinant equation is:
\begin{align*}
\det{\left(\begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}-\lambda \begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}\right)} &= 0 \\
\det{\left(\begin{bmatrix}
1 - \lambda & 2 \\
3 & 4 - \lambda
\end{bmatrix}\right)} &= 0 \\
\begin{vmatrix}
1 - \lambda & 2 \\
3 & 4-\lambda
\end{vmatrix} & = 0 \\
(1-\lambda)(4-\lambda)-2\times 3 &= 0 \\
4 - 4\lambda - 1 \lambda + \lambda^2 -6 &= 0 \\
-2-5\lambda + \lambda^2 &= 0 \\
\lambda^2 -5 \lambda - 2 &= 0
\end{align*}
Next we will use the completing the square technique to solve our quadratic equation:
\begin{align*}
\lambda^2 -5 \lambda + \left(\dfrac{5}{2}\right)^2 - \left(\dfrac{5}{2}\right)^2- 2 &= 0 \\
\left(\lambda-\dfrac{5}{2}\right)^2 - \dfrac{25}{4} - 2 &= 0 \\
\left(\lambda-\dfrac{5}{2}\right)^2 - \dfrac{33}{4} &= 0 \\
\left(\lambda-\dfrac{5}{2}\right)^2 &= \dfrac{33}{4} \\
\lambda - \dfrac{5}{2} &= \pm \dfrac{\sqrt{33}}{2} \\
\lambda &= \dfrac{5}{2} \pm \dfrac{\sqrt{33}}{2} \\
&= \dfrac{5\pm \sqrt{33}}{2}.
\end{align*}
Hence to find the eigenvectors we need to solve:
\begin{align*}
\begin{bmatrix}
1 - \dfrac{5\pm \sqrt{33}}{2} & 2 \\
3 & 4-\dfrac{5\pm \sqrt{33}}{2}
\end{bmatrix} \begin{bmatrix}
a \\
b
\end{bmatrix} &= \begin{bmatrix}
0 \\
0
\end{bmatrix} \\
\begin{bmatrix}
\left(1 - \dfrac{5\pm \sqrt{33}}{2}\right)a + 2b \\
3a + \left(4 - \dfrac{5\pm \sqrt{33}}{2}\right)b 
\end{bmatrix} &= \begin{bmatrix}
0 \\
0
\end{bmatrix}.
\end{align*}
It might seem like we need to solve two equations here, but with eigenvector problems for $2\times2$ matrices the two equations end up containing the same information. So in this case we only need to solve one. 
\begin{align*}
\left(1 - \dfrac{5\pm \sqrt{33}}{2}\right)a + 2b &= 0 \\
b &= -\dfrac{\left(1 - \dfrac{5\pm \sqrt{33}}{2}\right)a}{2} \\
&= \dfrac{\left(\dfrac{5\pm \sqrt{33}}{2}\right)a-1}{2}.
\end{align*}
So we are left to choose what value of $a$ we want. Say we choose $a=1$ for the sake of simplicity:
\begin{align*}
b &= \dfrac{\left(\dfrac{3\pm \sqrt{33}}{2}\right)}{2} \\
&= \dfrac{3\pm \sqrt{33}}{4}.
\end{align*}
Hence our eigenvectors are:
\begin{align*}
\mathbf{x} &= \begin{bmatrix}
1 \\
\dfrac{3\pm \sqrt{33}}{4}
\end{bmatrix} \\
&= \begin{bmatrix}
1 \\
\dfrac{3\pm \sqrt{33}}{4}
\end{bmatrix}
\end{align*}
We could, at this point, scale our eigenvectors so that its inner product with itself is one. This is called normalization.

\begin{align*}
\mathbf{x} \cdot \mathbf{x} &= 1^2 + \left(\dfrac{3\pm \sqrt{33}}{4}\right)^2 \\
&= 1 + \left(\dfrac{3}{4}\right)^2 + \left(\dfrac{\sqrt{33}}{4}\right)^2 \pm 2 \cdot \dfrac{3}{4} \cdot \dfrac{\sqrt{33}}{4} \\
&= 1 + \dfrac{9}{16} + \dfrac{33}{16} \pm \dfrac{3\sqrt{33}}{8} \\
&= 1 + \dfrac{42}{16} \pm \dfrac{3\sqrt{33}}{8} \\
&= \dfrac{58}{16} \pm \dfrac{3\sqrt{33}}{8} \\
&= \dfrac{1}{16} \left(58\pm 6\sqrt{33}\right).
\end{align*}
Hence our normalized eigenvectors, which we will denote as $\mathbf{\hat{x}}$ would be:
\begin{align*}
\mathbf{\hat{x}} &= \dfrac{\mathbf{x}}{\sqrt{\mathbf{x}\cdot \mathbf{x}}} \\
&= \dfrac{1}{\sqrt{\dfrac{1}{16} \left(58\pm 6\sqrt{33}\right)}} \begin{bmatrix}
1 \\
\dfrac{3\pm \sqrt{33}}{4}
\end{bmatrix} \\
&= \dfrac{4}{\sqrt{58\pm 6\sqrt{33}}} \begin{bmatrix}
1 \\
\dfrac{3\pm \sqrt{33}}{4}
\end{bmatrix} \\
&= \dfrac{1}{\sqrt{58\pm 6\sqrt{33}}} \begin{bmatrix}
4 \\
3\pm \sqrt{33}
\end{bmatrix}.
\end{align*}
