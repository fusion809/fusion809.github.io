#### Sturm-Liouville problems
Interestingly, our differential equation for $X_n$ is an example of a type of ODE problem known as a Sturm-Liouville problem. Sturm-Liouville problems, in general, are of the form:[^1]
\begin{align}
\dfrac{d}{dx}\left[p(x) \dfrac{dy}{dx}\right] + q(x)y &= -\lambda r(x) y \label{SLP}
\end{align}
Where $p>0$, $r\geq 0$ and $p$, $q$ and $r$ are well-defined and continuous functions throughout the domain in which we are solving this problem (we will denote this domain as $[a,b]$). There are also boundary conditions imposed as part of a Sturm-Liouville problem. $\lambda$ can take on a set of discrete values, which are called eigenvalues. Each eigenvalue corresponds to a different solution $y$ which is called the corresponding eigenfunction. And yes, as you have probably guessed, the "eigen" prefix is German in origin. I should say that eigenvalues also arise in the study of matrices. Namely, if you have the matrix equation:
\begin{align*}
\mathbf{A} \mathbf{x} &= \lambda \mathbf{x}
\end{align*}
then $\lambda$ is said to be an eigenvalue of $\mathbf{A}$ and $\mathbf{x}$ is an eigenvector of $\mathbf{A}$. We solve for $\lambda$ by solving the polynomial equation that results from:
\begin{align*}
\det{\left(\mathbf{A}-\lambda \mathbf{I}\right)} &= 0.
\end{align*}
Where $\det{\mathbf{B}}$ denotes the determinant of matrix $\mathbf{B}$ and $\mathbf{I}$ is the identity matrix. Then we figure out $\mathbf{x}$ by solving:
\begin{align*}
(\mathbf{A}-\lambda\mathbf{I})\mathbf{x} = \mathbf{0}.
\end{align*}
So in the example we solved for the heat equation, $X_n$ were the eigenfunctions and $\lambda_n$ were the eigenvalues. An interesting quality of eigenfunctions is that they are orthogonal. That means that:
\begin{align*}
\langle y_m | y_n \rangle &= \int_a^b y_m y_n r(x) dx \\
&= c_m \delta_{mn}.
\end{align*}
Where $c_m$ is some constant and $\langle y_m | y_n \rangle$ is called the inner product between $y_m$ and $y_n$. Interestingly, the dot product between two vectors is also referred to as the inner product. Just like the function inner product, the dot product also yields zero if two vectors are orthogonal. There are a lot of parallels between Sturm-Liouville theory and the study of matrices and vectors. 

This orthogonality quality of eigenfunctions is quite useful to us when we are solving linear partial differential equations using separation of variables and at least one of the equations we get from this is a Sturm-Liouville problem. As it often makes finding the coefficients of the series expansion of our solution easier. Another partial differential equation whose solution gives rise to Sturm-Liouville problems is the Schrödinger equation.