#### Sturm-Liouville problems
Interestingly, our differential equation for $X_n$ is an example of a type of ODE problem known as a Sturm-Liouville problem. Sturm-Liouville problems, in general, are of the form \cite{trench_132_2024}:
\begin{align}
\dfrac{d}{dx}\left[p(x) \dfrac{dy}{dx}\right] + q(x)y &= -\lambda r(x) y. \label{SLP}
\end{align}
Where $p>0$, $r\geq 0$. There are also boundary conditions imposed as part of a Sturm-Liouville problem, namely \cite{lebl_51_2025}:
\begin{align*}
\alpha_1 y(a) -\alpha_2 y'(a) &= 0 \\
\beta_1 y(b) + \beta_2 y'(b) &= 0.
\end{align*}
Where $\alpha_1$ and $\alpha_2$ cannot both be simultaneously zero and nor can $\beta_1$ and $\beta_2$. $\lambda$ can take on a set of discrete values, which are called eigenvalues. Each eigenvalue corresponds to a different solution $y$ which is called the corresponding eigenfunction. We usually order the eigenvalues, and the corresponding eigenfunctions, by the absolute value of the eigenvalue in ascending order. So the first eigenvalue is the one with the smallest absolute value and the first eigenfunction is the eigenfunction corresponding to the eigenvalue with the lowest absolute value. 

So in the example we solved for the heat equation, $X_n$ were the eigenfunctions and $\lambda_n$ were the eigenvalues. An interesting quality of eigenfunctions is that they are orthogonal (similarly to the eigenvectors of a matrix). That means that:
\begin{align*}
\langle y_m | y_n \rangle &= \int_a^b y_m y_n r(x) dx \\
&= c_m \delta_{mn}.
\end{align*}
Where $c_m$ is some constant and $\langle y_m | y_n \rangle$ is called the inner product between $y_m$ and $y_n$. Interestingly, the dot product between two vectors is also referred to as the inner product. Just like the function inner product, the dot product also yields zero if two vectors are orthogonal. There are a lot of parallels between Sturm-Liouville theory and the study of matrices and vectors. 

This orthogonality quality of eigenfunctions is quite useful to us when we are solving linear partial differential equations using separation of variables and at least one of the equations we get from this is a Sturm-Liouville problem. As it often makes finding the coefficients of the series expansion of our solution easier. Another partial differential equation whose solution gives rise to Sturm-Liouville problems is the Schrödinger equation.
