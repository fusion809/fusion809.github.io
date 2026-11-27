##### Eigenvalue problems
An important issue I should probably address is how do we solve Sturm-Liouville and other differential equations with eigenvalues? Well, the answer is that we construct an eigenvector equation of the form:
\begin{align*}
\mathbf{A} \mathbf{u} = \lambda \mathbf{R}\mathbf{u}.
\end{align*}
Where $\mathbf{u}$ is our vector of values of our function at our grid values and $\mathbf{R}=\mathbf{Diag}(r(x_i))$. These eigenvalue matrix problems can be solved by first multiplying by $\mathbf{R}^{-1}$ before applying the techniques previously mentioned. Typically, we omit the endpoints from our matrices to impose our boundary conditions when those boundary conditions are of the form $u(a)=u(b)=0$. We can impose other boundary conditions by including them in the first and final row of our matrix and if they are zero, the corresponding rows of $\mathbf{R}$ should be zero, too. 

Then we use established techniques for finding eigenvalues and eigenvectors of matrices. We can construct such an equation when we apply finite difference, spectral, finite element and spectral element methods to our linear eigenvalue differential equation problem. 

For example, say we have:
\begin{align}
u'' - xu &= \lambda u \label{airyeig}\\
u(0) &= u(\infty) && = 0. \label{airyeigbc}
\end{align}
Rearranging yields:
\begin{align*}
u'' - (x+\lambda)u &= 0
\end{align*}
If we define a new variable, say $y=x+\lambda$ then:
\begin{align*}
dy &= \dfrac{dy}{dx}dx \\
&= dx.
\end{align*}
Therefore:
\begin{align*}
\dfrac{d}{dx} &= \dfrac{d}{dy}.
\end{align*}
Hence our equation would become, if we rewrite it in terms of $y$:
\begin{align*}
\dfrac{d^2u}{dy^2} - yu &= 0.
\end{align*}
Which is the Airy equation, hence $u(y) = C_1 \mathrm{Ai}(y) + C_2 \mathrm{Bi}(y)$. Hence rewriting $u$ to be in terms of $x$ yields:
\begin{align*}
u(x) &= C_1 \mathrm{Ai}(x+\lambda) + C_2 \mathrm{Bi}(x+\lambda).
\end{align*}
Only $\mathrm{Ai}(x)$ goes to $0$ as $x\rightarrow\infty$. Its derivative also goes to $0$ in this limit. Imposing our other boundary condition, $u'(0)=0$, yields:
\begin{align*}
u'(x) &= C_1 \mathrm{Ai}'(x+\lambda) \\
u'(0) &= C_1 \mathrm{Ai}'(\lambda) \\
&= 0.
\end{align*}
Which implies that $\lambda_n$ will be the turning points of the Airy $\mathrm{Ai}(x)$ function. 
\codelisting{airyeigen}{Solution to Equations \eqref{airyeig} and \eqref{airyeigbc}.}{./code/airyeigen}{fidello/code/airyeigen.jl}
\output{./code/airyeigen}
\figure{airyeigenfig}{right}{600px}{/airyeigen.svg}{First eigenfunction for our Airy Sturm-Liouville problem.}