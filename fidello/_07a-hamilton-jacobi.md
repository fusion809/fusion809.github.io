# Hamilton-Jacobi equation
There is a final formulation of classical mechanics that I should probably cover. It is based around the so called Hamilton-Jacobi equation which is of the form:
\begin{align}
\dfrac{\partial S}{\partial t} + H\left(q_i, \dfrac{\partial S}{\partial q_i}, t\right) &= 0.\label{HJE}
\end{align}
Where $S$ is the action, which is also called Hamilton's principal function. In this approach, we essentially write our Hamiltonian in terms of $p_i$ and $q_i$ and then replace the $p_i$ with $\dfrac{\partial S}{\partial q_i}$. If our Hamiltonian does not explicitly depend on time, which is usually the case, we then use a separation of variables of the form:
\begin{align}
S(q_i, t) &= W(q_i)-Et.\label{SWE}
\end{align}
Where $W$ is called Hamilton's characteristic function or the abbreviated action and $E$ is the constant energy. 

This approach to mechanics actually lead directly to the Schrödinger equation. Schrödinger actually derived his equation from the Hamilton-Jacobi equation. Why? Well, out of the formulations of classical mechanics, the Hamilton-Jacobi equation allows us to describe mechanical systems in a way that is analogous to light waves. The wave-particle duality &mdash; the idea that every quantum particle, including electrons and photons, exhibit both wave and particle properties &mdash; is a fundamental part of quantum mechanics and hence the Hamilton-Jacobi equation provides a particularly natural starting point for the development of the central equation of quantum mechanics. For an explanation for how we can describe a mechanical system in a way analogous to light waves, I suggest you check out [this excellent video](https://www.youtube.com/watch?v=mkpKR_8z6GU) by Abide by Reason. 