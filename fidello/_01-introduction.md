\figure{dp}{right}{600px}{/Double_pendulum.gif}{Double pendulum animation.}

In this article, I explain to Fidello from Foci's Discord server some prerequisite knowledge he needs to effectively derive the equations of motion for a double rigid pendulum and then cover more advanced topics for him. In theory, the derivation of the equations of motion for a double rigid pendulum can be done with Newtonian mechanics. Newtonian mechanics is the formulation of classical mechanics that utilizes Newton's laws of motion, most importantly his second law $\textbf{F}=m\mathbf{a}$. Classical mechanics is, as its name suggests, physics based on old-fashioned (or classical) notions of forces, determinism and absolute time. It does not account for relativistic effects nor quantum effects.

Given our second pendulum is in a moving reference frame and Newtonian mechanics does not easily adapt to alternative coordinate systems, like polar coordinates, Lagrangian mechanics makes the equations of motion easier to derive. Lagrangian mechanics also leads nicely into other fields of physics such as Hamiltonian mechanics, the Hamilton-Jacobi equations, quantum mechanics and relativity. 

Figure \figref{dp} shows an animation of the double rigid pendulum. [This page](/doublePendulum/) does interactive simulations of this system if you would like to specify the parameter values yourself. Code Listing \coderef{doublepen} is the Julia code used to generate the above double pendulum animation. I will mention that the Julia code included in this article was largely generated with the help of ChatGPT.

\codelisting{doublepen}{Double pendulum simulation code.}{./code/doublepen}{fidello/code/doublepen.jl}
