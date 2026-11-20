~~~
<figure id="dp" style="float: right;">
<img src="/Double_pendulum.gif" width="500px"/>
<figcaption>Double pendulum animation.</figcaption>
</figure>
~~~
In this article, I explain to Fidello from Foci's Discord server some prerequisite knowledge he needs to effectively derive the equations of motion for a double rigid pendulum and then cover more advanced topics for him. In theory, the derivation of the equations of motion for a double rigid pendulum can be done with Newtonian mechanics with $\textbf{F}=m\mathbf{a}$, but given our second pendulum is in a moving reference frame and Newtonian mechanics does not easily adapt to alternative coordinate systems, like polar coordinates, Lagrangian mechanics makes the equations of motion easier to derive. If you would like to see my derivation of the double rigid pendulum equations, it is [here](/doublePendulum/derivation/). Lagrangian mechanics also leads nicely into Hamiltonian mechanics and the path integral formulation of quantum mechanics. 

Figure \figref{dp} shows an animation of the double rigid pendulum. [This page](/doublePendulum/) does interactive simulations of this system if you would like to specify the parameter values yourself. \coderef{doublepen} is the Julia code used to generate the above double pendulum animation. I will mention that the Julia code included in this article was largely generated with the help of ChatGPT.

\codelisting{doublepen}{Double pendulum simulation code.}{./code/doublepen}{fidello/code/doublepen.jl}
