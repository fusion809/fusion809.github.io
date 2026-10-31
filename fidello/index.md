@def hassim=false;
@def title="Advanced classical mechanics explanation for Fidello."
@def tag=["Mathematics", "Physics"]
@def mintoclevel=1

~~~
<figure id="dp" style="float: right;">
<img src="/fidello/Double pendulum.png" width="500px"/>
<figcaption>Double pendulum animation frame.</figcaption>
</figure>
~~~
In this article, I explain to Fidello from Foci's Discord server some prerequisite knowledge he needs to effectively derive the equations of motion for a double rigid pendulum and then cover Hamiltonian mechanics. In theory, this derivation can be done with Newtonian mechanics with $\textbf{F}=m\mathbf{a}$, but given our second pendulum is in a moving reference frame and Newtonian mechanics does not easily adapt to alternative coordinate systems, like polar coordinates, Lagrangian mechanics makes the equations of motion easier to derive. Lagrangian mechanics also leads nicely into Hamiltonian mechanics and the path integral formulation of quantum mechanics. Figure \figref{dp} shows a frame from the double pendulum animation [this page](/doublePendulum/) does by default.

\toc
# Partial differentiation
Partial derivatives are just a straightforward extension of derivatives as you understand them already. The key difference is that they are derivatives of multivariable functions instead of single-variable functions. Further, when you differentiate with respect to one variable, all other variables are kept constant. If you would like to build a graphical image of what the partial derivative is, well recall that the derivative of a single-variable function is just the gradient of that function. Consequently, it is also the gradient of the tangent line to that function. The partial derivative is similar to that. As when we have a function with more than a single independent variable, we have the freedom to choose which variable or variables we wish to determine our function's gradient with respect to. 

If you want to visualize this, imagine a surface defined by $z=f(x,y)$ for example. Because a tangent line is two dimensional and we're in three dimensional space, there can be more than one tangent line to a given point. As when constructing the tangent line, we must answer a question: in which direction do we want our line to be tangential to our surface? It could be tangential in the $x$ direction or the $y$ direction, or even in any arbitrary direction in the $x$-$y$ plane. If we want the tangent line in the $x$ direction, we will need the partial derivative with respect to $x$ to get its gradient, whereas if we want the tangent line in the $y$ direction, we will need the partial derivative with respect to $y$ to get its gradient. More generally, for an arbitrary direction in the $x$-$y$ plane we would use the directional derivative as it is called. 

Say $f$ is a function of two variables &mdash; $x$ and $y$. Then its partial derivative with respect to $x$ is formally defined as:

\begin{align*}
    \dfrac{\partial f}{\partial x} &= \lim_{h\rightarrow 0} \dfrac{f(x+h, y)-f(x, y)}{h}.
\end{align*}

This is analogous to the first principles definition of the derivative you would have learned for standard derivatives. Namely, if we are differentiating $f(x)$:
\begin{align*}
    \dfrac{df}{dx} &= \lim_{h\rightarrow 0}\dfrac{f(x+h) - f(x)}{h}.
\end{align*}
The only difference is that now we also have a variable $y$ that our function $f(x,y)$ depends on, but is kept constant when differentiating with respect to $x$.
Say $f(x,y)=2x^2y$. Then $\dfrac{\partial f}{\partial x}$ would be:
\begin{align*}
    \dfrac{\partial f}{\partial x} &= \lim_{h\rightarrow 0} \dfrac{2(x+h)^2y - 2x^2y}{h}.
\end{align*}
Did you notice that $2y$ is a common factor of both terms in our numerator? Let us take it out of our limit altogether, as it is effectively a constant (it doesn't depend on $h$, after all). 

\begin{align*}
\dfrac{\partial f}{\partial x} &= 2y\lim_{h\rightarrow 0} \dfrac{(x+h)^2-x^2}{h}
\end{align*}
Using the binomial theorem, we know that
\begin{align*}
(x+h)^2 &= x^2 + 2xh + h^2.
\end{align*}

Hence:
\begin{align*}
\dfrac{\partial f}{\partial x} &= 2y\lim_{h\rightarrow 0} \dfrac{x^2+2xh + h^2-x^2}{h} \\
&= 2y\lim_{h\rightarrow 0} \dfrac{2xh + h^2}{h} \\
&= 2y\lim_{h\rightarrow 0} 2x + h \\
&= 2y \cdot 2x \\
&= 4xy.
\end{align*}
Okay, so that is the long way of finding out partial derivatives. Fortunately, there is also a faster way. Which is that we apply normal differentiation rules, but any variable we are not differentiating with respect to is treated like it's a constant. So for $f(x,y)=2x^2y$ we would treat the $2y$ as just a constant number and only differentiate the $x^2$. Applying the power rule, we see that $\dfrac{d}{dx}x^2 = 2x$ and hence:
\begin{align*}
    \dfrac{\partial f}{\partial x} &= 2y \cdot 2x \\
    &= 4xy.
\end{align*}

# Heat equation example
One of the simplest partial differential equations is the heat equation. In its simplest form it is:
\begin{align}
    \dfrac{\partial u}{\partial t} &= \alpha \dfrac{\partial^2 u}{\partial x^2} \label{heat}
\end{align}
Where $u(x,t)$ is a function describing the temperature of the one-dimensional object (only one dimensional as the above heat equation only has a single spatial dimension, if we want more spatial dimensions we'll have to add second partial derivatives with respect to those other dimensions, too) we're modelling and $\alpha$ is some constant. I know this equation may look confusing at first glance, but it is actually saying something rather intuitive to you. 

Do you remember from your maths class that the second derivative of a function is its concavity? So if the second derivative is negative, that means our function is concave down &mdash; it is leading to some sort of maximum. Whereas if it is positive, the function is concave up, which means it is going up from some minimum. We have the second derivative with respect to space on the right-hand side of Equation \eqref{heat}, so it essentially returns positive if we have a localized colder area and neighbouring hotter points. Whereas it returns negative if we have a really hot point surrounded by cooler points. The left-hand side of Equation \eqref{heat}, however, represents the rate of change of our function in time.

Equation \eqref{heat} therefore essentially says that if we are in a region where there's a really cold point surrounded by hotter points, the temperature will increase in the short-term for that colder point. Whereas the equation says that the temperature will decrease in a hotter region surrounded by colder points. If you've ever brought something hot out of the oven, or put something of room temperature into the fridge or freezer, this should make intuitive sense to you. 

# Lagrangian mechanics
Lagrangian mechanics is one of three major formulations of classical mechanics with the others being Newtonian mechanics and Hamiltonian mechanics. Classical mechanics as a whole is a branch of physics most applicable to the situations that humans most often encounter, as it does not account for quantum randomness (hence cannot be applied to very small scales), nor does it account for relativistic effects (which are most relevant when dealing with strong gravitational fields or velocities close to the speed of light). In Lagrangian mechanics, we define a quantity known as the Lagrangian of our physical system:
\begin{align*}
    L &= T - V. \label{L}
\end{align*}
Where $T$ is kinetic energy and $V$ is potential energy of our entire physical system. It is important to note that only conservative forces contribute to $V$. Lagrangian mechanics **can** account for non-conservative (i.e. dissipative) forces, too, but we will get to how one handles them later. 

The distinction between conservative and non-conservative forces is whether they change the net energy of your system. For instance, gravity doesn't because any increase in kinetic energy it causes is balanced by a corresponding reduction in potential energy. Likewise, any loss of kinetic energy caused by gravity is balanced by an increase in potential energy. Consequently, it is a conservative force. Air resistance, however, is a non-conservative force because the losses it causes in kinetic energy are not balanced by increases in potential energy. 
## Euler-Lagrange equations
When there are no non-conservative forces, the equations of motion can be derived by first determining our system's Lagrangian and then applying a set of partial differential equations known as the Euler-Lagrange equations to this Lagrangian. The Euler-Lagrange equations are:
\begin{align}
    \hspace{4cm}\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i} - \dfrac{\partial L}{\partial q_i} &= 0, \hspace{2cm} i\in \{1,2,3,...,n\}.\label{ELE}
\end{align}
$q_i$ are the $n$ degrees of freedom of our problem and $\dot{q}_i$ are their time derivatives. The $\dfrac{d}{dt}$ is not a partial derivative and this is purposeful, as when we differentiate $\dfrac{\partial L}{\partial \dot{q}_i}$ with respect to $t$, we treat $q_i$ and $\dot{q}_i$ as functions of $t$ and wherever they appear in $\dfrac{\partial L}{\partial \dot{q}_i}$ we differentiate them against $t$. 

Next, I should explain what "degrees of freedom" are. They are also known as generalized coordinates and they are essentially our variables of interest. For simpler problems, the degrees of freedom in our problem are often the spatial dimensions of our problem. For instance, if we were to study an object in free fall and we were to try to simplify our problem as much as possible and not consider any movement in the horizontal directions then $n$ could be 1 and $q_1$ would just be the height of the object above ground level. 

Likewise, when we study the simple pendulum &mdash; which is a problem in which we have an idealized rigid pendulum whose rod has zero mass, zero friction and for which there is absolutely no air resistance &mdash; typically we only consider a single degree of freedom, which is the angle our pendulum's rod makes with either the negative $y$ axis or the positive $x$ axis, depending on the preferences of whoever is setting up the problem. 
## Derivation of Equations \eqref{ELE}
The derivation of Equations \eqref{ELE} requires knowledge of a field of mathematics known as the calculus of variations as well as a physics principle known as the stationary action principle. If you are interested in learning more, I'd suggest you check out Veritasium's video on this topic ([https://www.youtube.com/watch?v=Q10_srZ-pbs](https://www.youtube.com/watch?v=Q10_srZ-pbs)). In short, I will say that there's a quantity known as the action, which is denoted as $S$, that is defined as:
\begin{align}
    S &= \int_{t_0}^{t_f} L(q_i, \dot{q}_i, t) dt. \label{act}
\end{align}
One interesting property of the action is that nature consistently sets it to its stationary point, and this is regardless of context. Regardless of whether we're examining the motion of planets with Newton's second law, or studying the behaviour of black holes with general relativity, or checking out the behaviour of subatomic particles, there is always an action that is at its stationary point. Traditionally, this stationary point is typically thought of as a minimum (hence why this is also known as the principle of **least** action), but strictly speaking it is not always a minimum. When you use the calculus of variations to find the stationary point of Equation \eqref{act}, Equations \eqref{ELE} pop out. 

## Terminology
In Equation \eqref{ELE}:
\begin{align}
p_i =\dfrac{\partial L}{\partial \dot{q}_i}\label{canmom}
\end{align}
is known as the canonical (or generalized) momentum. This is not that relevant at this stage of your learning, but if you later learn Hamiltonian mechanics, it becomes very relevant. 

## Non-conservative forces
I guess next I should explain how we account for non-conservative forces. Essentially we replace the $0$ in Equations \eqref{ELE} with:

\begin{align}
    Q_i &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \textbf{r}_j}{\partial q_i}. \label{Qir}
\end{align}
Where $\textbf{F}_j$ represents the vector of the $j$th non-conservative force we wish to account for in our equations of motion. As you can see by the summation bounds, we have $m$ non-conservative forces we wish to account for. $\textbf{r}_j$ is the position vector for force $j$. $\textbf{F}_j\cdot\dfrac{\partial \textbf{r}_j}{\partial q_i}$ denotes the dot product between $\textbf{F}_j$ and $\dfrac{\partial \textbf{r}_j}{\partial q_i}$. Alternatively, to calculate $Q_i$ we can use velocities instead of positions:
\begin{align}
    Q_i &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \dot{\textbf{r}}_j}{\partial \dot{q}_i}.\label{Qirdot}
\end{align}
$\dot{\textbf{r}}_j$ is the time derivative of the position vector (that is, the velocity vector) for force $j$. Hence the more generalized form of Equations \eqref{ELE} is:
\begin{align}
    \hspace{4cm}\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i} - \dfrac{\partial L}{\partial q_i} &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \textbf{r}_j}{\partial q_i}, \hspace{2cm} i\in\{1,2,3,...,n\}. \label{ELENCF}
\end{align}
The velocity-based definition of $Q_i$ can also be used, which yields:
\begin{align}
    \hspace{4cm}\dfrac{d}{dt} \dfrac{\partial L}{\partial \dot{q}_i} - \dfrac{\partial L}{\partial q_i} &= \sum_{j=1}^m \textbf{F}_j \cdot \dfrac{\partial \dot{\textbf{r}}_j}{\partial \dot{q}_i}, \hspace{2cm} i\in\{1,2,3,...,n\}. \label{ELENCFrdot}
\end{align}

# Single pendulum
~~~
<figure id="fig1" style="float:right;">
<img src="/fidello/pendulum.svg" width="500px"/>
<figcaption>Single pendulum.</figcaption>
</figure>
~~~
So to punctuate how one applies Equations \eqref{ELE} and \eqref{ELENCF}, I thought I would give you a nice, simple example: the single pendulum. 
## Without air resistance
To apply Equations \eqref{ELE}, we will treat it as a simple pendulum with no friction, a massless rod and no air resistance. The pendulum will consist of a rod with a bob at one end and freely rotating around a central point. We will define this central point as being the origin in Cartesian coordinates. Say the rod is of fixed length $l$ metres. The pendulum bob will be of mass $m$ kilograms. As you can see in Figure \figref{fig1}, we have defined $\theta$ as the angle between the rod and the positive $x$-axis going counter-clockwise. Hence the $x$ and $y$ coordinates of the pendulum bobs are:
\begin{align*}
    x &= l\cos{\theta} \\
    y &= l\sin{\theta}.
\end{align*}
Hence: 
\begin{align*}
    \dot{x} &= -l\dot{\theta}\sin{\theta} \\
    \dot{y} &= l\dot{\theta}\cos{\theta}.
\end{align*}
As for velocity squared, it is defined as:
\begin{align*}
    v^2 &= \dot{x}^2+\dot{y}^2 \\
    &= l^2\dot{\theta}^2\sin^2{\theta} + l^2\dot{\theta}^2\cos^2{\theta} \\
    &= l^2\dot{\theta}^2.
\end{align*}
Hence $v=l\dot{\theta}$. The kinetic energy is:
\begin{align*}
    T &= \dfrac{mv^2}{2} \\
    &= \dfrac{ml^2\dot{\theta}^2}{2}.
\end{align*}
As for potential energy, given this bob is moving solely under the force of gravity:
\begin{align*}
    V &= mgy \\
    &= mgl\sin{\theta}.
\end{align*}
Hence:
\begin{align*}
    L &= T - V \\
    &= \dfrac{ml^2\dot{\theta}^2}{2} - mgl\sin{\theta}.
\end{align*}
$\dfrac{\partial L}{\partial \dot{\theta}}$ is hence:
\begin{align*}
    \dfrac{\partial L}{\partial \dot{\theta}} &= ml^2 \dot{\theta}.
\end{align*}
And its time derivative is:
\begin{align*}
    \dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{\theta}} = ml^2 \ddot{\theta}.
\end{align*}
As for the $\dfrac{\partial L}{\partial \theta}$ term it is:
\begin{align*}
    \dfrac{\partial L}{\partial \theta} &= -mgl\cos{\theta}.
\end{align*}
So Equations \eqref{ELE} becomes:
\begin{align*}
    ml^2 \ddot{\theta} + mgl\cos{\theta} &= 0.
\end{align*}
Dividing by $ml^2$ and taking the second term to the right-hand side yields:
\begin{align}
\ddot{\theta} &= -\dfrac{g}{l}\cos{\theta}.\label{simppen}
\end{align}
Technically speaking it is possible to analytically solve Equation \eqref{simppen}, but it requires the use of non-standard functions called elliptic integrals. If you would like to check out numerical solutions to this equation, see [this page](/simplePendulum/woErrorAnalysis).

## With air resistance
For the sake of simplicity, we will only consider the air resistance experienced by the pendulum bob itself, not by the rod. The air resistance force vector for the pendulum bob would be:
\begin{align*}
    \textbf{F}_D &= -\dfrac{\rho C_D A|v|\textbf{v}}{2}.
\end{align*}
Where $\textbf{v}$ is the velocity vector of the bob. It is equal to $(\dot{x}, \dot{y})$.
Hence:
\begin{align*}
    \textbf{F}_D &= -\dfrac{\rho C_D Al|\dot{\theta}|}{2}(\dot{x}, \dot{y}) \\
    &= -\dfrac{\rho C_D Al|\dot{\theta}|}{2} (-l\dot{\theta}\sin{\theta}, l\dot{\theta}\cos{\theta}) \\
    &= \dfrac{\rho C_D Al^2|\dot{\theta}|\dot{\theta}}{2} (\sin{\theta}, -\cos{\theta}).
\end{align*}
To find the generalized force term, we need to find $\dfrac{\partial \textbf{r}}{\partial \theta}$:
\begin{align*}
    \dfrac{\partial \textbf{r}}{\partial \theta} &= \left(\dfrac{\partial x}{\partial \theta}, \dfrac{\partial y}{\partial \theta}\right) \\
    &= (-l\sin{\theta}, l\cos{\theta}).
\end{align*}
Hence:
\begin{align*}
    Q_{\theta} &= \textbf{F}_D \cdot \dfrac{\partial \textbf{r}}{\partial \theta} \\
    &= \dfrac{\rho C_D Al^2|\dot{\theta}|\dot{\theta}}{2} (\sin{\theta}, -\cos{\theta}) \cdot (-l\sin{\theta}, l\cos{\theta}) \\
    &= \dfrac{\rho C_D Al^2|\dot{\theta}|\dot{\theta}}{2}(-l\sin^2{\theta}-l\cos^2{\theta}) \\
    &= -\dfrac{\rho C_D Al^3|\dot{\theta}|\dot{\theta}}{2}.
\end{align*}
Hence Equations \eqref{ELENCF} become:
\begin{align*}
    ml^2 \ddot{\theta} + mgl\cos{\theta} &= -\dfrac{\rho C_D Al^3|\dot{\theta}|\dot{\theta}}{2}.
\end{align*}
Dividing by $ml^2$ yields:
\begin{align*}
    \ddot{\theta} + \dfrac{g}{l}\cos{\theta} &= -\dfrac{\rho C_D Al|\dot{\theta}|\dot{\theta}}{2m}.
\end{align*}
Subtracting $\dfrac{g}{l}\cos{\theta}$ from both sides yields:
\begin{align*}
    \ddot{\theta} &= -\dfrac{g}{l}\cos{\theta} -\dfrac{\rho C_D Al|\dot{\theta}|\dot{\theta}}{2m}.
\end{align*}
This equation cannot be analytically solved, so numerical solutions are necessary. I numerically solve an even more complex version of the single pendulum equation [here](/singlePendulum/).

If you would like to see my derivation of the double rigid pendulum equations, it is [here](/doublePendulum/derivation/).

# Hamiltonian mechanics
In Hamiltonian mechanics, we define the Hamiltonian as:
\begin{align}
H &= T + V.\label{H}
\end{align}
Which is the total energy of our system. Another difference from Lagrangian mechanics is that we define $H$ in terms of $p_i$ and $q_i$ instead of $q_i$ and $\dot{q}_i$. The equations of motion with no non-conservative forces are:
\begin{align}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} \label{Heq1}\\
\dot{q}_i &= \dfrac{\partial H}{\partial p_i}.\label{Heq2}
\end{align}
Equations \eqref{Heq1} and \eqref{Heq2} are often called Hamilton's equations. To account for non-conservative forces, we merely add the generalized force terms for them to the right-hand side of Equation \eqref{Heq1}. In other words:
\begin{align}
\dot{p}_i &= -\dfrac{\partial H}{\partial q_i} + Q_i. \label{HNConseq1}
\end{align}
This creates a new kind of space that we can visualize our solution on that has dimensions $q_i$ and $p_i$. This type of space is known as phase space. 
## Example: simple pendulum
\begin{align}
H &= \dfrac{ml^2 \dot{\theta}^2}{2} + mgl\sin{\theta}. \label{Hsimp}
\end{align}
Problem is, however, that this is written in terms of $q_i$ and $\dot{q}_i$ not $q_i$ and $p_i$. 

\begin{align*}
p_{\theta} &= \dfrac{\partial L}{\partial \dot{\theta}}
\end{align*}
and so:
\begin{align}
p_{\theta} &= ml^2\dot{\theta}. \label{pthetasimp}
\end{align}
Rewriting Equation \eqref{pthetasimp} with $\dot{\theta}$ as its subject:
\begin{align}
\dot{\theta} &= \dfrac{p_{\theta}}{ml^2}. \label{thetadotsimp}
\end{align}
Substituting Equation \eqref{thetadotsimp} into Equation \eqref{Hsimp} yields:
\begin{align*}
H &= \dfrac{ml^2}{2} \left( \dfrac{p_{\theta}}{ml^2} \right)^2 + mgl\sin{\theta} \\
&= \dfrac{ml^2}{2} \dfrac{p_{\theta}^2}{m^2l^4} + mgl\sin{\theta}.
\end{align*}
So our final Hamiltonian is:
\begin{align}
H &= \dfrac{p_{\theta}^2}{2ml^2} + mgl\sin{\theta}. \label{Hsimpf}
\end{align}
Hence the equations of motion are:
\begin{align}
\dot{p}_{\theta} &= -\dfrac{\partial H}{\partial \theta} \nonumber \\
 &= - mgl\cos{\theta}\label{pdotthetasimp}
\end{align}
and
\begin{align}
\dot{\theta} &= \dfrac{\partial H}{\partial p_{\theta}} \nonumber\\
&= \dfrac{p_{\theta}}{ml^2}.\label{thetadotsimp}
\end{align}
As we can see, Hamiltonian mechanics yields a series of first-order ordinary differential equations (ODEs) instead of second-order ODEs. Which is important because Lagrangian mechanics and Newtonian mechanics both yield second-order ODE systems. That being said, if you rewrite Equation \eqref{thetadotsimp} to make $p_{\theta}$ its subject and substitute it into Equation \eqref{pdotthetasimp} and rearrange it to make $\ddot{\theta}$ the subject, you will obtain Equation \eqref{simppen}.

## Poisson brackets
Poisson brackets provide us an elegant and powerful way to rewrite Hamilton's equations. If we have $n$ degrees of freedom and two functions $f$ and $g$ that depend on $p_i$ and $q_i$:[^1] 

\begin{align}
\{f, g\} &= \sum_{i=1}^n \left( \dfrac{\partial f}{\partial q_i} \dfrac{\partial g}{\partial p_i} - \dfrac{\partial f}{\partial p_i}\dfrac{\partial g}{\partial q_i} \right). \label{Poissonfg}
\end{align}

Consequently:

\begin{align}
\dot{q}_k &= \{q_k, H\}\label{qdotk}\\
\dot{p}_k &= -\{p_k, H\}. \label{pdotk}
\end{align}

We can establish Equations \eqref{qdotk} with
\begin{align}
\{q_k, H\} &= \sum_{i=1}^n \left(\dfrac{\partial q_k}{\partial q_i}\dfrac{\partial H}{\partial p_i} - \dfrac{\partial q_k}{\partial p_i}\dfrac{\partial H}{\partial q_i}\right).\label{qkHPoisson}
\end{align}
$\dfrac{\partial q_k}{\partial q_i}=\delta_{ik}$, where $\delta_{ik}$ is the Kronecker delta function which is is equal to one if $i=k$ and zero otherwise. This equation holds because we assume our degrees of freedom or coordinates are independent of each other. For instance, $y$ is independent of $x$ in Cartesian coordinates and $\theta$ is independent of $r$ in polar coordinates. As for $\dfrac{\partial q_k}{\partial p_i}$ it is zero. Hence Equation \eqref{qkHPoisson} becomes:
\begin{align}
\{q_k, H\} &= \sum_{i=1}^n \delta_{ik} \dfrac{\partial H}{\partial p_i}.
\end{align}
Summing over $i$ yields only one term, thanks to the Kronecker delta, namely the $i=k$ term, as all other terms have a zero multiplier:
\begin{align}
    \{q_k, H\} &= \dfrac{\partial H}{\partial p_k}.\label{QkHPoisson}
\end{align}
Thanks to Hamilton's equations, we know the right-hand side of Equation \eqref{QkHPoisson} is equal to $\dot{q}_k$.

We can establish Equations \eqref{pdotk} with:
\begin{align}
\{p_k, H\} &= \sum_{i=1}^n \left(\dfrac{\partial p_k}{\partial q_i}\dfrac{\partial H}{\partial p_i} - \dfrac{\partial p_k}{\partial p_i}\dfrac{\partial H}{\partial q_i}\right)
\end{align}
$\dfrac{\partial p_k}{\partial q_i}$ is zero and $\dfrac{\partial p_k}{\partial p_i}$ is $\delta_{ik}$. Hence:
\begin{align}
\{p_k, H\} &= \sum_{i=1}^n \left(-\delta_{ik} \dfrac{\partial H}{\partial q_i}\right) \\
&= -\dfrac{\partial H}{\partial q_k}. \label{pkHR}
\end{align}
According to Hamilton's equations, the right-hand side of Equations \eqref{pkHR} should equal $\dot{p}_k$.

Poisson brackets also have several interesting properties, such as antisymmetry:
\begin{align}
\{A, B\} &= -\{B, A\}. \label{asym}
\end{align}
Further:
\begin{align}
\{q_k, p_k\} &= \sum_{i=1}^n \left( \dfrac{\partial q_k}{\partial q_i} \dfrac{\partial p_k}{\partial p_i} - \dfrac{\partial q_k}{\partial p_i} \dfrac{\partial p_k}{\partial q_i} \right) \nonumber\\
&= \sum_{i=1}^n \delta_{ik}\delta_{ik} \nonumber\\
&= 1 \label{qkpk} \\
\{q_j, q_k\} &= \sum_{i=1}^n \left( \dfrac{\partial q_j}{\partial q_i} \dfrac{\partial q_k}{\partial p_i} - \dfrac{\partial q_j}{\partial p_i} \dfrac{\partial q_k}{\partial q_i} \right) \nonumber\\
&= 0 \label{qjqk}\\
\{p_j, p_k\} &= \sum_{i=1}^n \left(\dfrac{\partial p_j}{\partial q_i}\dfrac{\partial p_k}{\partial p_i} - \dfrac{\partial p_j}{\partial p_i}\dfrac{\partial p_k}{\partial q_i}\right) \nonumber\\
&= 0. \label{pjpk}
\end{align}
The quality of two functions, say $f$ and $g$, having a Poisson bracket that equals zero is phrased as $f$ and $g$  **commuting** under the Poisson bracket. 

If a function $f(q_i, p_i)$ is not explicitly dependent on $t$ and only depends on $t$ indirectly through $q_i$ and $p_i$ then, via the chain rule:
\begin{align}
\dfrac{df}{dt} &= \sum_{i=1}^n \dfrac{\partial f}{\partial q_i}\dfrac{d q_i}{dt} + \dfrac{\partial f}{\partial p_i}\dfrac{dp_i}{dt}.
\end{align}
We can use Hamilton's equations to replace $\dfrac{dq_i}{dt}$ and $\dfrac{dp_i}{dt}$:
\begin{align}
\dfrac{df}{dt} &= \sum_{i=1}^n \dfrac{\partial f}{\partial q_i}\dfrac{\partial H}{\partial p_i} - \dfrac{\partial f}{\partial p_i}\dfrac{\partial H}{\partial q_i}.
\end{align}
Interestingly, the right-hand side is equal to the Poisson bracket of $f$ with $H$, in other words:
\begin{align}
\dfrac{df}{dt} &= \{f, H\}.
\end{align}
Hence commuting with the Hamiltonian under the Poisson bracket is a sufficient condition for a quantity that does not explicitly depend on time to have a zero time derivative and hence be conserved. 

# Noether's theorem
Noether's theorem essentially states that for every continuous symmetry in the Lagrangian, there is a corresponding conservation law. For instance, time translation symmetry of the Lagrangian yields the conservation of energy. Say we increase time according to the transformation: $t\rightarrow t+\epsilon$. 

\begin{align}
L(q_i(t+\epsilon), \dot{q}_i(t+\epsilon), t+\epsilon) &\approx L(q_i, p_i, t) + \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \dfrac{dq_i}{dt} \epsilon + \dfrac{\partial L}{\partial \dot{q}_i} \dfrac{d\dot{q}_i}{dt} \epsilon \nonumber\\
L(q_i(t+\epsilon), \dot{q}_i(t+\epsilon), t+\epsilon) - L(q_i, \dot{q}_i, t) &\approx \sum_{i=1}^n\dfrac{\partial L}{\partial q_i} \dfrac{dq_i}{dt} \epsilon + \dfrac{\partial L}{\partial \dot{q}_i} \dfrac{d\dot{q}_i}{dt} \epsilon. \label{diffL}
\end{align}
We will call the left-hand side of Equation \eqref{diffL} $\Delta L$, hence Equation \eqref{diffL} becomes:
\begin{align}
\Delta L &\approx \sum_{i=1}^n \dfrac{\partial L}{\partial q_i} \dfrac{dq_i}{dt} \epsilon + \dfrac{\partial L}{\partial \dot{q}_i} \dfrac{d\dot{q}_i}{dt} \epsilon. \label{DeltaL1}
\end{align}
The Euler-Lagrange equation yields that:
\begin{align}
\dfrac{\partial L}{\partial q_i} = \dfrac{d}{dt}\left(\dfrac{\partial L}{\partial \dot{q}_i}\right). \label{ELEnh}
\end{align}
If you are wondering what about if we have non-conservative forces? Wouldn't that mean Equations \eqref{ELEnh} would not hold as there would be a generalized force term? Yes and in that event the Lagrangian would indeed change in time and conservation of energy would not hold, as our system would be losing energy to those non-conservative forces. 

Substituting Equations \eqref{ELEnh} into Equation \eqref{DeltaL1} yields:
\begin{align}
\Delta L &\approx \sum_{i=1}^n \epsilon\left(\dfrac{d}{dt}\dfrac{\partial L}{\partial \dot{q}_i}\dot{q}_i + \dfrac{\partial L}{\partial \dot{q}_i}\ddot{q}_i\right). \label{DeltaL2}
\end{align}
Substituting Equations \eqref{canmom} into Equation \eqref{DeltaL2} yields:
\begin{align}
\Delta L &\approx \epsilon\sum_{i=1}^n \left(\dfrac{dp_i}{dt} \dot{q}_i + p_i \ddot{q}_i\right)\\
&= \epsilon\sum_{i=1}^n \left(\dot{p}_i \dot{q}_i + p_i \ddot{q}_i\right). \label{DeltaL3}
\end{align}
We can rewrite the bracketed expression on the right-hand side of Equation \eqref{DeltaL3} using the product rule of differentiation as:
\begin{align}
\Delta L &\approx \epsilon\sum_{i=1}^n \dfrac{d}{dt}\left(p_i \dot{q}_i\right) \\
&= \epsilon\dfrac{d}{dt}\sum_{i=1}^n p_i \dot{q}_i.
\end{align}
Hence if $L$ does not change in time, we get:
\begin{align}
\epsilon \dfrac{d}{dt} \sum_{i=1}^n p_i \dot{q}_i &= 0.
\end{align}
$\epsilon$ should not matter, so we will divide both sides by it:
\begin{align}
\dfrac{d}{dt} \sum_{i=1}^n p_i \dot{q}_i &= 0.
\end{align}
$L$ does not change, so $\dfrac{dL}{dt}=0$, too. Hence:
\begin{align}
\dfrac{d}{dt} \left(\sum_{i=1}^n p_i \dot{q}_i - L\right) &= 0.\label{constE}
\end{align}
Interestingly, $\displaystyle\sum_{i=1}^n p_i \dot{q}_i=2T$, that is it is twice the kinetic energy of our system. Hence if we substitute Equation \eqref{L} into Equation \eqref{constE}, we obtain:
\begin{align}
\dfrac{d}{dt} \left(2T-(T-V)\right) &= 0 \nonumber\\
\dfrac{d}{dt} \left(T+V\right) &= 0.\label{constE2}
\end{align}
Hence substituting Equation \eqref{H} into Equation \eqref{constE2} yields:
\begin{align}
\dfrac{dH}{dt} &= 0.\label{constE3}
\end{align}
In other words, we obtain that the Hamiltonian &mdash; the total energy of the system &mdash; does not change, hence it is conserved.

# Reference list
[^1]: Cline D. (2026). [15.2: Poisson bracket Representation of Hamiltonian Mechanics](https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/15%3A_Advanced_Hamiltonian_Mechanics/15.02%3A_Poisson_bracket_Representation_of_Hamiltonian_Mechanics). In *Variational Principles in Classical Mechanics*. LibreTexts.
