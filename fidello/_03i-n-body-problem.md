### Gravitational $n$-body problem
In this section, we will try to generalize our Lagrangian and system of equations from the previous section so that they are applicable to a gravitational system with an arbitrary number of bodies. Let $M$ denote the static central mass, $m_i$ the $i$th moving mass, $r_i$ the distance from the central mass to the $i$th moving mass, $\theta_i$ is the angle between the line formed by $r_i$ and the positive $x$-axis, measured counter-clockwise. Which is essentially how we defined things in Figure \figref{threebody}.

So kinetic energy for this system should be pretty straight forward, namely it should simply be the sum of the kinetic energy of each individual body, in other words it should be:
\begin{align*}
T &= \sum_{i=1}^{n-1} \dfrac{m_i(\dot{r}_i^2+r_i^2\dot{\theta}_i^2)}{2}.
\end{align*}
Where $m_i$ is the $i$th moving body. As for potential energy, well there will be terms for the interaction between each of the $n$ bodies. For interactions between the central mass and the moving bodies, these will simply be of the form $-\dfrac{GMm_i}{r_i}$. As for the interactions between the moving bodies, their potential energy terms will be of the form $-\dfrac{Gm_i m_j}{r_{ij}}$ where $r_{ij}$ is the distance between masses $i$ and $j$. Hence:
\begin{align*}
x_i &= r_i \cos{\theta_i},\hspace{0.2cm}x_j &&= r_j \cos{\theta_j} \\
y_i &= r_i \sin{\theta_i},\hspace{0.2cm}y_j &&= r_j \sin{\theta_j}.
\end{align*}
Let us define the change in $x$ between masses $i$ and $j$ as $\Delta x_{ij}$ and the corresponding change in $y$ is $\Delta y_{ij}$. 
\begin{align*}
\Delta x_{ij} &= r_j \cos{\theta_j} - r_i \cos{\theta_i} \\
\Delta y_{ij} &= r_j \sin{\theta_j} - r_i \sin{\theta_i}.
\end{align*}
Hence their squares would be:
\begin{align*}
\Delta x_{ij}^2 &= r_j^2 \cos^2{\theta_j} + r_i^2\cos^2{\theta_i} - 2r_ir_j \cos{\theta_i}\cos{\theta_j} \\
\Delta y_{ij}^2 &= r_j^2 \sin^2{\theta_j} + r_i^2\sin^2{\theta_i} - 2r_ir_j \sin{\theta_i}\sin{\theta_j}.
\end{align*}
Hence the squared distance between $m_i$ and $m_j$ is, by the Pythagorean theorem:
\begin{align*}
r_{ij}^2 &= \Delta x_{ij}^2 + \Delta y_{ij}^2 \\
&= r_j^2 \cos^2{\theta_j} + r_i^2\cos^2{\theta_i} - 2r_ir_j \cos{\theta_i}\cos{\theta_j} + r_j^2 \sin^2{\theta_j} + r_i^2\sin^2{\theta_i} - 2r_ir_j \sin{\theta_i}\sin{\theta_j} \\
&= r_i^2 (\cos^2{\theta_i}+\sin^2{\theta_i}) + r_j^2 (\cos^2{\theta_j}+\sin^2{\theta_j}) - 2r_i r_j (\cos{\theta_i}\cos{\theta_j} + \sin{\theta_i}\sin{\theta_j}) \\
&= r_i^2 + r_j^2 - 2r_i r_j \cos{(\theta_j-\theta_i)}.
\end{align*}
Square rooting yields:
\begin{align*}
r_{ij} &= \sqrt{r_i^2 + r_j^2 - 2r_i r_j \cos{(\theta_j-\theta_i)}}.
\end{align*}
As we can see, $r_{ij} = r_{ji}$ due to cosine being an even function. Hence the potential energy is:
\begin{align*}
V &= -\sum_{i=1}^{n-1} \left(\dfrac{GMm_i}{r_i} + \sum_{j>i}^{n-1} \dfrac{Gm_i m_j}{r_{ij}}\right).
\end{align*}
Why the $j>i$ on the second sum? Well, to ensure we do not count the same terms but with a different order of their indices (e.g. $\dfrac{Gm_1 m_2}{r_{12}}$ and $\dfrac{Gm_2 m_1}{r_{21}}$ are equivalent despite the different order of indices) twice. 
Hence:
\begin{align*}
L &= \sum_{i=1}^{n-1} \left(\dfrac{m_i (\dot{r}_i^2 + r_i^2 \dot{\theta}_i^2)}{2} + \dfrac{GMm_i}{r_i} + \sum_{j>i}^{n-1} \dfrac{Gm_i m_j}{r_{ij}}\right)
\end{align*}
Next we will calculate $\dfrac{\partial r_{ij}}{\partial q_k}$:
\begin{align*}
\dfrac{\partial r_{ij}}{\partial r_k} &= \dfrac{2r_i \delta_{ik} + 2r_j \delta_{jk} - 2(\delta_{ik}r_j + r_i \delta_{jk})\cos{(\theta_j-\theta_i)}}{2r_{ij}} \\
&= \dfrac{r_i \delta_{ik} + r_j \delta_{jk} - (\delta_{ik}r_j + r_i \delta_{jk})\cos{(\theta_j-\theta_i)}}{r_{ij}} \\
\dfrac{\partial r_{ij}}{\partial \theta_k} &= \dfrac{2r_i r_j\sin{(\theta_j-\theta_i)}(\delta_{jk}-\delta_{ik})}{2r_{ij}} \\
&= \dfrac{r_i r_j\sin{(\theta_j-\theta_i)}(\delta_{jk}-\delta_{ik})}{r_{ij}}.
\end{align*}
Hence the canonical momentum for $r_k$ is:
\begin{align*}
p_{r_k} &= \dfrac{\partial L}{\partial \dot{r}_k} \\
&= \sum_{i=1}^{n-1} \left(\dfrac{m_i (2\dot{r}_i \delta_{ik})}{2}\right) \\
&= m_k \dot{r}_k.
\end{align*}
Hence its time derivative is:
\begin{align*}
\dot{p}_{r_k} &= m_k \ddot{r}_k.
\end{align*}
As for the canonical momentum for $\theta_k$, it is:
\begin{align*}
p_{\theta_k} &= \dfrac{\partial L}{\partial \dot{\theta}_k} \\
&= \sum_{i=1}^{n-1} \left(\dfrac{m_i (2r_i^2 \dot{\theta}_i \delta_{ik})}{2}\right) \\
&= m_k r_k^2 \dot{\theta}_k.
\end{align*}
Differentiating with respect to $t$ yields:
\begin{align*}
\dot{p}_{\theta_k} &= 2m_k r_k \dot{r_k} \dot{\theta}_k + m r_k^2 \ddot{\theta}_k.
\end{align*}
As for $\dfrac{\partial L}{\partial q_k}$:
\begin{align*}
\dfrac{\partial L}{\partial r_k} &= \sum_{i=1}^{n-1} \left(\dfrac{m_i (2r_i\delta_{ik} \dot{\theta}_i^2)}{2} - \dfrac{GMm_i\delta_{ik}}{r_i^2} - \sum_{j>i}^{n-1} \dfrac{Gm_i m_j\left(r_i \delta_{ik} + r_j \delta_{jk} - (\delta_{ik}r_j + r_i \delta_{jk})\cos{(\theta_j-\theta_i)}\right)}{r_{ij}^3}\right)\\
&= m_k r_k \dot{\theta}_k^2 - \dfrac{GMm_k}{r_k^2} - \sum_{i=1}^{n-1}\sum_{j>i}^{n-1} \dfrac{Gm_i m_j\left(r_i \delta_{ik} + r_j \delta_{jk} - (\delta_{ik}r_j + r_i \delta_{jk})\cos{(\theta_j-\theta_i)}\right)}{r_{ij}^3}.
\end{align*}
Let us work on the final double sum. Factoring out the Kronecker deltas with like indices:
\begin{align*}
\sum_{i=1}^{n-1}\sum_{j>i}^{n-1} \dfrac{Gm_i m_j\left(r_i \delta_{ik} + r_j \delta_{jk} - (\delta_{ik}r_j + r_i \delta_{jk})\cos{(\theta_j-\theta_i)}\right)}{r_{ij}^3} &= \sum_{i=1}^{n-1}\sum_{j>i}^{n-1} \dfrac{Gm_i m_j\left[\delta_{ik}(r_i - r_j \cos{(\theta_j-\theta_i)} )+ \delta_{jk}(r_j -r_i \cos{(\theta_j-\theta_i)})\right]}{r_{ij}^3}.
\end{align*}
The $\delta_{ik}$ imposes that $i=k$. Since $j>i$ we get $j=k+1, k+2, k+3, ...,  n-1$. So we are left only with the $j$ sum for terms multiplied by $\delta_{ik}$. As for $\delta_{jk}$, it imposes that $j=k$ which then requires, given $j>i$, that $i=1,2,3,...,k-1$. 
\begin{align*}
\sum_{i=1}^{n-1}\sum_{j>i}^{n-1} \dfrac{Gm_i m_j\left[\delta_{ik}(r_i - r_j \cos{(\theta_j-\theta_i)} )+ \delta_{jk}(r_j -r_i \cos{(\theta_j-\theta_i)})\right]}{r_{ij}^3} &= \sum_{j=k+1}^{n-1} \dfrac{Gm_{k}m_j (r_k-r_j\cos{(\theta_j-\theta_k)})}{r_{kj}^3} + \sum_{i=1}^{k-1} \dfrac{Gm_{i}m_k (r_k-r_i\cos{(\theta_k-\theta_i)})}{r_{ik}^3}.
\end{align*}
We can combine these sums into one, as $r_{ij}=r_{ji}$ and $\cos{(\theta_j-\theta_i)}=\cos{(\theta_i-\theta_j)}$:
\begin{align*}
\sum_{j=k+1}^{n-1} \dfrac{Gm_{k}m_j (r_k-r_j\cos{(\theta_j-\theta_k)})}{r_{kj}^3} + \sum_{i=1}^{k-1} \dfrac{Gm_{i}m_k (r_k-r_i\cos{(\theta_k-\theta_i)})}{r_{ik}^3} &= \sum_{\substack{i=1 \\ i\neq k}}^{n-1} \dfrac{Gm_{i}m_k (r_k-r_i\cos{(\theta_k-\theta_i)})}{r_{ik}^3}.
\end{align*}
Hence $\dfrac{\partial L}{\partial r_k}$ becomes:
\begin{align*}
\dfrac{\partial L}{\partial r_k} &= m_k r_k \dot{\theta}_k^2 - \dfrac{GMm_k}{r_k^2} - \sum_{\substack{i=1 \\ i\neq k}}^{n-1} \dfrac{Gm_{i}m_k (r_k-r_i\cos{(\theta_k-\theta_i)})}{r_{ik}^3}.
\end{align*}
Hence the $r_k$ equation of motion is:
\begin{align*}
m_k\ddot{r}_k &= m_k r_k \dot{\theta}_k^2 - \dfrac{GMm_k}{r_k^2} - \sum_{\substack{i=1 \\ i\neq k}}^{n-1} \dfrac{Gm_{i}m_k (r_k-r_i\cos{(\theta_k-\theta_i)})}{r_{ik}^3}.
\end{align*}
Dividing by $m_k$ yields:
\begin{align*}
\ddot{r}_k &= r_k \dot{\theta}_k^2 - \dfrac{GM}{r_k^2} - \sum_{\substack{i=1 \\ i\neq k}}^{n-1} \dfrac{Gm_{i} (r_k-r_i\cos{(\theta_k-\theta_i)})}{r_{ik}^3}.
\end{align*}
As for $\theta_k$, the canonical momentum and its time derivative are:
\begin{align*}
p_{\theta_k} &= \dfrac{\partial L}{\partial \dot{\theta}_k} \\
&= \sum_{i=1}^{n-1} \left(\dfrac{2m_i r_i^2 \dot{\theta}_i \delta_{ik}}{2}\right) \\
&= m_k r_k^2 \dot{\theta}_k.
\end{align*}
Hence the time derivative of $p_{\theta_k}$ is:
\begin{align*}
\dot{p}_{\theta_k} &= 2m_kr_k \dot{r}_k \dot{\theta}_k + mr_k^2 \ddot{\theta}_k.
\end{align*}
As for the other component of Equations \eqref{ELE}:
\begin{align*}
\dfrac{\partial L}{\partial \theta_k} &= -\sum_{i=1}^{n-1} \sum_{j>i}^{n-1} \dfrac{Gm_im_j}{r_{ij}^2} \dfrac{\partial r_{ij}}{\partial \theta_k} \\
&= -\sum_{i=1}^{n-1} \sum_{j>i}^{n-1} \dfrac{Gm_im_jr_i r_j\sin{(\theta_j-\theta_i)}(\delta_{jk}-\delta_{ik})}{r_{ij}^3} \\
&= \sum_{j=k+1}^{n-1}\dfrac{Gm_k m_jr_k r_j\sin{(\theta_j-\theta_k)}}{r_{kj}^3} - \sum_{i=1}^{k-1} \dfrac{Gm_i m_kr_i r_k\sin{(\theta_k-\theta_i)}}{r_{ik}^3} \\
&= \sum_{\substack{i=1\\ i\neq k}}^{n-1} \dfrac{Gm_i m_kr_i r_k\sin{(\theta_i-\theta_k)}}{r_{ik}^3}.
\end{align*}
Hence the Euler-Lagrange equations for $\theta_k$ can be expressed as:
\begin{align*}
2m_kr_k \dot{r}_k \dot{\theta}_k + mr_k^2 \ddot{\theta}_k &= \sum_{\substack{i=1\\i\neq k}}^{n-1} \dfrac{Gm_i m_kr_i r_k\sin{(\theta_i-\theta_k)}}{r_{ik}^3}.
\end{align*}
Dividing both sides by $mr_k^2$ yields:
\begin{align*}
\dfrac{2\dot{r}_k\dot{\theta}_k}{r_k} + \ddot{\theta}_k &= \sum_{\substack{i=1\\i\neq k}}^{n-1} \dfrac{Gm_i r_i \sin{(\theta_i-\theta_k)}}{r_k r_{ik}^3}.
\end{align*}
Moving the first term on the left-hand side to the right-hand side yields:
\begin{align*}
\ddot{\theta}_k &= -\dfrac{2\dot{r}_k\dot{\theta}_k}{r_k} + \sum_{\substack{i=1\\i\neq k}}^{n-1} \dfrac{Gm_i r_i \sin{(\theta_i-\theta_k)}}{r_k r_{ik}^3}.
\end{align*}

An important verification of our working is that if $n=3$, do we get the same equations as we obtained in [gravitational three-body problem](#gravitational_three-body_problem). If $n=3$, we get for $\ddot{r}_1$:

\begin{align*}
\ddot{r}_1 &= r_1\dot{\theta}_1^2 - \dfrac{GM}{r_1^2} - \dfrac{Gm_2(r_1 - r_2\cos{(\theta_1-\theta_2)})}{r_{12}^3} \\
&= r_1\dot{\theta}_1^2 - \dfrac{GM}{r_1^2} - \dfrac{Gm_2(r_1 - r_2\cos{(\theta_2-\theta_1)})}{r_{12}^3}.
\end{align*}
This matches our equation for $\ddot{r}_1$ for the three-body problem. As for $\ddot{r}_2$:
\begin{align*}
\ddot{r}_2 &= r_2\dot{\theta}_2^2 - \dfrac{GM}{r_2^2} - \dfrac{Gm_2(r_2 - r_1\cos{(\theta_2-\theta_1)})}{r_{12}^3}.
\end{align*}
Which also matches our three-body system equation for $\ddot{r}_2$. As for $\ddot{\theta}_1$:
\begin{align*}
\ddot{\theta}_1 &= -\dfrac{2\dot{r}_1\dot{\theta}_1}{r_1} + \dfrac{Gm_2r_2 \sin{(\theta_2-\theta_1)}}{r_1 r_{12}^3}.
\end{align*}
Which also matches our three-body system equation for $\ddot{\theta}_1$. As for $\ddot{\theta}_2$:
\begin{align*}
\ddot{\theta}_2 &= -\dfrac{2\dot{r}_2\dot{\theta}_2}{r_2} + \dfrac{Gm_1r_1 \sin{(\theta_1-\theta_2)}}{r_2 r_{12}^3} \\
&= -\dfrac{2\dot{r}_2\dot{\theta}_2}{r_2} - \dfrac{Gm_1r_1 \sin{(\theta_2-\theta_1)}}{r_2 r_{12}^3}.
\end{align*}
Which also matches our three-body system equation for $\ddot{\theta}_2$.

\codelisting{nbody}{N-body problem simulation code.}
```julia:./code/nbody
\juliacode{fidello/code/nbody.jl}
```

\coderef{nbody} produces the animation in Figure \figref{nbodyan}. As you can see, the orbits majorly and seemingly unpredictably change over time. This is pretty typical in the $n$-body problem: either the orbits are fairly stable and unchanging (like the orbits of the planets in the Solar System), or they are unstable and prone to objects being ejected from the system or colliding with one another. 

~~~
<figure id="nbodyan" style="float: right; width: 100%">
<img src="/n_body_system.gif" width="800px"/>
<figcaption>A four-body system simulation.</figcaption>
</figure>
~~~
