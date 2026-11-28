### Integrating factor method
Say we have a differential equation of the form:
\begin{align}
y'(x) + P(x)y(x) &= Q(x). \label{IntDE}
\end{align}
Then we can use the integrating factor method, where:
\begin{align*}
I(x) &= e^{\int P(x) dx}.
\end{align*}
Say we multiply Equation \eqref{IntDE} by $I(x)$:
\begin{align*}
I(x)y'(x) + I(x)P(x)y(x) &= I(x)Q(x). \label{IntIDE}
\end{align*}
Keeping in mind that $I'(x)=P(x)I(x)$. Then, using product rule of differentiation:
\begin{align*}
\dfrac{d}{dx}(I(x)y(x)) &= I(x)y'(x) + I'(x)y(x) \\
&= I(x) y'(x) + I(x)P(x) y(x).
\end{align*}
Hence:
\begin{align*}
\dfrac{d}{dx}(I(x)y(x)) &= I(x)Q(x).
\end{align*}
Integrating both sides yields:
\begin{align*}
I(x)y(x) &= \int I(x) Q(x) dx \\
y(x) &= \dfrac{1}{I(x)} \int I(x) Q(x) dx.
\end{align*}
We can also combine this technique with the rearranging and integrating technique previously covered. For example, say we have:
\begin{align*}
y'' + y' = \cos{x}.
\end{align*}
We could integrate both sides:
\begin{align*}
y' + y &= \sin{x} + C_1
\end{align*}
Hence our integrating factor must be:
\begin{align*}
I(x) &= e^{\int 1 dx} \\
&= e^{x+C_0} \\
&= C_0 e^{x}.
\end{align*}
The $C_0$ between the last two lines are not the same; as they are arbitrary constants we do not need to keep them strictly the same between these lines. 
\begin{align*}
\dfrac{d}{dx} (C_0 e^x y(x)) &= C_0 e^x \sin{x} + C_1 C_0 e^{x} \\
C_0 e^xy(x) &= C_0 \int (e^x \sin{x} + C_1 e^{x}) dx.
\end{align*}
Dividing both sides by $C_0$ we get:
\begin{align*}
e^xy(x) &= \int (e^x \sin{x} + C_1 e^{x}) dx.
\end{align*}
Via Euler's formula, we obtain $\sin{x} = \dfrac{e^{ix}-e^{-ix}}{2i}$. Hence:
\begin{align*}
\int e^x \sin{x} dx &= \int e^x \dfrac{e^{ix}-e^{-ix}}{2i} dx \\
&= -\dfrac{i}{2} \int e^{(i+1)x}-e^{(1-i)x} dx \\
&= -\dfrac{i}{2} \left[ \dfrac{e^{(i+1)x}}{i+1}-\dfrac{e^{(1-i)x}}{1-i} \right] + C_2 \\
&= -\dfrac{i}{2} \left[ \dfrac{e^{(i+1)x}(1-i)}{2}-\dfrac{e^{(1-i)x}(1+i)}{2} \right] + C_2 \\
&= \dfrac{e^x}{4i} \left[(e^{ix} - e^{-ix}) -i(e^{ix}+e^{ix})\right] + C_2 \\
&= \dfrac{e^x \sin{x}}{2} - \dfrac{e^x \cos{x}}{2} + C_2. 
\end{align*}
Alternatively, we could have tackled this integral via integration by parts. Hence our differential equation becomes:
\begin{align*}
e^x y(x)  &= \dfrac{e^x \sin{x}}{2} - \dfrac{e^x \cos{x}}{2} + C_2 + C_1e^{x} \\
y(x) &= \dfrac{\sin{x}}{2} - \dfrac{\cos{x}}{2} + C_1 + C_2 e^{-x}.
\end{align*}