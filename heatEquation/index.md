@def title = "Numerically solving the 1D heat equation"
@def hassim = true
@def pde = true
@def ids = ["tableOutputs", "3DPlot"];
@def params = (alpha=(val=1, desc="Heat propogation parameter."), x0=(val=0, desc="Start of spatial domain."), x1=(val=0, desc="End of spatial domain."), N=(val=1000, desc="Number of spatial grid point values used."))
@def funcs = ["generateTable()", "removeTable()", "generate3DPlot()", "remove3DPlot()"];
@def labels = ["Tabulate the solution", "Remove the solution table", "Generate temperature distribution against time plot"];

The heat equation is one of the simplest partial differential equations, it can be expressed as

\begin{align*}
\dfrac{\partial u}{\partial t} &= \Delta u
\end{align*}

where $\Delta$ is the Laplacian. 

~~~
    {{ insert template.html}}
~~~