@def title = "FunctionIntegrator.jl"
@def tags = ["packages", "julia"]

![CompatHelper](https://github.com/fusion809/FunctionIntegrator.jl/workflows/CompatHelper/badge.svg)
![TagBot](https://github.com/fusion809/FunctionIntegrator.jl/workflows/TagBot/badge.svg)
![Travis](https://travis-ci.com/fusion809/FunctionIntegrator.jl.svg?branch=master)

[FunctionIntegrator.jl](https://github.com/fusion809/FunctionIntegrator.jl) should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not.

This package provides the following functions:

* `adaptive_simpsons_rule(f::Function, a::Number, b::Number, ε::Float64)`
* `chebyshev_quadrature(f::Function, N::Number, k::Integer, a::Number, b::Number)`
* `hermite_quadrature(f::Function, N::Number, k::Integer)`
* `jacobi_quadrature(f::Function, N::Number, α::Number, β::Number, a::Number, b::Number)`
* `laguerre_quadrature(f::Function, N::Number, k::Integer)`
* `legendre_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `lobatto_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `radau_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `rectangle_rule_left(f::Function, N::Number, a::Number, b::Number)`
* `rectangle_rule_midpoint(f::Function, N::Number, a::Number, b::Number)`
* `rectangle_rule_right(f::Function, N::Number, a::Number, b::Number)`
* `rombergs_method(f::Function, N::Number, a::Number, b::Number)`
* `simpsons_rule(f::Function, N::Number, a::Number, b::Number)`
* `simpsons38_rule(f::Function, N::Number, a::Number, b::Number)`
* `trapezoidal_rule(f::Function, N::Number, a::Number, b::Number)`

use Julia's help function (e.g. by typing `?chebyshev_quadrature`) to find out usage information, should you need it. The choice of function table below also explains some of the details of each of these functions, such as their arguments. 

This package is currently in Julia's General registry, thus to install it one merely needs to run:

```julia-repl
(v1.4) pkg> add FunctionIntegrator
```

and import it using:

```julia-repl
julia> using FunctionIntegrator
```

\toc

## Choice of function
As a general rule of thumb, `adaptive_simpsons_rule` should be the function you use when you are unsure which function to use, as its accuracy is more predictable. The main time when `adaptive_simpsons_rule` should be avoided is when there are unremovable singularities at the endpoints of the domain of integration, in which case using `chebyshev_quadrature` with $k=1$ or using `legendre_quadrature` is likely best.

Each of the functions whose name ends with `_quadrature` uses [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature) is the specifics of which differ between functions. Out of them, `legendre_quadrature` is perhaps the best to go with when you are uncertain which out of the `_quadrature` functions to go with as it can be applied to any grid and any problem and *usually* arrives at a result that is accurate to at least 7 decimal places provided $N \geq 1 \times 10^{4}$. `legendre_quadrature` has the added benefit of being perhaps the fastest `_quadrature` function after controlling for accuracy of the result obtained.

The [test/](https://github.com/fusion809/FunctionIntegrator.jl/tree/master/test/) folder has test scripts that approximate various different integrals (each file, except [test/runtests.jl](https://github.com/fusion809/FunctionIntegrator.jl/blob/master/test/runtests.jl), pertains to a different integral); the N values given are the smallest possible to pass each of the tests listed (except when the test involves a less than (<) sign). If you want to know which function to use for which integral is these tests may be useful as a rough guide.

~~~
<figure float="center">
    <img src="/assets/Root_mean_square_computation_time_FunctionIntegrator.jl.png" width="100%">
    <figcaption><b>Figure 1: a bar graph to compare the computation times for each fo the flexible-domain methods provided by the package except for <code>rectange_rule_left</code>, <code>rectangle_rule_right</code> and <code>rombergs_method</code>.<sup>1</sup> Data is from <a href="https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/359262615" link="_blank">this build</a>.</b></caption>
</figure>
~~~

| Function               | Domain~~~<sup>2</sup>~~~ | Weight | Arguments | Notes |
|------------------------|--------|--------|-----------|-------|
| `adaptive_simpsons_rule` | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`ε` is the relative tolerance in the solution (that is, the maximum amount of relative error in the solution you are willing to tolerate).~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [adaptive Simpson's rule](https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method) with the Lyness criterion. |
| `chebyshev_quadrature` | $[a,b]$~~~<br/>~~~$-1\leq x \leq 1$ | $k=1$: $\dfrac{1}{\sqrt{1-x^2}}$~~~<br/><br/>~~~$k=2$: $\sqrt{1-x^2}$~~~<br/><br/>~~~$k=3$: $\sqrt{\dfrac{1+x}{1-x}}$~~~<br/><br/>~~~$k=4$: $\sqrt{\dfrac{1-x}{1+x}}$ | `f` is the function being integrated. ~~~<br/>~~~`N` is the number of grid points.~~~<br/>~~~`k` is the type of Chebyshev quadrature being used. 1, 2, 3, and 4 refer to the Chebyshev $T_n$, $U_n$, $V_n$ and $W_n$ polynomials respectively.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [Chebyshev-Gauss quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) (note this article does not mention 3rd and 4th quadrature types, corresponding to $k=3$ and $k=4$, respectively). If there are unremovable singularities at the endpoints, it with $k=1$ or `legendre_quadrature` are preferred. |
| `hermite_quadrature`   | $[-\infty,\infty]$ | $e^{-x^2}$ | `f` is the function being integrated.~~~<br/>~~~`N` is the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x^2}$ is assumed to be part of the integrand ($k=2$) or not). | Uses [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature). Only use this if your integration domain is $[-\infty,\infty]$ and your integrand rapidly goes to zero as the absolute value of $x$ gets larger. |
| `jacobi_quadrature`    | $[a,b]$~~~<br/>~~~$-1\leq x \leq 1$ | $(1-x)^{\alpha}(1+x)^{\beta}$ | `f` is the function being integrated.~~~<br/>~~~`N` is the number of grid points.~~~<br/>~~~`α` and `β` are parameters of the weighting function.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature). After controlling for `N` it is one of the slowest quadrature methods. A condition of Gauss-Jacobi quadrature is that $\alpha, \beta \gt -1$. When $\alpha = \beta = 0$, this reduces to `legendre_quadrature`. |
| `laguerre_quadrature`  | $[0,\infty]$ | $e^{-x}$ | `f` is the function being integrated.~~~<br/>~~~`N` is the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x}$ is assumed to be part of the integrand ($k=2$) or not). | Uses [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature). Only use this if your integration domain is $[0,\infty]$ and your integrand rapidly goes to $0$ as $x$ gets larger. |
| `legendre_quadrature`  | $[a,b]$~~~<br/>~~~$-1\leq x \leq 1$ | 1 | `f` is the function being integrated.~~~<br/>~~~`N` is the number of grid points.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature). Generally, this is the best `_quadrature` function to go with when you are otherwise unsure which to go with. |
| `lobatto_quadrature`   | $[a,b]$~~~<br/>~~~$-1\leq x \leq 1$ | 1 | `f` is the function being integrated.~~~<br/>~~~`N` is the number of grid points.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [Gauss-Lobatto quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Lobatto_rules). This function includes, in the calculation is the values of the integrand at one of the endpoints. Consequently, if there are unremovable singularities at the endpoints, this function may fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `radau_quadrature`     | $[a,b]$~~~<br/>~~~$-1\leq x \leq 1$ | 1 | `f` is the function being integrated.~~~<br/>~~~`N` is the number of grid points.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses Gauss–Radau quadrature, for which there is no Wikipedia article is the best article (simplest) I could find on it are [these lecture notes](https://web.archive.org/web/20200628202423/http://www.math.hkbu.edu.hk/ICM/LecturesAndSeminars/08OctMaterials/2/Slide3.pdf). This function includes, in the calculation is the values of the function at the endpoints. Consequently, if there are unremovable singularities at either or both of the endpoints, this function will fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `rectangle_rule_left`  | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`N`. $N$ is the number of grid points used in the integration.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses the rectangle rule, specifically the [left Riemann sum](https://en.wikipedia.org/wiki/Riemann_sum#Left_Riemann_sum). Usually this or `rectangle_rule_right` is the least accurate method. In fact, many of the tests in the FunctionIntegrator.jl repository fail to get accuracy to 7 significant figures with `rectangle_rule_left` with any practically viable value of `N`. |
| `rectangle_rule_midpoint`  | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`N`. $N$ is the number of grid points used in the integration.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses the rectangle rule, specifically the [Riemann midpoint rule](https://en.wikipedia.org/wiki/Riemann_sum#Midpoint_rule). Usually this is more accurate than `rectangle_rule_left` and `rectangle_rule_right` and sometimes rivals `trapezoidal_rule` for accuracy. Interestingly, going by my Travis tests it appears to be even more efficient than `simpsons_rule`. |
| `rectangle_rule_right`  | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`N`. $N$ is the number of grid points used in the integration.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses the rectangle rule, specifically the [right Riemann sum](https://en.wikipedia.org/wiki/Riemann_sum#Right_Riemann_sum). Usually this or `rectangle_rule_left` is the least accurate method. In fact, many of the tests in the FunctionIntegrator.jl repository fail to get accuracy to 7 significant figures with `rectangle_rule_right` with any practically viable value of `N`. |
| `rombergs_method`      | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`N`. Equal to $n$ in [this article](https://en.wikipedia.org/wiki/Romberg's_method) that describes the method.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [Romberg's method](https://en.wikipedia.org/wiki/Romberg's_method). If your PC has only 16 GB RAM, you should not exceed $N=30$, as otherwise you will likely run out of RAM. For most problems, $N\leq20$ should suffice.
| `simpsons_rule`        | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included. $N$ must be even, otherwise this function will throw an error. ~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule). It is one of the best functions to use when you are unsure which to use, provided there are no unremovable singularities within the integration domain, including the endpoints. |
| `simpsons38_rule`        | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included. $N$ must be divisible by three, as otherwise this function will throw an error. ~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses [Simpson's 3/8 rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule). It is less robust than `simpsons_rule`. |
| `trapezoidal_rule`     | $[a,b]$ | N/A | `f` is the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a` is the start of the domain of integration.~~~<br/>~~~`b` is the end of the domain of integration. | Uses the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule). It has the same caveats as `simpsons_rule`. |

**Notes**:
1. `rectangle_rule_left`, `rectangle_rule_midpoint` and `rectangle_rule_right` are not included because they failed to provide the required level of accuracy for many tests with all practically viable `N` values.
2. The $-1\leq x \leq 1$ refers to the quadrature nodes, which are also referred to in the weighting function column.

## Examples
### adaptive\_simpsons\_rule
In the following example, the following integral (henceforth called integral 1) is being computed:

$$ \int_0^{\frac{\pi}{2}} \cos{x} \hspace{0.1cm} dx$$

and the result compared to the analytical solution of $1$.

```julia:./code/adaptive_simpsons
using FunctionIntegrator
a = abs(adaptive_simpsons_rule(x -> cos(x), 0, pi/2)-1);
show(a)
```

\output{./code/adaptive_simpsons}

### chebyshev_quadrature
In the following examples, integral 1 is being computing and the result compared to the analytical solution of $1$. The value being printed is difference between the computed solution and the analytical solution.

$k=1$:
```julia:./code/chebyshev1
using FunctionIntegrator
a = abs(chebyshev_quadrature(x -> cos(x), 1000, 1, 0, pi/2) - 1);
show(a)
```
\output{./code/chebyshev1}

$k=2$:
```julia:./code/chebyshev2
using FunctionIntegrator
a = abs(chebyshev_quadrature(x -> cos(x), 1000, 2, 0, pi/2) - 1);
show(a)
```
\output{./code/chebyshev2}

$k=3$:
```julia:./code/chebyshev3
using FunctionIntegrator
a = abs(chebyshev_quadrature(x -> cos(x), 1000, 3, 0, pi/2) - 1);
show(a)
```
\output{./code/chebyshev3}

$k=4$:
```julia:./code/chebyshev4
using FunctionIntegrator
a = abs(chebyshev_quadrature(x -> cos(x), 1000, 4, 0, pi/2) - 1);
show(a)
```
\output{./code/chebyshev4}

### hermite_quadrature
In the following examples is the integral:

$$ \int_{-\infty}^{\infty} x^2 e^{-x^2} dx.$$

the exact solution of which is $\dfrac{\sqrt{\pi}}{2}$, will be approximated using `hermite_quadrature` and the result will be compared to the analytical solution. The plot of the integrand is (truncated to the domain $x\in[-10,10]$):

```julia:./code/hermite_plot
using PyPlot
x = LinRange(-10,10,1001);
clf()
y = (x.^2).*exp.(-x.^2);
PyPlot.plot(x,y)
PyPlot.savefig(joinpath(@OUTPUT, "hermite_plot.png"), dpi=80)
```
\fig{./code/output/hermite_plot.png}

~~~<br>~~~
$k=1$:

```julia:./code/hermite1
using FunctionIntegrator
a = abs(hermite_quadrature(x -> x^2*exp(-x^2), 100, 1)-sqrt(pi)/2)
show(a)
```
\output{./code/hermite1}

$k=2$:

```julia:./code/hermite2
using FunctionIntegrator
a = abs(hermite_quadrature(x -> x^2, 100, 2)-sqrt(pi)/2)
show(a)
```
\output{./code/hermite2}

### jacobi_quadrature
In this example, we will approximate integral 1 and compare the result to the analytical result of $1$. For this example, we are setting both $\alpha$ and $\beta$ to 1.

```julia:./code/jacobi1
using FunctionIntegrator
a = abs(jacobi_quadrature(x -> cos(x), 1000, 1, 1, 0, pi/2)-1);
show(a)
```
\output{./code/jacobi1}

### laguerre_quadrature
In this section is the integral:

$$ \int_0^{\infty} xe^{-x} dx$$

is being approximated and the result compared to the analytical result of 1. The integrand has the following curve:

```julia:./code/laguerre_plot
using PyPlot
x = LinRange(0,15,1001);
clf()
y = x.*exp.(-x);
PyPlot.plot(x,y)
PyPlot.savefig(joinpath(@OUTPUT, "laguerre_plot.png"), dpi=80)
```
\fig{./code/output/laguerre_plot.png}

~~~<br>~~~
$k=1$:
```julia:./code/laguerre1
using FunctionIntegrator
a = abs(laguerre_quadrature(x -> x*exp(-x), 100, 1)-1)
show(a)
```
\output{./code/laguerre1}

$k=2$:
```julia:./code/laguerre2
using FunctionIntegrator
a = abs(laguerre_quadrature(x -> x, 100, 2)-1)
show(a)
```
\output{./code/laguerre2}

### legendre_quadrature
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/legendre
using FunctionIntegrator
a = abs(legendre_quadrature(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/legendre}

### lobatto_quadrature
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/lobatto
using FunctionIntegrator
a = abs(lobatto_quadrature(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/lobatto}

### radau_quadrature
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/radau
using FunctionIntegrator
a = abs(radau_quadrature(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/radau}

### rectangle\_rule\_left
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/rectangle_left
using FunctionIntegrator
a = abs(rectangle_rule_left(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/rectangle_left}

### rectangle\_rule\_midpoint
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/rectangle_midpoint
using FunctionIntegrator
a = abs(rectangle_rule_midpoint(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/rectangle_midpoint}

### rectangle\_rule\_right
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/rectangle_right
using FunctionIntegrator
a = abs(rectangle_rule_right(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/rectangle_right}

### rombergs_method
In this section, integral is being approximated and the result compared to the analytical result.

```julia:./code/rombergs
using FunctionIntegrator
a = abs(rombergs_method(x -> cos(x), 10, 0, pi/2)-1)
show(a)
```

\output{./code/rombergs}

### simpsons_rule
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/simpsons
using FunctionIntegrator
a = abs(simpsons_rule(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/simpsons}

### simpsons38_rule
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/simpsons38
using FunctionIntegrator
a = abs(simpsons38_rule(x -> cos(x), 1002, 0, pi/2)-1)
show(a)
```
\output{./code/simpsons38}

### trapezoidal_rule
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/trapezoidal
using FunctionIntegrator
a = abs(trapezoidal_rule(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/trapezoidal}

```julia
function f
```

## Error analysis
### Code
```julia:./code/error_analysis
using FunctionIntegrator, PyPlot
nc = round.(exp.(1/4*log(10)*(1:24)));
chebyshev1_error         = zeros(length(nc));
chebyshev2_error         = zeros(length(nc));
chebyshev3_error         = zeros(length(nc));
chebyshev4_error         = zeros(length(nc));
jacobi_error             = zeros(length(nc));
laguerre_error           = zeros(length(nc));
legendre_error           = zeros(length(nc));
lobatto_error            = zeros(length(nc));
radau_error              = zeros(length(nc));
rectangle_left_error     = zeros(length(nc));
rectangle_midpoint_error = zeros(length(nc));
rectangle_right_error    = zeros(length(nc));
rombergs_error           = zeros(30);
simpsons_error           = zeros(length(nc));
simpsons38_error         = zeros(length(nc));
trapezoidal_error        = zeros(length(nc));

for i=1:length(nc)
    chebyshev1_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 1, 0, 100)-pi/2);
    chebyshev2_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 2, 0, 100)-pi/2);
    chebyshev3_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 3, 0, 100)-pi/2);
    chebyshev4_error[i]         = abs(chebyshev_quadrature(x -> sech(x), nc[i], 4, 0, 100)-pi/2);
    jacobi_error[i]             = abs(jacobi_quadrature(x -> sech(x), nc[i], 1, 1, 0, 100)-pi/2);
    laguerre_error[i]           = abs(laguerre_quadrature(x -> sech(x), nc[i], 1)-pi/2);
    legendre_error[i]           = abs(legendre_quadrature(x -> sech(x), nc[i], 0, 100)-pi/2);
    lobatto_error[i]            = abs(lobatto_quadrature(x -> sech(x), nc[i], 0, 100)-pi/2);
    radau_error[i]              = abs(radau_quadrature(x -> sech(x), nc[i], 0, 100)-pi/2);
    rectangle_left_error[i]     = abs(rectangle_rule_left(x -> sech(x), nc[i], 0, 100)-pi/2);
    rectangle_midpoint_error[i] = abs(rectangle_rule_midpoint(x -> sech(x), nc[i], 0, 100)-pi/2);
    rectangle_right_error[i]    = abs(rectangle_rule_right(x -> sech(x), nc[i], 0, 100)-pi/2);
    simpsons_error[i]           = abs(simpsons_rule(x -> sech(x), 2*round(nc[i]/2), 0, 100)-pi/2);
    simpsons38_error[i]         = abs(simpsons38_rule(x -> sech(x), 3*round(nc[i]/3), 0, 100)-pi/2);
    trapezoidal_error[i]        = abs(trapezoidal_rule(x -> sech(x), nc[i], 0, 100)-pi/2);
end

for i=1:30
    rombergs_error[i]           = abs(rombergs_method(x -> sech(x), i, 0, 100)-pi/2);
end

PyPlot.figure(1)
PyPlot.clf()
PyPlot.plot(nc, chebyshev1_error);
PyPlot.title("Chebyshev k=1")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev1_error_plot.png"), dpi=80)
PyPlot.figure(2)
PyPlot.clf()
PyPlot.plot(nc, chebyshev2_error); 
PyPlot.title("Chebyshev k=2")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev2_error_plot.png"), dpi=80)
PyPlot.figure(3)
PyPlot.clf()
PyPlot.plot(nc, chebyshev3_error); 
PyPlot.title("Chebyshev k=3")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev3_error_plot.png"), dpi=80)
PyPlot.figure(4)
PyPlot.clf()
PyPlot.plot(nc, chebyshev4_error)
PyPlot.title("Chebyshev k=4")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "chebyshev4_error_plot.png"), dpi=80)
PyPlot.figure(5)
PyPlot.clf()
PyPlot.plot(nc, jacobi_error)
PyPlot.title("Jacobi")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "jacobi_error_plot.png"), dpi=80)
PyPlot.figure(6)
PyPlot.clf()
PyPlot.plot(nc, laguerre_error)
PyPlot.title("Laguerre")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "laguerre_error_plot.png"), dpi=80)
PyPlot.figure(7)
PyPlot.clf()
PyPlot.plot(nc, legendre_error)
PyPlot.title("Legendre")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "legendre_error_plot.png"), dpi=80)
PyPlot.figure(8)
PyPlot.clf()
PyPlot.plot(nc, lobatto_error)
PyPlot.title("Lobatto")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "lobatto_error_plot.png"), dpi=80)
PyPlot.figure(9)
PyPlot.clf()
PyPlot.plot(nc, radau_error)
PyPlot.title("Radau")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "radau_error_plot.png"), dpi=80)
PyPlot.figure(10)
PyPlot.clf()
PyPlot.plot(nc, rectangle_left_error)
PyPlot.title("Rectangle rule (left)")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rectangle_left_error_plot.png"), dpi=80)
PyPlot.figure(11)
PyPlot.clf()
PyPlot.plot(nc, rectangle_midpoint_error)
PyPlot.title("Rectangle rule (midpoint)")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rectangle_midpoint_error_plot.png"), dpi=80)
PyPlot.figure(12)
PyPlot.clf()
PyPlot.plot(nc, rectangle_right_error)
PyPlot.title("Rectangle rule (right)")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rectangle_right_error_plot.png"), dpi=80)
PyPlot.figure(13)
PyPlot.clf()
PyPlot.plot(1:30, rombergs_error)
PyPlot.title("Romberg's method")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "rombergs_error_plot.png"), dpi=80)
PyPlot.figure(14)
PyPlot.clf()
PyPlot.plot(2*round.(nc/2), simpsons_error)
PyPlot.title("Simpson's rule")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "simpsons_error_plot.png"), dpi=80)
PyPlot.figure(15)
PyPlot.clf()
PyPlot.plot(3*round.(nc/3), simpsons38_error)
PyPlot.title("Simpson's 3/8 rule")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "simpsons38_error_plot.png"), dpi=80)
PyPlot.figure(16)
PyPlot.clf()
PyPlot.plot(nc, trapezoidal_error)
PyPlot.title("Trapezoidal rule")
xscale("log")
xlabel("N")
yscale("log")
ylabel("Error")
PyPlot.savefig(joinpath(@OUTPUT, "trapezoidal_error_plot.png"), dpi=80)
```

### Plots
\fig{./code/output/chebyshev1_error_plot.png}
\fig{./code/output/chebyshev2_error_plot.png}
\fig{./code/output/chebyshev3_error_plot.png}
\fig{./code/output/chebyshev4_error_plot.png}
\fig{./code/output/jacobi_error_plot.png}
\fig{./code/output/laguerre_error_plot.png}
\fig{./code/output/legendre_error_plot.png}
\fig{./code/output/lobatto_error_plot.png}
\fig{./code/output/radau_error_plot.png}
\fig{./code/output/rectangle_left_error_plot.png}
\fig{./code/output/rectangle_midpoint_error_plot.png}
\fig{./code/output/rectangle_right_error_plot.png}
\fig{./code/output/rombergs_error_plot.png}
\fig{./code/output/simpsons_error_plot.png}
\fig{./code/output/simpsons38_error_plot.png}
\fig{./code/output/trapezoidal_error_plot.png}

## Acknowledgements
I would like to thank the Julia discourse community for their help through my Julia journey, and I would also like to thank the developers of Julia, as without their work this package would not even be possible (I know, obviously, as this is a Julia package). Likewise, I would also like to thank the developers of [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl), as without their efficient algorithms for finding the nodes and weights for various Gaussian quadrature techniques, several of the functions in this repository would be far less efficient (especially `legendre_quadrature`), or may not even exist. I would also like to thank the developers of [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), as some of their functions are useful for testing the functions in this package.
