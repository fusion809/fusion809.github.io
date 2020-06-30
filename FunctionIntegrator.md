@def title = "FunctionIntegrator.jl"
@def tags = ["packages", "julia"]

![CompatHelper](https://github.com/fusion809/FunctionIntegrator.jl/workflows/CompatHelper/badge.svg)
![TagBot](https://github.com/fusion809/FunctionIntegrator.jl/workflows/TagBot/badge.svg)
![Travis](https://travis-ci.com/fusion809/FunctionIntegrator.jl.svg?branch=master)

[FunctionIntegrator.jl](https://github.com/fusion809/FunctionIntegrator.jl) should be treated as a second-rate alternative to the excellent [QuadGK](https://github.com/JuliaMath/QuadGK.jl) package. QuadGK provides more accurate integration for many problems, and also provides an error estimate which functions in this package do not. Likewise this package is less user-friendly as it requires you to decide which $N$ value you are going to go with.

This package provides the following functions:

* `chebyshev_quadrature(f::Function, N::Number, k::Integer, a::Number, b::Number)`
* `hermite_quadrature(f::Function, N::Number, k::Integer)`
* `jacobi_quadrature(f::Function, N::Number, α::Number, β::Number, a::Number, b::Number)`
* `laguerre_quadrature(f::Function, N::Number, k::Integer)`
* `legendre_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `lobatto_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `radau_quadrature(f::Function, N::Number, a::Number, b::Number)`
* `rectangle_rule(f::Function, N::Number, a::Number, b::Number)`
* `simpsons_rule(f::Function, N::Number, a::Number, b::Number)`
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
As a general rule of thumb, `simpsons_rule` should be the function you use when you are unsure which function to use as its approximations with large N (e.g. $1 \times 10^{4}$) are nearly always accurate to at least 6 digits. Despite this, for many problems some of the `_quadrature` functions may provide more accurate results with a smaller N value. The main time when `simpsons_rule` should be avoided is when there are unremovable singularities at the endpoints of the domain of integration, in which case using `chebyshev_quadrature` with $k=1$ or using `legendre_quadrature` is likely best.

Each of the functions whose name ends with `_quadrature` uses [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature), the specifics of which differ between functions. Out of them, `legendre_quadrature` is perhaps the best to go with when you are uncertain which out of the `_quadrature` functions to go with as it can be applied to any grid and any problem and *usually* arrives at a result that is accurate to at least 7 decimal places provided $N \geq 1 times 10^{4}$. `legendre_quadrature` has the added benefit of being perhaps the fastest `_quadrature` function after controlling for accuracy of the result obtained.

The [test/](https://github.com/fusion809/FunctionIntegrator.jl/tree/master/test/) folder has test scripts that approximate various different integrals (each file, except [test/runtests.jl](https://github.com/fusion809/FunctionIntegrator.jl/blob/master/test/runtests.jl), pertains to a different integral); the N values given are the smallest possible to pass each of the tests listed (except when the test involves a less than (<) sign). If you want to know which function to use for which integral, these tests may be useful as a rough guide.

~~~
<figure float="center">
    <img src="/assets/Root_mean_square_computation_time_FunctionIntegrator.jl.png">
    <figcaption><b>Figure 1: a bar graph to compare the computation times for each fo the general-purpose methods provided by the package, except <code>rectange_rule</code>.<sup>1</sup> Data is from <a href="https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/354648732" link="_blank">this build</a>.</b></caption>
</figure>
~~~

| Function               | Domain~~~<sup>2</sup>~~~ | Weight | Arguments | Notes |
|------------------------|--------|--------|-----------|-------|
| `chebyshev_quadrature` | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | $k=1$: $\dfrac{1}{\sqrt{1-x^2}}$~~~<br/><br/>~~~$k=2$: $\sqrt{1-x^2}$~~~<br/><br/>~~~$k=3$: $\sqrt{\dfrac{1+x}{1-x}}$~~~<br/><br/>~~~$k=4$: $\sqrt{\dfrac{1-x}{1+x}}$ | `f`, the function being integrated. ~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k`, the type of Chebyshev quadrature being used. 1, 2, 3, and 4 refer to the Chebyshev $T_n$, $U_n$, $V_n$ and $W_n$ polynomials respectively.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Chebyshev-Gauss quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) (note this article does not mention 3rd and 4th quadrature types, corresponding to $k=3$ and $k=4$, respectively). If there are unremovable singularities at the endpoints, it with $k=1$ or `legendre_quadrature` are preferred. |
| `hermite_quadrature`   | $[-\infty,\infty]$ | $e^{-x^2}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x^2}$ is assumed to be part of the integrand ($k=2$) or not). | Uses [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature). Only use this if your integration domain is $[-\infty,\infty]$ and your integrand rapidly goes to zero as the absolute value of $x$ gets larger. |
| `jacobi_quadrature`    | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | $(1-x)^{\alpha}(1+x)^{\beta}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`α` and `β` are parameters of the weighting function.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature). After controlling for `N` it is one of the slowest quadrature methods. A condition of Gauss-Jacobi quadrature is that $\alpha, \beta \gt -1$. When $\alpha = \beta = 0$, this reduces to `legendre_quadrature`. |
| `laguerre_quadrature`  | $[0,\infty]$ | $e^{-x}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x}$ is assumed to be part of the integrand ($k=2$) or not). | Uses [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature). Only use this if your integration domain is $[0,\infty]$ and your integrand rapidly goes to $0$ as $x$ gets larger. |
| `legendre_quadrature`  | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature). Generally, this is the best `_quadrature` function to go with when you are otherwise unsure which to go with. |
| `lobatto_quadrature`   | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Gauss-Lobatto quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Lobatto_rules). This function includes, in the calculation, the values of the integrand at one of the endpoints. Consequently, if there are unremovable singularities at the endpoints, this function may fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `radau_quadrature`     | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses Gauss–Radau quadrature, for which there is no Wikipedia article, the best article (simplest) I could find on it are [these lecture notes](https://web.archive.org/web/20200628202423/http://www.math.hkbu.edu.hk/ICM/LecturesAndSeminars/08OctMaterials/2/Slide3.pdf). This function includes, in the calculation, the values of the function at the endpoints. Consequently, if there are unremovable singularities at either or both of the endpoints, this function will fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `rectangle_rule`       | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N$ is the number of grid points used in the integration.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses the rectangle rule, specifically the [left Riemann sum](https://en.wikipedia.org/wiki/Riemann_sum#Left_Riemann_sum). Usually this is the least accurate method. In fact, many of the tests in the FunctionIntegrator.jl repository fail to get accuracy to 7 significant figures with `rectangle_rule` with any practically viable value of `N`. |
| `simpsons_rule`        | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule). The best function to use when you are unsure which to use, provided there are no unremovable singularities within the integration domain, including the endpoints. It is the best in the sense of being, on average, the most efficient (providing the most accurate result for the least amount of computing time). |
| `trapezoidal_rule`     | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Uses the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule). This is the second most efficient integration function after `simpsons_rule`. It has the same caveats. |

**Notes**:
1. `rectangle_rule` is not included because it failed to provide the required level of accuracy for many tests with all practically viable `N` values.
2. The $x\in[-1,1]$ refers to the quadrature nodes, which are also referred to in the weighting function column.

## Examples
### chebyshev_quadrature
In the following examples, the following integral (henceforth called integral 1) is being computed:

$\displaystyle \int_0^{\frac{\pi}{2}} \cos{x} \hspace{0.1cm} dx$

and the result compared to the analytical solution of $1$. The value being printed is difference between the computed solution and the analytical solution.

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
In the following examples the integral:

$\displaystyle \int_{-\infty}^{\infty} x^2 e^{-x^2} dx.$

The exact solution of which is $\dfrac{\sqrt{\pi}}{2}$. The plot of the integrand is (truncated to the domain $x\in[-10,10]$):

```julia:./code/hermite_plot
using PyPlot
x = LinRange(-10,10,1001);
clf()
y = (x.^2).*exp.(-x.^2);
PyPlot.plot(x,y)
PyPlot.savefig(joinpath(@OUTPUT, "hermite_plot.png"), dpi=50)
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
In this section, the integral:

$\displaystyle \int_0^{\infty} xe^{-x} dx$

is being approximated and the result compared to the analytical result of 1. The integrand has the following curve:

```julia:./code/laguerre_plot
using PyPlot
x = LinRange(0,15,1001);
clf()
y = x.*exp.(-x);
PyPlot.plot(x,y)
PyPlot.savefig(joinpath(@OUTPUT, "laguerre_plot.png"), dpi=50)
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

### rectangle_rule
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/rectangle
using FunctionIntegrator
a = abs(rectangle_rule(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/rectangle}

### simpsons_rule
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/simpsons
using FunctionIntegrator
a = abs(simpsons_rule(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/simpsons}

### trapezoidal_rule
In this section, integral 1 is being approximated and the result compared to the analytical result. 

```julia:./code/trapezoidal
using FunctionIntegrator
a = abs(trapezoidal_rule(x -> cos(x), 1000, 0, pi/2)-1)
show(a)
```
\output{./code/trapezoidal}

## Acknowledgements
I would like to thank the Julia discourse community for their help through my Julia journey, and I would also like to thank the developers of Julia, as without their work this package would not even be possible (I know, obviously, as this is a Julia package). Likewise, I would also like to thank the developers of [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl), as without their efficient algorithms for finding the nodes and weights for various Gaussian quadrature techniques, several of the functions in this repository would be far less efficient (especially `legendre_quadrature`), or may not even exist. I would also like to thank the developers of [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), as some of their functions are useful for testing the functions in this package.
