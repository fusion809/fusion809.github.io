@def title = "FunctionIntegrator.jl"
@def tags = ["packages", "julia"]

![Travis](https://travis-ci.com/fusion809/FunctionIntegrator.jl.svg?branch=master)
![CompatHelper](https://github.com/fusion809/FunctionIntegrator.jl/workflows/CompatHelper/badge.svg)

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

use Julia's help function (e.g. by typing `?chebyshev_quadrature`) to find out usage information, should you need it.

This package is currently in Julia's General registry, thus to install it one merely needs to run:

```julia-repl
(v1.4) pkg> add FunctionIntegrator
```

and import it using:

```julia-repl
julia> using FunctionIntegrator
```

## Choice of function
As a general rule of thumb, `simpsons_rule` should be the function you use when you are unsure which function to use as its approximations with large N (e.g. $1\times10^{4}$) are nearly always accurate to at least 6 digits. Despite this, for many problems some of the `_quadrature` functions may provide more accurate results with a smaller N value. The main time when `simpsons_rule` should be avoided is when there are singularities at the endpoints of the domain of integration, in which case using `chebyshev_quadrature` with $k=1$ or using `legendre_quadrature` is likely best.

Each of the functions whose name ends with `_quadrature` uses [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature), the specifics of which differ between functions. Out of them, `legendre_quadrature` is perhaps the best to go with when you are uncertain which out of the `_quadrature` functions to go with as it can be applied to any grid and any problem and *usually* arrives at a result that is accurate to at least 6 decimal places provided $N\geq1\times10^{4}$. `legendre_quadrature` has the added benefit of being perhaps the fastest `_quadrature` function after controlling for accuracy of the result obtained. 

The [test/](https://github.com/fusion809/FunctionIntegrator.jl/tree/master/test/) folder has test scripts that approximate various different integrals (each file, except [test/runtests.jl](https://github.com/fusion809/FunctionIntegrator.jl/blob/master/test/runtests.jl), pertains to a different integral); the N values given are the smallest possible to pass each of the tests listed (except when the test involves a less than (<) sign). If you want to know which function to use for which integral, these tests may be useful as a rough guide.

~~~
<figure float="center">
    <img src="/assets/Root_mean_square_computation_time_FunctionIntegrator.jl.png">
    <figcaption><b>Figure 1: a bar graph to compare the computation times for each fo the general-purpose methods provided by the package. Data is from <a href="https://travis-ci.com/github/fusion809/FunctionIntegrator.jl/jobs/354648732" link="_blank">this build</a>.</b><sup>1</sup></caption>
</figure>
~~~

| Function               | Domain~~~<sup>2</sup>~~~ | Weight | Arguments | Notes |
|------------------------|--------|--------|-----------|-------|
| `chebyshev_quadrature` | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | $k=1$: $\dfrac{1}{\sqrt{1-x^2}}$~~~<br/><br/>~~~$k=2$: $\sqrt{1-x^2}$~~~<br/><br/>~~~$k=3$: $\sqrt{\dfrac{1+x}{1-x}}$~~~<br/><br/>~~~$k=4$: $\sqrt{\dfrac{1-x}{1+x}}$ | `f`, the function being integrated. ~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k`, the type of Chebyshev quadrature being used. 1, 2, 3, and 4 refer to the Chebyshev $T_n$, $U_n$, $V_n$ and $W_n$ polynomials respectively.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | If there are unremovable singularities at the endpoints, it with $k=1$ or `legendre_quadrature` are preferred. |
| `hermite_quadrature`   | $[-\infty,\infty]$ | $e^{-x^2}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x^2}$ is assumed to be part of the integrand ($k=2$) or not). | Only use this if your integration domain is $[-\infty,\infty]$ and your integrand converges rapidly. |
| `jacobi_quadrature`    | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | $(1-x)^{\alpha}(1+x)^{\beta}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`α` and `β` are parameters of the weighting function.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | When $\alpha = \beta = 0$, this reduces to `legendre_quadrature`. |
| `laguerre_quadrature`  | $[0,\infty]$ | $e^{-x}$ | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`k` determines the problem being solved (whether $e^{-x}$ is assumed to be part of the integrand ($k=2$) or not). | Only use this if your integration domain is $[0,\infty]$ and your integrand converges rapidly. |
| `legendre_quadrature`  | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Generally, this is the best `_quadrature` function to go with when you're otherwise unsure which to go with. |
| `lobatto_quadrature`   | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | This function includes, in the calculation, the values of the function at the endpoints. Consequently, if there are unremovable singularities at either or both of the endpoints, this function will fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `radau_quadrature`     | $[a,b]$~~~<br/>~~~$x\in[-1,1]$ | 1 | `f`, the function being integrated.~~~<br/>~~~`N`, the number of grid points.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | This function includes, in the calculation, the values of the function at the endpoints. Consequently, if there are unremovable singularities at either or both of the endpoints, this function will fail to give an accurate result even if you adjust the endpoints slightly to avoid the singularities. |
| `rectangle_rule`       | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | Usually this is the least accurate method. |
| `simpsons_rule`        | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | The best function to use when you're unsure which to use, provided there are no unremovable singularities within the integration domain, including the endpoints. It is the best in the sense of being, on average, the most efficient (providing the most accurate result for the least amount of computing time). |
| `trapezoidal_rule`     | $[a,b]$ | N/A | `f`, the function being integrated.~~~<br/>~~~`N`. $N+1$ is the number of grid points, if endpoints are included.~~~<br/>~~~`a`, the start of the domain of integration.~~~<br/>~~~`b`, the end of the domain of integration. | This is the second most efficient integration function after `simpsons_rule`. |

**Notes**:
1. `rectangle_rule` is not included because it failed to provide the required level of accuracy for many tests with all practically viable `N` values.
2. The $x\in[-1,1]$ refers to the quadrature nodes, which are also referred to in the weighting function column.

## Acknowledgements
I'd like to thank the Julia discourse community for their generous help through my Julia journey, and I would also like to thank the developers of  Julia, as without their work this package would not even be possible (I know, obviously, as this is a Julia package). Likewise, I'd also like to thank the developers of [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl), as without their efficient algorithms for finding the nodes and weights for various Gaussian quadrature techniques, several of the functions in this repository would be far less efficient (especially `legendre_quadrature`), or may not even exist. I'd also like to thank the developers of [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl), as some of their functions are useful for testing the functions in this package.
