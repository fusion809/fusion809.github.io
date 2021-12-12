@def title="Root-finding without a calculator"
@def tags = ["maths", "numerical methods"]
@def mintoclevel=1

In this article, I will introduce techniques for approximating square and cube roots. The methods used to derive these techniques are very general and can be used to find the roots of any continuous real-valued function.

\tableofcontents
# Square root
A popular technique for approximating the square root of a number is to use the [tangent line approximation](https://en.wikipedia.org/wiki/Linear_approximation). Say we wish to find $a_s=\sqrt{a}$, and we know that the square root of a neighbouring number $b$ is $b_s$, then we would use the approximation: $\sqrt{a} \approx \sqrt{b} + \dfrac{a-b}{2\sqrt{b}}$. This approximation is called the tangent line approximation because it is derived by using the tangent line for the square root function $\sqrt{x}$ to approximate its values. This approximation can also be derived by using the [Taylor series](https://en.wikipedia.org/wiki/Taylor_series) of the square root function. 

## Derivation of the tangent line approximation formula
If we have a function $f(x)$ that we wish to approximate using a tangent line to at $x=x_0$, then the tangent line equation must be of the form:

\[y = mx + c\]

as it is a straight line. At the point $x=x_0$ the following equations must be satisfied:

\begin{eqnarray}
    y &=& f(x) \\\\
    y' &=& f'(x) \\\\
\end{eqnarray}

or in other words:

\begin{eqnarray}
    mx_0 + c &=& f(x_0) \\\\
    m &=& f'(x_0)
\end{eqnarray}

substituting the second of these equations into the first yields:

\begin{eqnarray}
    f'(x_0)x_0 + c &=& f(x_0) \\\\
    c &=& f(x_0) - f'(x_0) x_0 \\\\
    \therefore y &=& f'(x_0) x + f(x_0) - f'(x_0) x_0 \\\\
    &=& f'(x_0)(x-x_0) + f(x_0).
\end{eqnarray}

Which is the equation of the tangent line. If we let $f(x) = \sqrt{x}$ then:

\begin{eqnarray}
    f'(x) &=& \dfrac{1}{2\sqrt{x}} \\\\
    y &=& \dfrac{1}{2\sqrt{x_0}} (x-x_0) + \sqrt{x_0} \\\\
    &=& \sqrt{x_0} + \dfrac{x-x_0}{2\sqrt{x_0}}.
\end{eqnarray}

If we let $x_0 = b$, a neighbouring number we know the square root of, and $x=a$ then this gives us:

\begin{eqnarray}
    y &=& \sqrt{b} + \dfrac{a-b}{2\sqrt{b}}.
\end{eqnarray}

Which is our approximation for $a_s$. 

## Limitations of the tangent line approximation
While its approximation is often satisfactory, it can be substantially off if $b$ and $a$ differ by a number that is fairly large relative to $b$. 

## Newton's method
In this case, you can use [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method) to refine the approximation. Newton's method is a technique in which we essentially use the tangent line to approximate the zeros of a function. If one of that function's roots happens to be the square root you're searching for, applying Newton's method will have the result of giving you ever improving approximations to the square root. One function that is easy to exactly compute that's root is $a_s$ is $f(x) = x^2-a$. Newton's method then gives us this scheme for approximating $a_s$:

\begin{eqnarray}
x_{n+1} &=& x_n - \dfrac{f(x_n)}{f'(x_n)} \\\\
&=& x_n - \dfrac{x_n^2 - a}{2x_n} \\\\
&=& x_n - \dfrac{x_n}{2} + \dfrac{a}{2x_n} \\\\
&=& \dfrac{x_n}{2} + \dfrac{a}{2x_n}.
\end{eqnarray}

Where $x_{n}$ is our $n$th Newton's method approximation to $a_s$.

## Example 1, finding $\sqrt{2}$ to five decimal places
A classic example where the tangent line approximation is pretty poor is when we are trying to find $\sqrt{2}$. The tangent line approximation gives:

\begin{eqnarray}
    \sqrt{2} &\approx& \sqrt{1} + \dfrac{2-1}{2\sqrt{1}} \\\\
    &=& 1 + \dfrac{1}{2} \\\\
    &=& \dfrac{3}{2} \\\\
    &=& 1.5.
\end{eqnarray}

When of course $\sqrt{2}$ is to five decimal places $1.41421$, and the above approximation is only accurate to one significant figure (or not even that if you round it to the nearest integer). We can improve on this approximation using Newton's method. Our first iteration is:

\begin{eqnarray}
    x_1 &=& \dfrac{x_0}{2} + \dfrac{2}{2x_0} \\\\
    &=& \dfrac{\dfrac{3}{2}}{2} + \dfrac{2}{2\times \dfrac{3}{2}} \\\\
    &=& \dfrac{3}{4} + \dfrac{2}{3} \\\\
    &=& \dfrac{3\times 3}{4\times 3} + \dfrac{2\times 4}{3\times 4} \\\\
    &=& \dfrac{9}{12} + \dfrac{8}{12} \\\\
    &=& \dfrac{17}{12} \\\\
    &\approx& 1.417.
\end{eqnarray}

Now this is clearly not accurate to five decimal places, so another iteration of Newton's method is required. 

\begin{eqnarray}
    x_2 &=& \dfrac{x_1}{2} + \dfrac{2}{2x_1} \\\\
    &=& \dfrac{\dfrac{17}{12}}{2} + \dfrac{1}{\dfrac{17}{12}} \\\\
    &=& \dfrac{17}{24} + \dfrac{12}{17} \\\\
    &=& \dfrac{17 \times 17}{24 \times 17} + \dfrac{12 \times 24}{17 \times 24} \\\\
    &=& \dfrac{289}{408} + \dfrac{288}{408} \\\\
    &=& \dfrac{577}{408} \\
    &\approx& 1.4142157.
\end{eqnarray}

This is accurate to four decimal places, the fifth is not as we would round it up to two, not one. Running Newton's method again yields:

\begin{eqnarray}
    x_3 &=& \dfrac{x_2}{2} + \dfrac{2}{2x_2} \\\\
    &=& \dfrac{\dfrac{577}{408}}{2} + \dfrac{1}{\dfrac{577}{408}} \\\\
    &=& \dfrac{577}{816} + \dfrac{408}{577} \\\\
    &=& \dfrac{577 \times 577}{816 \times 577} + \dfrac{408 \times 816}{577 \times 816} \\\\
    &=& \dfrac{332929}{470832} + \dfrac{332928}{470832} \\\\
    &=& \dfrac{665857}{470832} \\
    &\approx& 1.41421356237469.
\end{eqnarray}

Which is an accurate approximation of $\sqrt{2}$ to 11 decimal places. 

## Example 2, approximating $\sqrt{76}$ to five decimal places
The tangent line approximation yields:

\begin{eqnarray}
    \sqrt{76} &\approx& \sqrt{81} - \dfrac{5}{2\sqrt{81}} \\\\
    &=& 9 - \dfrac{5}{2\times 9} \\\\
    &=& 9 - \dfrac{5}{18} \\\\
    &=& \dfrac{9\times 18}{18} - \dfrac{5}{18} \\\\
    &=& \dfrac{162}{18} - \dfrac{5}{18} \\\\
    &=& \dfrac{157}{18} \\\\
    &\approx& 8.722.
\end{eqnarray}

The actual value of $\sqrt{76}$ is $8.717797887081348$, so this approximation will not suffice. Applying Newton's method yields:

\begin{eqnarray}
    x_1 &=& \dfrac{x_0}{2} + \dfrac{76}{2x_0} \\\\
    &=& \dfrac{157}{18 \times 2} + \dfrac{76}{2 \times \dfrac{157}{18}} \\\\
    &=& \dfrac{157}{36} + \dfrac{38}{\dfrac{157}{18}} \\\\
    &=& \dfrac{157}{36} + \dfrac{38 \times 18}{157} \\\\
    &=& \dfrac{157 \times 157}{36 \times 157} + \dfrac{684 \times 36}{157 \times 36} \\\\
    &=& \dfrac{24649}{5652} + \dfrac{24624}{5652} \\\\
    &=& \dfrac{49273}{5652} \\\\
    &\approx& 8.717799.
\end{eqnarray}

Which is accurate to five decimal places. 

## Example 3, approximating $\sqrt{5}$ to five decimal places
The tangent line approximation yields:

\begin{eqnarray}
    \sqrt{5} &\approx& \sqrt{4} + \dfrac{1}{2\sqrt{4}} \\\\
    &=& 2 + \dfrac{1}{4} \\\\
    &=& \dfrac{9}{4} \\\\
    &=& 2.25.
\end{eqnarray}

The actual value of $\sqrt{5}$ is approximately $2.23606797749979$. The first iteration of Newton's method yields:

\begin{eqnarray}
    x_1 &=& \dfrac{x_0}{2} + \dfrac{5}{2x_0} \\\\
    &=& \dfrac{\dfrac{9}{4}}{2} + \dfrac{5}{2 \times \dfrac{9}{4}} \\\\
    &=& \dfrac{9}{8} + \dfrac{10}{9} \\\\
    &=& \dfrac{9 \times 9 + 8 \times 10}{9 \times 8} \\\\
    &=& \dfrac{161}{72} \\\\
    &=& 2.236\overline{11}.
\end{eqnarray}

So this is only accurate to three decimal places. The next iteration of Newton's method yields:

\begin{eqnarray}
    x_2 &=& \dfrac{x_1}{2} + \dfrac{5}{2x_1} \\\\
    &=& \dfrac{\dfrac{161}{72}}{2} + \dfrac{5}{2 \times \dfrac{161}{72}} \\\\
    &=& \dfrac{161}{144} + \dfrac{5 \times 36}{161} \\\\
    &=& \dfrac{161 \times 161 + 5 \times 36 \times 144}{161 \times 144} \\\\
    &=& \dfrac{25921 + 25920}{23184} \\\\
    &=& \dfrac{51841}{23184} \\\\
    &\approx& 2.23606797791.
\end{eqnarray}

Which is accurate to nine decimal places, or eight if we account for the effect of rounding.

# Cube root
Similarly the cube root of a number can be approximated by a tangent line approximation and refined using Newton's method. If we let $f(x) = \sqrt[3]{x}$, then $f'(x) = \dfrac{1}{3(\sqrt[3]{x})^2}$ and:

\begin{eqnarray}
    y &=& f'(x_0)(x-x_0) + f(x_0) \\\\
    &=& \dfrac{1}{3(\sqrt[3]{x_0})^2} (x-x_0) + \sqrt[3]{x_0} \\\\
    &=& \sqrt[3]{x_0} + \dfrac{x-x_0}{3(\sqrt[3]{x_0})^2}.
\end{eqnarray}

Letting $x_0 = b$, some neighbouring point we know the cube root of exactly, and $x = a$, the number we want to find the cube root of we get (if $a_c = \sqrt[3]{a}$):

\begin{eqnarray}
    a_c &\approx& \sqrt[3]{b} + \dfrac{a-b}{3(\sqrt[3]{b})^2}.
\end{eqnarray}

And to refine this approximation we can use Newton's method with $f(x) = x^3 - a$ and hence $f'(x) = 3x^2$. Therefore we use the formula:

\begin{eqnarray}
    x_{n+1} &=& x_n - \dfrac{f(x_n)}{f'(x_n)} \\\\
    &=& x_n - \dfrac{x_n^3-a}{3x_n^2} \\\\
    &=& x_n - \dfrac{x_n}{3} + \dfrac{a}{3x_n^2} \\\\
    &=& \dfrac{2x_n}{3} + \dfrac{a}{3x_n^2}.
\end{eqnarray}

## Example 4, approximating $\sqrt[3]{63}$ to five decimal places
Here we use $b=64$, as $4^3 = 64$ and hence we know the cube root of 64 is 4. Therefore our tangent line approximation to the cube root is:

\begin{eqnarray}
    \sqrt[3]{63} &\approx& \sqrt[3]{64} - \dfrac{1}{3(\sqrt[3]{64})^2} \\\\
    &=& 4 - \dfrac{1}{3\times 16} \\\\
    &=& 4 - \dfrac{1}{48} \\\\
    &=& \dfrac{191}{48} \\\\
    &\approx& 3.9791\overline{6}.
\end{eqnarray}

The correct value of $\sqrt[3]{63}$ is approximately $3.9790572078963917$, so this is not accurate to five decimal places. So to improve this estimate we apply Newton's method:

\begin{eqnarray}
    x_1 &=& \dfrac{2x_0}{3} + \dfrac{63}{3x_0^2} \\\\
    &=& \dfrac{2x_0}{3} + \dfrac{21}{x_0^2} \\\\
    &=& \dfrac{2 \times \dfrac{191}{48}}{3} + \dfrac{21}{\left(\dfrac{191}{48}\right)^2} \\\\
    &=& \dfrac{191}{72} + \dfrac{21 \times 48^2}{191^2} \\\\
    &=& \dfrac{191}{72} + \dfrac{21 \times 2304}{36481} \\\\
    &=& \dfrac{191 \times 36481 + 21 \times 2304 \times 72}{72 \times 36481} \\\\
    &=& \dfrac{6967871 + 3483648}{2626632} \\\\
    &=& \dfrac{10451519}{2626632} \\\\
    &\approx& 3.97905721.
\end{eqnarray}

Which is accurate to seven decimal places, or eight if we round off the remaining digits from there.

# General technique
Newton's method can be used for finding the roots of any continuous real-valued function, or even system of continuous real-valued functions, although it has some shortcomings. One is that it requires an initial guess as to the root, and that its results can depend heavily on this initial guess. Additionally, it can sometimes fail to converge. For polynomial equations that have real roots this is uncommon provided the initial guess is reasonably close to the true value of the root. 

Getting the initial guesses for the roots, especially of polynomials that have multiple real roots, can be a challenge. One technique is to use the [bisection method](https://en.wikipedia.org/wiki/Bisection_method). In this technique one takes an interval over which the sign of the function changes from positive to negative or vice versa (which implies that at least one root must lie within this interval) and continuously subdivides this interval and re-evaluates the function at its endpoints until we find the root. This method is very slow, so usually one would only apply it a few times, and then once our interval is acceptably small we would apply Newton's method to get the root. One can also evaluate the function at a long list of points within an interval and find where in the interval the sign changes, as that is where a root will be. If the interval is large enough and the function has multiple real roots there may be multiple sign change points and hence multiple roots we can converge to using Newton's method. 

These methods can be done by hand, but for many problems they can be the quite tedious; as such I wrote this Julia script that implements the bisection method and Newton's method to find the roots of a given function. Here is the script:

```julia
#!/usr/bin/env julia
"""
    bisection(f::Function, N::Integer, a::Number, b::Number)

Uses the bisection method with N subdivisions of the specified domain [a, b] 
to find the roots of the function f within said domain. Used to find the 
initial guess that Newton's method then uses to get a more precise estimate of 
the roots.
"""
function bisection(f, N, a, b)
    x = LinRange(a, b, N+1)
    change = zeros(size(x))
    fval = f.(x)
    for i=2:N+1
        if (sign(fval[i]) + sign(fval[i-1]) == 0)
            change[i] = 1
        end
    end

    initGuess = x[change.==1]
    return initGuess
end

"""
    newtons(f::Function, h::Float, tol::Float, itMax::Integer, 
    initGuess::Vector{Float})

Uses Newton's method to refine initGuess, our initial guess of the root(s) of 
f, until either itMax iterations has been performed or the relative magnitude 
of the update Newton's provides to the estimate is below tol. h is the step 
size used to approximate the derivative.  
"""
function newtons(f, h, tol, itMax, initGuess)
    function fd(x, h)
        return (f(x+h)-f(x-h))/(2*h)
    end
    sol = initGuess
    count = Int64.(zeros(size(initGuess)))
    for j=1:length(initGuess)
        diff = 1
        while (abs(diff/sol[j]) > tol && count[j] < itMax)
            diff = f(sol[j])/fd(sol[j], h)
            sol[j] -= diff
            count[j] += 1
        end
        if (count[j] == itMax)
            println("Maximum iterations exceeded and the amount by which 
            Newton's last updated the solution was: ", diff)
        end
    end
    return sol, count
end

"""
    findRoot(f::Function, h::Float, tol::Float, itMax::Integer, a::Number, b::Number, N::Integer)

Uses the bisection method to get an initial guess of the root(s) of f on the 
domain [a, b] with N subdivisions, then uses Newton's method with a maximum of 
itMax iterations and a relative error tolerance of tol. h is the step size 
used to approximate the derivative. 
"""
function findRoot(f, h, tol, itMax, a, b, N)
    initGuess = bisection(f, N, a, b)
    sol, count = newtons(f, h, tol, itMax, initGuess);
    return sol, count
end

# This is where you specify the function you want to
# find the root of
function f(x)
    return x^4 + x^3 - 10x^2 - 4x + 16
end

h = 1e-10
tol = 1e-15
itMax = 1000
a = -100
b = 100
N = 100000
sol, count = findRoot(f, h, tol, itMax, a, b, N)
for i=1:length(sol)
    println("The $(i)th root = ", sol[i], " count = ", count[i])
end
```

.