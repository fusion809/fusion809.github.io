
using LinearAlgebra

function newton_system()
    # Initial guesses
    x = 2.0
    y = sqrt(5)

    # Iteration parameters
    tol = 1e-12
    maxiter = 100

    for n in 1:maxiter
        s = x + y

        # Residual vector
        F = [
            x*y - log(s),
            x^2 + y^2 - 9
        ]

        # Jacobian matrix
        J = [
            y - 1/s   x - 1/s
            2x        2y
        ]

        # Newton correction
        Δ = J \ F

        # Update solution
        x -= Δ[1]
        y -= Δ[2]

        # Correction magnitude
        iteration_error = norm(Δ, Inf)

        println("Iteration $n: correction magnitude = $iteration_error")

        if iteration_error < tol
            println("\nNewton iteration converged.")
            break
        end

        if n == maxiter
            @warn "Newton iteration did not converge."
        end
    end

    println("\nNumerical solution:")
    println("x = $x")
    println("y = $y")

    return x, y
end

x, y = newton_system()
