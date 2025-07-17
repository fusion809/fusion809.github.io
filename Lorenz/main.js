/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param objectOfInputs Interaction parameter.
 * @param t              Time (seconds).
 * @param vars           Solution variables.
 * @return               [dx/dt, dy/dt, dz/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {sigma, rho, beta} = objectOfInputs;
    var [x, y, z] = vars;
    return [dt*sigma*(y-x), dt*(x*(rho-z) - y), dt*(x*y-beta*z)];
}