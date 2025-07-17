/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param objectOfInputs An object containing all problem parameters.
 * @param t              Time (seconds).
 * @param x              x value.
 * @param y              y value.
 * @param z              z value.
 * @return               [dx/dt, dy/dt, dz/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {a, b, c} = objectOfInputs;
    var [x, y, z] = vars;
    return [dt*a*(y-x), dt*(x*(c-a-z) + c*y), dt*(x*y-b*z)];
}