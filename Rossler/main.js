/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param a          Interaction parameter.
 * @param b          Interaction parameter.
 * @param c          Interaction parameter.
 * @param t          Time (seconds).
 * @param x          x value.
 * @param y          y value.
 * @param z          z value.
 * @return           [dx/dt, dy/dt, dz/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {a, b, c} = objectOfInputs;
    var [x, y, z] = vars;
    return [dt*(-y-z), dt*(x + a*y), dt*(b + z*(x-c))];
}