/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param alpha      Interaction parameter.
 * @param gamma      Interaction parameter.
 * @param t          Time (seconds).
 * @param x          x value.
 * @param y          y value.
 * @param z          z value.
 * @return           [dx/dt, dy/dt, dz/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {alpha, gamma} = objectOfInputs;
    var [x, y, z] = vars;
    return [dt*(y*(z-1+x**2) + gamma*x), dt*(x*(3*z+1-x**2)+gamma*y), dt*(-2*z*(alpha+x*y))];
}