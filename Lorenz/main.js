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

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions from object and write to 2d array
    var {x0, y0, z0} = objectOfInputs;
    var vars0 = [[x0, y0, z0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generates a 3D phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generate3DPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y, z] = vars;

    gen3DPlot(x, y, z, "phasePlotXYZ", "Lorenz attractor phase plot.", [2, 0, 0])
}

/**
 * Generates a XY phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateXYPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y] = vars;

    // Generate 2D plot
    gen2DPlot(x, y, "phasePlotXY", "Lorenz system: y against x phase plot.", "x", "y");
}

/**
 * Generates a XZ phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateXZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var z = vars[2];
    
    // Generate 2D plot
    gen2DPlot(x, z, "phasePlotXZ", "Lorenz system: z against x phase plot.", "x", "z");
}

/**
 * Generates a YZ phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateYZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var y = vars[1];
    var z = vars[2];

    // Generate 2D plot
    gen2DPlot(y, z, "phasePlotYZ", "Lorenz system: z against y phase plot.", "y", "z");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "y", "z"], "timePlot", "Lorenz system: plot of x, y and z against time.")
}

/**
 * Generate five plots:
 * - The first is a 3D phase plot of x, y and z.
 * - The second is a 2D phase plot of y against x.
 * - The third is a 2D phase plot of z against x.
 * - The fourth is a 2D phase plot of z against y.
 * - The fifth is a plot of x, y and z against time.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generate3DPhasePlot(solution);
    generateXYPhasePlot(solution);
    generateXZPhasePlot(solution);
    generateYZPhasePlot(solution);
    generateTimePlot(solution);
}

/**
 * Remove animation
 * @return nothing.
 */
function removeAnimation() {
    rmPlot("animation");
}

/**
 * Generate animation
 * @return nothing.
 */
function generateAnimation() {
    var solution = solveProblem(RKF45, readInputs());
    animate3D({solution=solution, view=[2, 0, 0], title="Lorenz system: X, Y and Z phase plot."} = {});
}