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
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generate3DPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y, z] = vars;

    // Generate 3D phase plot
    gen3DPlot(x, y, z, "phasePlotXYZ", "3D phase plot");
}

/**
 * Generates a XY phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXYPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var y = vars[1];

    // Generate phase plot
    gen2DPlot(x, y, "phasePlotXY", "y against x phase plot");
}

/**
 * Generates a XZ phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var z = vars[2];
    
    // Generate phase plot
    gen2DPlot(x, z, "phasePlotXZ", "z against x phase plot");
}

/**
 * Generates a YZ phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateYZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var y = vars[1];
    var z = vars[2];

    // Generate 2D phase plot
    gen2DPlot(y, z, "phasePlotYZ", "z against y phase plot");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "y", "z"], "timePlot", "Plot of x, y and z against time");
}

/**
 * Generate five plots:
 * - The first is a 3D phase plot of x, y and z.
 * - The second is a 2D phase plot of y against x.
 * - The third is a 2D phase plot of z against x.
 * - The fourth is a 2D phase plot of z against y.
 * - The fifth is a plot of x, y and z against time.
 * 
 * @param objectOfInputs An object containing the parameters in the form.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generate3DPhasePlot(solution);
    generateXYPhasePlot(solution);
    generateXZPhasePlot(solution);
    generateYZPhasePlot(solution);
    generateTimePlot(solution);
}