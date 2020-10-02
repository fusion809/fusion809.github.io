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
function f(objectOfInputs, t, vars) {
    var {a, b, c} = objectOfInputs;
    var [x, y, z] = vars;
    return [a*(y-x), x*(c-a-z) + c*y, x*y-b*z];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions from object and enter it into RKF45Body
    var {x0, y0, z0} = objectOfInputs;
    var vars0 = [[x0, y0, z0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generates a 3D phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generate3DPhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var vars = solution.vars;
    var [x, y, z] = vars;

    gen3DPlot(x, y, z, "phasePlotXYZ", "Phase plot of the solution to the Chen equations")
}

/**
 * Generates a XY phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateXYPhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var y = vars[1];

    // Generate plot
    gen2DPlot(x, y, "phasePlotXY", "y against x phase plot");
}

/**
 * Generates a XZ phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateXZPhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);
    
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var z = vars[2];

    // Generate plot
    gen2DPlot(x, z, "phasePlotXZ", "z against x phase plot");
}

/**
 * Generates a YZ phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateYZPhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {vars} = solution;
    var y = vars[1];
    var z = vars[2];

    // Generate plot
    gen2DPlot(y, z, "phasePlotYZ", "z against y phase plot");
}

/**
 * Generates a time plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateTimePlot(objectOfInputs) {
    // Run solveProblem if unrun
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    genMultPlot(solution, ["x", "y", "z"], "timePlot", "Time plots of the solution to the problem")
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
    generate3DPhasePlot(objectOfInputs);
    generateXYPhasePlot(objectOfInputs);
    generateXZPhasePlot(objectOfInputs);
    generateYZPhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
}