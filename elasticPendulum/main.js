/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [x, xDot, theta, thetaDot]
 * @return               [dx/dt, d2x/dt2, dtheta/dt, d2theta/dt2]
 */
function f(objectOfInputs, t, vars) {
    var {g, l0, k, m} = objectOfInputs;
    var [x, xDot, theta, thetaDot] = vars;
    var xDDot = (l0+x)*thetaDot**2 - k*x/m + g*Math.sin(theta);
    var thetaDDot = -g*Math.cos(theta)/(l0+x)-2*xDot*thetaDot/(l0+x);
    return [xDot, xDDot, thetaDot, thetaDDot];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {x0, xDot0, theta0, thetaDot0} = objectOfInputs;
    var vars0 = [[x0, xDot0, theta0, thetaDot0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0); 
    return [t, vars];
}

/**
 * Generates a 2D phase plot of x against theta
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXThetaPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var theta = vars[2];

    // Generate 2D plot
    gen2DPlot(x, theta, "phasePlotXTheta", "Phase plot of theta against x");
}

/**
 * Generates a xdot against x phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXXDotPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var xDot = vars[1];

    // Generate 2D plot
    gen2DPlot(x, xDot, "phasePlotXXDot", "Phase plot of x dot against x");
}

/**
 * Generates a theta dot against theta phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateThetaThetaDotPhasePlot(solution) {    
    // Extract solution data from solution object
    var {vars} = solution;
    var theta = vars[2];
    var thetaDot = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta, thetaDot, "phasePlotThetaThetaDot", "Phase plot of theta dot against theta");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "x dot", "theta", "theta dot"], "timePlot", "Plot of x, x dot, theta and theta dot against time");
}

/**
 * Generate four plots:
 * - The first is a phase plot of theta against x.
 * - The second is a phase plot of x dot against x.
 * - The third is a phase plot of theta dot against theta.
 * - The fourth is a plot of x, x dot, theta and theta dot against time.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Plot solution
    generateXThetaPhasePlot(solution);
    generateXXDotPhasePlot(solution);
    generateThetaThetaDotPhasePlot(solution);
    generateTimePlot(solution);
}
