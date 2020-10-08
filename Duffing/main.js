/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param objectOfInputs Problem parameter.
 * @param t              Time (seconds).
 * @param x              x coordinate.
 * @param xDot           dx/dt
 * @return               [dx/dt, d2x/dt2]
 */
function f(objectOfInputs, t, vars, dt) {
    var {alpha, beta, gamma, delta, omega} = objectOfInputs;
    var [x, xDot] = vars;
    return [dt*xDot, dt*(- delta*xDot - alpha*x - beta*x**3 + gamma * Math.cos(omega*t))];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions and enter into 2d array
    var {x0, xDot0} = objectOfInputs;
    var vars0 = [[x0, xDot0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generate phase plot of x dot against x
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, xDot] = vars;

    // Generate plot
    gen2DPlot(x, xDot, "phasePlot", "Phase plot of x dot against x");
}

/**
 * Generate plot of x and x dot against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(solution) {
    // Plot
    genMultPlot(solution, ["x", "x dot"], "timePlot", "Plot of x dot and x against time");
}

/**
 * Generate two plots:
 * - one of xDot and x against t; and
 * - a phase plot of xDot against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    
    // Generate plots
    generateTimePlot(solution);
    generatePhasePlot(solution);
}