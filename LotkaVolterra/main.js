/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param alpha      Interaction parameter.
 * @param beta       Interaction parameter.
 * @param gamma      Interaction parameter.
 * @param delta      Interaction parameter.
 * @param t          Time (seconds).
 * @param x          Prey population.
 * @param y          Predator population.
 * @return           [dx/dt, dy/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {alpha, beta, gamma, delta} = objectOfInputs;
    var [x, y] = vars;
    return [dt*(alpha*x - beta*x*y), dt*(delta*x*y-gamma*y)];
}

/** 
 * Solve the problem using RKF45
 *
 * @param solution       An object containing solution data.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions from object and write to 2d array
    var {x0, y0} = objectOfInputs;
    var vars0 = [[x0, y0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generate phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the phase plot.
 */
function generatePhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y] = vars;

    // Generate 2D phase plot
    gen2DPlot(x, y, "phasePlot", "Phase plot of y against x")
}

/**
 * Generate a plot of x and y against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates a plot of x and y against time.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "y"], "timePlot", "Plot of x and y against time");
}

/**
 * Generate two plots:
 * - one of y and x against t; and
 * - a phase plot of y against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generatePhasePlot(solution);
    generateTimePlot(solution);
}