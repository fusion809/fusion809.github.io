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
    gen2DPlot(x, y, "phasePlot", "Phase plot of predator vs prey animals.", "Prey", "Predator")
}

/**
 * Generate a plot of x and y against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates a plot of x and y against time.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "y"], "timePlot", "Plot of prey and predator animal numbers against time");
}

/**
 * Generate two plots:
 * - one of y and x against t; and
 * - a phase plot of y against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }

    // Generate plots
    generatePhasePlot(solution);
    generateTimePlot(solution);
}

/**
 * Generate animation
 * @return nothing.
 */
function generateAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {timer: [0.98, 0.98], varnames: ["Prey animals", "Predator animals"], title: "Phase plot of predator vs prey animals in system."});
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    fillTable(objectOfInputs, ['x', 'y'], solution)
}