/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param objectOfInputs Problem parameter.
 * @param t              Time (seconds).
 * @return               [dx, d2x, dy, d2y]
 */
function f(objectOfInputs, t, vars, dt) {
    var {alpha, beta, gamma, delta, omega, kappa} = objectOfInputs;
    var [x, dx, y, dy] = vars;
    var d2x = - delta*dx - alpha*x - beta*x**3 + gamma * Math.cos(omega*t) + kappa * (y-x);
    var d2y = - delta*dy - alpha*y - beta*y**3 + kappa * (x-y);
    return [dt*dx, dt*d2x, dt*dy, dt*d2y];
}

/**
 * Generate phase plot of x dot against x
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateXPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, dx] = vars;

    // Generate plot
    gen2DPlot(x, dx, "xPhasePlot", "Phase plot of dx/dt against x", "x", "dx/dt");
}

/**
 * Generate phase plot of y dot against y
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateYPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, dx, y, dy] = vars;

    // Generate plot
    gen2DPlot(y, dy, "yPhasePlot", "Phase plot of dy/dt against y", "y", "dy/dt");
}

/**
 * Generate phase plot of y dot against y
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateXYPlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, dx, y, dy] = vars;

    // Generate plot
    gen2DPlot(x, y, "XYPlot", "Phase plot of y against x", "x", "y");
}

/**
 * Generate plot of x and x dot against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(solution) {
    // Plot
    genMultPlot(solution, ["x", "dx/dt", "y", "dy/dt"], "timePlot", "Plot of dx/dt, x, dy/dt and y against time");
}

/**
 * Generate two plots:
 * - one of dx and x against t; and
 * - a phase plot of dx against x.
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
    generateTimePlot(solution);
    generatePhasePlot(solution);
}

/**
 * Generate animation
 * @return nothing
 */
function generateXAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["x", "dx/dt"], IdSuffix: "X", nos: [0, 1], title: "Coupled Duffing system: phase plot of dx/dt against x."});
}

/**
 * Generate animation
 * @return nothing
 */
function generateYAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["y", "dy/dt"], IdSuffix: "Y", nos: [2, 3], title: "Coupled Duffing system: phase plot of dy/dt against y."});
}

/**
 * Generate animation
 * @return nothing
 */
function generateAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["x", "y"], nos: [0, 2], IdSuffix: "XY", title: "Coupled Duffing system: plot of y against x."});
}

function generateAnimations(objectOfInputs=undefined, solution=undefined) {
    generateXAnimation(objectOfInputs, solution);
    generateYAnimation(objectOfInputs, solution);
    generateAnimation(objectOfInputs, solution);
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    fillTable(objectOfInputs, ['x', 'dx/dt', 'y', 'dy/dt'], solution)
}