/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param objectOfInputs Problem parameter.
 * @param t              Time (seconds).
 * @param x              x coordinate.
 * @param dx           dx/dt
 * @return               [dx/dt, d2x/dt2]
 */
function f(objectOfInputs, t, vars, dt) {
    var {alpha, beta, gamma, delta, omega} = objectOfInputs;
    var [x, dx] = vars;
    return [dt*dx, dt*(- delta*dx - alpha*x - beta*x**3 + gamma * Math.cos(omega*t))];
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
    var [x, dx] = vars;

    // Generate plot
    gen2DPlot(x, dx, "phasePlot", "Phase plot of dx/dt against x", "x", "dx/dt");
}

/**
 * Generate plot of x and x dot against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(solution) {
    // Plot
    genMultPlot(solution, ["x", "dx/dt"], "timePlot", "Plot of dx/dt and x against time");
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
function generateAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["x", "dx/dt"], nos: [0, 1], title: "Duffing system: phase plot of dx/dt against x."});
}

/**
 * Remove animation
 * @return nothing
 */
function removeAnimation() {
  rmPlots("animation");
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    fillTable(objectOfInputs, ['x', 'dx/dt'], solution)
}