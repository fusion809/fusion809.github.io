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
function f(objectOfInputs, t, vars) {
    var {alpha, beta, gamma, delta, omega} = objectOfInputs;
    var [x, xDot] = vars;
    return [xDot, - delta*xDot - alpha*x - beta*x**3 + gamma * Math.cos(omega*t)];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract data from object
    var {x0, xDot0} = objectOfInputs;

    // Initialize the arrays used and loop variables
    var vars0 = [[x0, xDot0]];
    var [t, vars] = RKF45Body(objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generate phase plot of x dot against x
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {vars} = solution;
    var [x, xDot] = vars;

    // Height and width of plot
    adjustPlotHeight("phasePlot");

    // Characteristics of the phase plot
    var plot = {
        x: x,
        y: xDot,
        type: 'scatter',
        name: "Phase plot"
    };
    var layout = {
        title: "Phase plot of x dot against x",
        xaxis: {
            title: "x (metres)"
        },
        yaxis: {
            title: "x dot (metres per second)"
        }
    };
    var data = [plot];

    // Generate plot
    Plotly.newPlot('phasePlot', data, layout);
}

/**
 * Generate plot of x and x dot against time
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(objectOfInputs) {
    // Run solveProblem() if previously unrun
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {t, vars} = solution;
    var [x, xDot] = vars;

    // Height and width of plots
    adjustPlotHeight("timePlot");

    // Characteristics of the x and x dot against time plot
    var plotx = {
        x: t,
        y: x,
        type: 'scatter',
        name: "x (metres)"
    };
    var plotxDot = {
        x: t,
        y: xDot,
        type: 'scatter',
        name: "x dot (metres per second)"
    };
    var layout = {
        title: 'x and x dot against time plots',
        xaxis: {
            title: 'Time (seconds)'
        }
    };
    var data = [plotx, plotxDot];

    // Generate plots
    Plotly.newPlot('timePlot', data, layout);
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
    generateTimePlot(objectOfInputs);
    generatePhasePlot(objectOfInputs);
}