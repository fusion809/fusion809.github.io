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
function f(objectOfInputs, t, vars) {
    var {alpha, beta, gamma, delta} = objectOfInputs;
    var [x, y] = vars;
    return [alpha*x -beta*x*y, delta*x*y-gamma*y];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract data from object
    var {x0, y0} = objectOfInputs;

    // Initialize the arrays used and loop variables
    var vars0 = [[x0, y0]];
    var [t, vars] = RKF45Body(objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generate phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the phase plot.
 */
function generatePhasePlot(objectOfInputs) {
    // Solve the problem 
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y] = vars;

    // Height and width of plot
    adjustPlotHeight("phasePlot");

    // Characteristics of the phase plot
    var plot = {
        x: x,
        y: y,
        type: 'scatter',
        name: "Phase plot"
    };
    var layout = {
        title: "Phase plot of y against x",
        xaxis: {
            title: "x"
        },
        yaxis: {
            title: "y"
        }
    };
    var data = [plot];

    // Generate plot
    Plotly.newPlot('phasePlot', data, layout);    
}

/**
 * Generate a plot of x and y against time
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates a plot of x and y against time.
 */
function generateTimePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {t, vars} = solution;
    var [x, y] = vars;

    // Height and width of plots
    adjustPlotHeight("timePlot");

    // Characteristics of the x and y against time plot
    var plotx = {
        x: t,
        y: x,
        type: 'scatter',
        name: "x"
    };
    var ploty = {
        x: t,
        y: y,
        type: 'scatter',
        name: "y"
    };
    var layout = {
        title: 'y and x against time plots',
        xaxis: {
           title: 'Time (seconds)'
        }
    };
    var data = [plotx, ploty];

    // Generate plots
    Plotly.newPlot('timePlot', data, layout);
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
    generatePhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
}