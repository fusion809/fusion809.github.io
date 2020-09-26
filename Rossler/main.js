/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param a          Interaction parameter.
 * @param b          Interaction parameter.
 * @param c          Interaction parameter.
 * @param t          Time (seconds).
 * @param x          x value.
 * @param y          y value.
 * @param z          z value.
 * @return           [dx/dt, dy/dt, dz/dt]
 */
function f(objectOfInputs, t, vars) {
    var {a, b, c} = objectOfInputs;
    var [x, y, z] = vars;
    return [-y-z, x + a*y, b + z*(x-c)];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract data from object
    var {x0, y0, z0} = objectOfInputs;

    // Initialize the arrays used and loop variables
    var vars0 = [[x0, y0, z0]];
    var [t, vars] = RKF45Body(objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generates a 3D phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generate3DPhasePlot(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y, z] = vars;

    // Height and width of plot
    adjustPlotHeight("phasePlotXYZ");

    // Plot object and data object array
    var plotXYZ = {
        x: x,
        y: y,
        z: z,
        type: 'scatter3d',
        mode: 'lines',
        opacity: 1,
        line: {
            width: 6,
            reversescale: false
        }
    };
    var dataXYZ = [plotXYZ];

    // layout object
    var layoutXYZ = {
        title: 'Phase plot of the solution to the Rossler equations'
    };

    // Generate plot
    Plotly.newPlot('phasePlotXYZ', dataXYZ, layoutXYZ);
}

/**
 * Generates a XY phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateXYPhasePlot(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var y = vars[1];

    // Height and width of plot
    adjustPlotHeight("phasePlotXY");

    // Plot object and data object array
    var plotXY = {
        x: x,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXY = [plotXY];

    // layout object
    var layoutXY = {
        title: "xy phase plot"
    };

    // Generate plot
    Plotly.newPlot('phasePlotXY', dataXY, layoutXY);
}

/**
 * Generates a XZ phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateXZPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var z = vars[2];
    
    // Height and width of plot
    adjustPlotHeight("phasePlotXZ");
    
    // Plot object and data object array
    var plotXZ = {
        x: x,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXZ = [plotXZ];
    
    // layout object
    var layoutXZ = {
        title: "xz phase plot"
    };
    
    // Generate plot
    Plotly.newPlot('phasePlotXZ', dataXZ, layoutXZ);
}

/**
 * Generates a YZ phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateYZPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {vars} = solution;
    var y = vars[1];
    var z = vars[2];

    // Height and width of plot
    adjustPlotHeight("phasePlotYZ");

    // Plot object and data object array
    var plotYZ = {
        x: y,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataYZ = [plotYZ];

    // layout object
    var layoutYZ = {
        title: "yz phase plot"
    };

    // Generate plot
    Plotly.newPlot('phasePlotYZ', dataYZ, layoutYZ);
}

/**
 * Generates a time plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateTimePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Extract solution data from solution object
    var {t, vars} = solution;
    var [x, y, z] = vars;

    // Height and width of plot
    adjustPlotHeight("timePlot");

    // Plot object and data object array
    var plotTX = {
        x: t,
        y: x,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'x'
    };
    var plotTY = {
        x: t,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'y'
    };
    var plotTZ = {
        x: t,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'z'
    };
    var dataTimePlot = [plotTX, plotTY, plotTZ];

    // layout object
    var layoutTimePlot = {
        title: "Time plots of the solution to the problem"
    };

    // Generate plot
    Plotly.newPlot('timePlot', dataTimePlot, layoutTimePlot);
}

/**
 * Generate five plots:
 * - The first is a 3D phase plot of x, y and z.
 * - The second is a 2D phase plot of y against x.
 * - The third is a 2D phase plot of z against x.
 * - The fourth is a 2D phase plot of z against y.
 * - The fifth is a plot of x, y and z against time.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generate3DPhasePlot(objectOfInputs);
    generateXYPhasePlot(objectOfInputs);
    generateXZPhasePlot(objectOfInputs);
    generateYZPhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
}