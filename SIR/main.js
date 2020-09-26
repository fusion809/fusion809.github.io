/**
 * Right-hand side of our ODE written as a simple of first-order ODEs.
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of dependent variables.
 * @return               An array of derivatives.
 */
function f(objectOfInputs, t, vars) {
    var [S, I, R] = vars;
    var {beta, gamma, delta} = objectOfInputs;
    // Determine N
    var N = S+I+R;
    // Calculate derivatives
    var dSdt = - beta * S * I * (1-delta)/N;
    var dIdt = beta * S * I * (1-delta)/N - gamma * I;
    var dRdt = gamma*I;
    // Put into return value
    return [dSdt, dIdt, dRdt];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions from object and write to 2d array
    var {S0, I0, R0} = objectOfInputs;
    var vars0 = [[S0, I0, R0]];
    var [t, vars] = RKF45Body(objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generates a 3D phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generate3DPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var [S, I, R] = vars;

    // Height and width of plot
    adjustPlotHeight("phasePlotXYZ");

    // Plot object and data object array
    var plotXYZ = {
        x: S,
        y: I,
        z: R,
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
        title: 'Phase plot of the solution to the SIR equations. x = S, y = I and z = R'
    };

    // Generate plot
    Plotly.newPlot('phasePlotXYZ', dataXYZ, layoutXYZ);
}

/**
 * Generates a XY phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateXYPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var S = vars[0];
    var I = vars[1];

    // Height and width of plot
    adjustPlotHeight("phasePlotXY");

    // Plot object and data object array
    var plotXY = {
        x: S,
        y: I,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXY = [plotXY];

    // layout object
    var layoutXY = {
        title: "SI phase plot, x = S and y = I"
    };

    // Generate plot
    Plotly.newPlot('phasePlotXY', dataXY, layoutXY);
}

/**
 * Generates a XZ phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateXZPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var S = vars[0];
    var R = vars[2];
    
    // Height and width of plot
    adjustPlotHeight("phasePlotXZ");
    
    // Plot object and data object array
    var plotXZ = {
        x: S,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXZ = [plotXZ];
    
    // layout object
    var layoutXZ = {
        title: "SR phase plot, x = S and y = R"
    };
    
    // Generate plot
    Plotly.newPlot('phasePlotXZ', dataXZ, layoutXZ);
}

/**
 * Generates a YZ phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateYZPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var I = vars[1];
    var R = vars[2];

    // Height and width of plot
    adjustPlotHeight("phasePlotYZ");

    // Plot object and data object array
    var plotYZ = {
        x: I,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataYZ = [plotYZ];

    // layout object
    var layoutYZ = {
        title: "IR phase plot, x = I and y = R"
    };

    // Generate plot
    Plotly.newPlot('phasePlotYZ', dataYZ, layoutYZ);
}

/**
 * Generates a time plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateTimePlot(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);
    var {t, vars} = solution;
    var [S, I, R] = vars;

    // Height and width of plot
    adjustPlotHeight("timePlot");

    // Plot object and data object array
    var plotTX = {
        x: t,
        y: S,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'S'
    };
    var plotTY = {
        x: t,
        y: I,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'I'
    };
    var plotTZ = {
        x: t,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'R'
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
 * - The first is a 3D phase plot of S, I and R.
 * - The second is a 2D phase plot of I against S.
 * - The third is a 2D phase plot of R against S.
 * - The fourth is a 2D phase plot of R against I.
 * - The fifth is a plot of S, I and R against time.
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generate3DPhasePlot(objectOfInputs);
    generateXYPhasePlot(objectOfInputs);
    generateXZPhasePlot(objectOfInputs);
    generateYZPhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
};