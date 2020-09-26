/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param a          Inverse of the incubation period.
 * @param beta       Infectivity parameter.
 * @param gamma      Parameter that measures the rate of recovery.
 * @param delta      Parameter that measures the effect of quarantine.
 * @param lambda     Birth rate measured in people per day.
 * @param mu         Death rate measured in people per day. 
 * @param t          Time (seconds).
 * @param S          Susceptible person count.
 * @param E          Exposed person count.
 * @param I          Infectious person count.
 * @param R          Recovered person count.
 * @return           [dS/dt, dE/dt, dI/dt, dR/dt]
 */
function f(objectOfInputs, t, vars) {
    var {a, beta, gamma, delta, lambda, mu} = objectOfInputs;
    var [S, E, I, R] = vars;
    // Determine N
    var N = S + E + I + R;
    // Calculate derivatives
    var exposure = (beta * S * I * (1-delta))/N;
    var dSdt = lambda*N - mu * S - exposure;
    var dEdt = exposure - (mu + a)*E;
    var dIdt = a*E - (mu + gamma) * I;
    var dRdt = gamma*I - mu*R;
    // Put into return value
    return [dSdt, dEdt, dIdt, dRdt];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions from object and write to 2d array
    var {S0, E0, I0, R0} = objectOfInputs;
    var vars0 = [[S0, E0, I0, R0]];
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
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var S = vars[0];
    var I = vars[2];
    var R = vars[3];

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
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var S = vars[0];
    var I = vars[2];

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
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var S = vars[0];
    var R = vars[3];

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
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var I = vars[2];
    var R = vars[3];

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
    // Solve problem and extract relevant solution variables
    var solution = solveProblem(RKF45, objectOfInputs);
    var {t, vars} = solution;
    var [S, E, I, R] = vars;

    // Height and width of plot
    adjustPlotHeight("timePlot");

    // Plot object and data object array
    var plotTS = {
        x: t,
        y: S,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'S'
    };
    var plotTE = {
        x: t,
        y: E,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'E'
    };
    var plotTI = {
        x: t,
        y: I,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'I'
    };
    var plotTR = {
        x: t,
        y: R,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'R'
    };
    var dataTimePlot = [plotTS, plotTE, plotTI, plotTR];

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
 * - The fifth is a plot of S, I and R.
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