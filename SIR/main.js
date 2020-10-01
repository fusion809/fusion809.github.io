/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
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
    // Solve the problem and extract relevant data
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var [S, I, R] = vars;

    // Generate 3D plot
    gen3DPlot(S, I, R, "phasePlotXYZ", "Phase plot of the solution to the SIR equations. x = S, y = I and z = R");
}

/**
 * Generates a XY phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateXYPhasePlot(objectOfInputs) {
    // Solve the problem and extract relevant data
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var S = vars[0];
    var I = vars[1];

    // Generate 2D plot
    gen2DPlot(S, I, "phasePlotXY", "SI phase plot, x = S and y = I")
}

/**
 * Generates a XZ phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateXZPhasePlot(objectOfInputs) {
    // Solve the problem and extract relevant data
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var S = vars[0];
    var R = vars[2];
    
    // Generate 2D plot
    gen2DPlot(S, R, "phasePlotXZ", "SR phase plot, x = S and y = R");
}

/**
 * Generates a YZ phase plot
 * 
 * @param objectOfInputs An object containing all the form parameters. 
 * @return               Nothing.
 */
function generateYZPhasePlot(objectOfInputs) {
    // Solve the problem and extract relevant data
    var solution = solveProblem(RKF45, objectOfInputs);
    var {vars} = solution;
    var I = vars[1];
    var R = vars[2];

    // Generate 2D plot
    gen2DPlot(I, R, "phasePlotYZ", "IR phase plot, x = I and y = R");
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

    // Generate time plot
    genMultPlot(solution, ["S", "I", "R"], "timePlot", "Plot of SIR against time");
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