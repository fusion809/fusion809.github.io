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
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generates a 3D phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generate3DPhasePlot(solution) {
    // Extract relevant solution values
    var {vars} = solution;
    var S = vars[0];
    var I = vars[2];
    var R = vars[3];

    // Generate 3D phase plot
    gen3DPlot(S, I, R, "phasePlotXYZ", "3D phase plot. x: S, y: I and z: R");
}

/**
 * Generates a XY phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXYPhasePlot(solution) {
    // Extract relevant solution values
    var {vars} = solution;
    var S = vars[0];
    var I = vars[2];

    // Generate 2D plot
    gen2DPlot(S, I, "phasePlotXY", "Number of infectious persons against number of susceptible persons");
}

/**
 * Generates a XZ phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXZPhasePlot(solution) {
    // Extract relevant solution variables
    var {vars} = solution;
    var S = vars[0];
    var R = vars[3];

    // Generate 2D plot
    gen2DPlot(S, R, "phasePlotXZ", "Number of recovered persons against number of susceptible persons");
}

/**
 * Generates a YZ phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateYZPhasePlot(solution) {
    // Extract relevant solution values
    var {vars} = solution;
    var I = vars[2];
    var R = vars[3];

    // Generate 2D plot
    gen2DPlot(I, R, "phasePlotYZ", "Number of recovered persons against number of infectious persons");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["S", "E", "I", "R"], "timePlot", "Plot of S, E, I and R against time");
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
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generate3DPhasePlot(solution);
    generateXYPhasePlot(solution);
    generateXZPhasePlot(solution);
    generateYZPhasePlot(solution);
    generateTimePlot(solution);
};