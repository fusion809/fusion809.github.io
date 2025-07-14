/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of dependent variables.
 * @return               An array of derivatives.
 */
function f(objectOfInputs, t, vars, dt) {
    var [S, I, R] = vars;
    var {beta, gamma, delta} = objectOfInputs;
    // Determine N
    var N = S+I+R;
    // Calculate derivatives
    var dSdt = - beta * S * I * (1-delta)/N;
    var dIdt = beta * S * I * (1-delta)/N - gamma * I;
    var dRdt = gamma*I;
    // Put into return value
    return [dt*dSdt, dt*dIdt, dt*dRdt];
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
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generates a 3D phase plot
 * 
 * @param solution       An object containing solution data. 
 * @return               Nothing.
 */
function generateSIRPhasePlot(solution) {
    // Extract solution data
    var {vars} = solution;
    var [S, I, R] = vars;

    // Generate 3D plot
    gen3DPlot(S, I, R, "phasePlotSIR", "Phase plot of the solution to the SIR equations.", {xtitle: "Susceptible", ytitle: "Infectious", ztitle: "Recovered"});
}

/**
 * Generates a SI phase plot
 * 
 * @param solution       An object containing solution data. 
 * @return               Nothing.
 */
function generateSIPhasePlot(solution) {
    // Extract solution data
    var {vars} = solution;
    var S = vars[0];
    var I = vars[1];

    // Generate 2D plot
    gen2DPlot(S, I, "phasePlotSI", "SI phase plot, x = S and y = I", "S", "I")
}

/**
 * Generates a SR phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateSRPhasePlot(solution) {
    // Extract solution data
    var {vars} = solution;
    var S = vars[0];
    var R = vars[2];
    
    // Generate 2D plot
    gen2DPlot(S, R, "phasePlotSR", "SR phase plot, x = S and y = R", "S", "R");
}

/**
 * Generates a IR phase plot
 * 
 * @param solution       An object containing solution data. 
 * @return               Nothing.
 */
function generateIRPhasePlot(solution) {
    // Extract solution data
    var {vars} = solution;
    var I = vars[1];
    var R = vars[2];

    // Generate 2D plot
    gen2DPlot(I, R, "phasePlotIR", "IR phase plot, x = I and y = R", "S", "R");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data. 
 * @return               Nothing.
 */
function generateTimePlot(solution) {
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
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generateSIRPhasePlot(solution);
    generateSRPhasePlot(solution);
    generateSRPhasePlot(solution);
    generateIRPhasePlot(solution);
    generateTimePlot(solution);
};

/**
 * Remove animation
 * @return nothing.
 */
function removeAnimation() {
    rmPlot("animation");
}

/**
 * Generate animation
 * @return nothing. 
 */
function generateAnimation() {
    var solution = solveProblem(RKF45, readInputs());
    animate3D(solution, {varnames: ["Susceptible", "Infectious", "Recovered"], title: "Susceptible, infectious and recovered population phase plot for the SIR infectious disease model."});
}

function generateTable() {
    fillTable(readInputs(), ['Susceptible', 'Infectious', 'Recovered'])
}

function removeSRPhasePlot() {
    rmPlot("phasePlotSR");
}
function removeSIPhasePlot() {
    rmPlot("phasePlotSI");
}

function removeIRPhasePlot() {
    rmPlot("phasePlotIR");
}

function removeSIRPhasePlot() {
    rmPlot("phasePlotSIR");
}
function removePlots() {
    removeSIRPhasePlot()
    removeSRPhasePlot()
    removeSIPhasePlot()
    removeIRPhasePlot()
    removeTimePlot();
}