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
function f(objectOfInputs, t, vars, dt) {
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
    return [dt*dSdt, dt*dEdt, dt*dIdt, dt*dRdt];
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
    gen3DPlot(S, I, R, "phasePlotXYZ", "Phase plot of susceptible, infectious and recovered populations.", undefined, "Susceptible", "Infectious", "Recovered");
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
    gen2DPlot(S, I, "phasePlotXY", "Number of infectious persons against number of susceptible persons", "Susceptible", "Infectious");
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
    gen2DPlot(S, R, "phasePlotXZ", "Number of recovered persons against number of susceptible persons", "Susceptible", "Recovered");
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
    gen2DPlot(I, R, "phasePlotYZ", "Number of recovered persons against number of infectious persons", "Infectious", "Recovered");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["Susceptible", "Exposed", "Infectious", "Recovered"], "timePlot", "Plot of susceptible, exposed, infectious and recovered population against time");
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

/**
 * Remove SIR animation
 * @return nothing
 */
function removeAnimationSIR() {
    rmPlot("animationSIR");
}

function generateAnimationBaseSIR(solution) {
    animate3D(solution, [0, 0, 0], ["Susceptible", "Infectious", "Recovered"], [0, 2, 3], "SIR", "Susceptible, infectious and recovered phase plot for the SEIR infectious disease model.");
}
/**
 * Generate SIR animation
 * @return nothing
 */
function generateAnimationSIR() {
    var solution = solveProblem(RKF45, readInputs());
    generateAnimationBaseSIR(solution);
}

/**
 * Remove SEI animation
 * @return nothing
 */
function removeAnimationSEI() {
    rmPlot("animationSEI");
}

function generateAnimationBaseSEI(solution) {
    animate3D(solution, [0, 0, 0], ["Susceptible", "Exposed", "Infectious"], [0, 1, 2], "SEI", "Susceptible, exposed and recovered population phase plot for the SEI system.");
}
/**
 * Generate SEI animation
 * @return nothing
 */
function generateAnimationSEI() {
    var solution = solveProblem(RKF45, readInputs());
    generateAnimationBaseSEI(solution);
}

/**
 * Remove SER animation
 * @return nothing
 */
function removeAnimationSER() {
    rmPlot("animationSER");
}

function generateAnimationBaseSER(solution) {
    animate3D(solution, [0, 0, 0], ["Susceptible", "Exposed", "Recovered"], [0, 1, 3], "SER", "Susceptible, exposed and recovered population phase plot for the SEIR infectious disease model.");
}
/**
 * Generate SER animation
 * @return nothing
 */
function generateAnimationSER() {
    var solution = solveProblem(RKF45, readInputs());
    generateAnimationBaseSER(solution);
}

/**
 * Remove EIR animation
 * @return nothing
 */
function removeAnimationEIR() {
    rmPlot("animationEIR");
}

function generateAnimationBaseEIR(solution) {
    animate3D(solution, [0, 0, 0], ["Exposed", "Infectious", "Recovered"], [1, 2, 3], "EIR", "Exposed, infectious and recovered population phase plot for the SEIR infectious disease model.");
}

/**
 * Generate EIR animation
 * @return nothing
 */
function generateAnimationEIR() {
    var solution = solveProblem(RKF45, readInputs());
    generateAnimationBaseEIR(solution);
}

function generateAnimations() {
    var objectOfInputs = readInputs();
    var solution = solveProblem(RKF45, objectOfInputs);
    generateAnimationBaseSIR(solution);
    generateAnimationBaseSEI(solution);
    generateAnimationBaseSER(solution);
    generateAnimationBaseEIR(solution);
}