/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [x, xDot, theta, dtheta]
 * @return               [dx/dt, d2x/dt2, dtheta/dt, d2theta/dt2]
 */
function f(objectOfInputs, t, vars, dt) {
    var {g, l0, k, m} = objectOfInputs;
    var [z, dz, theta, dtheta] = vars;
    var zDDot = (l0+z)*dtheta**2 - k*z/m + g*Math.sin(theta);
    var thetaDDot = -g*Math.cos(theta)/(l0+z)-2*dz*dtheta/(l0+z);
    return [dt*dz, dt*zDDot, dt*dtheta, dt*thetaDDot];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {z0, dz0, theta0, dtheta0} = objectOfInputs;
    var vars0 = [[z0, dz0, theta0, dtheta0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0); 
    return [t, vars];
}

/**
 * Generates a 2D phase plot of x against theta
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateZThetaPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var z = vars[0];
    var theta = vars[2];

    // Generate 2D plot
    gen2DPlot(z, theta, "phasePlotZTheta", "Phase plot of $\\theta$ against $z$", "z", "θ");
}

/**
 * Generates a xdot against x phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateZZDotPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var z = vars[0];
    var dz = vars[1];

    // Generate 2D plot
    gen2DPlot(z, dz, "phasePlotZZDot", "Phase plot of $\\dfrac{dz}{dt}$ against $z$", "z", "dz/dt");
}

/**
 * Generates a theta dot against theta phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateThetaThetaDotPhasePlot(solution) {    
    // Extract solution data from solution object
    var {vars} = solution;
    var theta = vars[2];
    var dtheta = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta, dtheta, "phasePlotThetaThetaDot", "Phase plot of $\\dfrac{d\\theta}{dt}$ against $\\theta$", "θ", "dθ/dt");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["z", "dz/dt", "θ", "dθ/dt"], "timePlot", "Plot of $z$, $\\dfrac{dz}{dt}$, $\\theta$ and $\\dfrac{d\\theta}{dt}$ against time");
}

/**
 * Generate four plots:
 * - The first is a phase plot of theta against x.
 * - The second is a phase plot of x dot against x.
 * - The third is a phase plot of theta dot against theta.
 * - The fourth is a plot of x, x dot, theta and theta dot against time.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Plot solution
    generateZThetaPhasePlot(solution);
    generateZZDotPhasePlot(solution);
    generateThetaThetaDotPhasePlot(solution);
    generateTimePlot(solution);
    generatePendulumPlots(objectOfInputs, solution);
}

function removePendulumPlot() {
    rmPlot("pendulumPlot");
}
function removePendulumTimePlot() {
    rmPlot("pendulumTimePlot");
}
function removePendulumPlots() {
    removePendulumPlots();
    removePendulumTimePlot();
}

/**
 * Remove animation
 * @return nothing. 
 */
function removeAnimation() {
    rmPlot("animation");
}

function generateAnimationBase(objectOfInputs, solution) {
    animatePendulum(objectOfInputs, solution, "Elastic pendulum");
}

/**
 * Generate animation
 * @return nothing.
 */
function generateAnimation() {
    var objectOfInputs = readInputs();
    var solution = solveProblem(RKF45, objectOfInputs);
    generateAnimationBase(objectOfInputs, solution);
}

function generateThetaPhaseAnimationBase(solution) {
    animate2D(solution, {varnames: ["θ", "dθ/dt"], timer: [0.9, 0.98], IdSuffix: "ThetaPhase", nos: [2, 3], title: "Elastic pendulum: phase plot of $\\dfrac{d\\theta}{dt}$ against $\\theta$."});
}

function generateThetaPhaseAnimation() {
    var objectOfInputs = readInputs();
    var solution = solveProblem(RKF45, objectOfInputs);
    generateThetaPhaseAnimationBase(solution);
}

function removeThetaPhaseAnimation() {
    rmPlot("animationThetaPhase");
}

function generateZPhaseAnimationBase(solution) {
    animate2D(solution, {varnames: ["z", "dz/dt"], timer: [0.9, 0.98], IdSuffix: "ZPhase", title: "Elastic pendulum: phase plot of $\\dfrac{dz}{dt}$ against $z$."});
}

function generateZPhaseAnimation() {
    var objectOfInputs = readInputs();
    var solution = solveProblem(RKF45, objectOfInputs);
    generateZPhaseAnimationBase(solution);    
}

function generateAnimations() {
    var objectOfInputs = readInputs();
    var solution = solveProblem(RKF45, objectOfInputs);
    generateAnimationBase(objectOfInputs, solution);
    generateThetaPhaseAnimationBase(solution);
    generateZPhaseAnimationBase(solution); 
}

function removeZPhaseAnimation() {
    rmPlot("animationZPhase");
}

function removeThetaDthetaPhasePlot() {
    rmPlot("phasePlotThetaThetaDot");
}

function removeZThetaPhasePlot() {
    rmPlot("phasePlotZTheta");
}

function removeZZDotPhasePlot() {
    rmPlot("phasePlotZZDot");
}

function removeThetaPhaseAnimation() {
    rmPlot("animationThetaPhase");
}

function removeAnimations() {
    removeAnimation();
    removeThetaPhaseAnimation();
    removeZPhaseAnimation();
}