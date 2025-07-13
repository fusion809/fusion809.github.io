/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [theta1, p1, theta2, p2]
 * @return               [dtheta1/dt, dp1/dt, dtheta2/dt, dp2/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {g, l1, l2, m1b, m1r, m2b, m2r, b1b, b1r, c1b, c1r, b2b, b2r, c2b, c2r} = objectOfInputs;
    var [theta1, dtheta1, theta2, dtheta2] = vars;
    outercoef = 1/((m2r/12+m2b)*l2 - (m2b**2*l2*Math.cos(theta1-theta2)**2)/(m1r/12 + m1b+m2b));
    mass = m1r/2 + m2r + m1b + m2b; 
    v1b = l1*dtheta1;
    v1r = v1b/2;
    v2b = Math.sqrt(l1**2*dtheta1**2+l2**2*dtheta2**2+2*l1*l2*dtheta1*dtheta2*Math.cos(theta1-theta2));
    v2r = Math.sqrt(l1**2*dtheta1**2+(l2**2*dtheta2**2)/4 + l2*dtheta1*dtheta2*Math.cos(theta1-theta2));
    drag1b = (b1b + c1b*v1b)*v1b;
    drag1r = (b1r + c1r*v1r)*v1r/2;
    drag2b = (b2b + c2b*v2b)*(l1*dtheta1 + l2*dtheta2*Math.cos(theta1-theta2));
    drag2b2 = (b2b + c2b*v2b)*(l1*dtheta1*Math.cos(theta1-theta2) + l2*dtheta2);
    drag2r = (b2r + c2r*v2r)*(l1*dtheta1 + l2*dtheta2*Math.cos(theta1-theta2)/2);
    drag2l2 = 1/4*(b2r + c2r*v2r)*(2*l1*dtheta1*Math.cos(theta1-theta2) + l2*dtheta2);
    innercoef = -m2b*Math.cos(theta1-theta2)/(m1r/12+m1b+m2b);
    inner = innercoef*(-m2b*l2*dtheta2**2*Math.sin(theta1-theta2) - g*Math.cos(theta1)*(m1r/2+m2r+m1b+m2b)-drag1b-drag2b-drag1r-drag2r);
    extra = m2b*(l1*dtheta1**2*Math.sin(theta1-theta2)-g*Math.cos(theta2))-drag2l2-drag2b2;
    d2theta2 = outercoef*(inner+extra);
    outercoef1 = 1/((m1r/12+m1b+m2b)*l1);
    innel11 = -m2b*l2*(d2theta2*Math.cos(theta1-theta2)+dtheta2**2*Math.sin(theta1-theta2));
    innel12 = -g*Math.cos(theta1)*mass - drag1b - drag2b - drag1r - drag2r;
    d2theta1 = outercoef1*(innel11 + innel12);

    
    // Return statement
    return [dt*dtheta1, dt*d2theta1, dt*dtheta2, dt*d2theta2];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {theta10, dtheta10, theta20, dtheta20} = objectOfInputs;
    var vars0 = [[theta10, dtheta10, theta20, dtheta20]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0); 
    return [t, vars];
}

/**
 * Generates a 2D phase plot of theta2 against theta1
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta1Theta2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var theta2 = vars[2];

    // Generate 2D plot
    gen2DPlot(theta1, theta2, "phasePlotTheta1Theta2", "Phase plot of θ<sub>2</sub> against θ<sub>1</sub>", "θ<sub>1</sub>", "θ<sub>2</sub>");
}

/**
 * Generates a dtheta2 against dtheta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateDtheta1Dtheta2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var dtheta1 = vars[1];
    var dtheta2 = vars[3];

    // Generate 2D plot
    gen2DPlot(dtheta1, dtheta2, "phasePlotDtheta1Dtheta2", "Phase plot of dθ<sub>2</sub>/dt against dθ<sub>1</sub>/dt", "dθ<sub>1</sub>/dt", "dθ<sub>2</sub>/dt");
}

/**
 * Generates a dtheta1 against theta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta1Dtheta1PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var dtheta1 = vars[1];
    
    // Generate 2D plot
    gen2DPlot(theta1, dtheta1, "phasePlotTheta1Dtheta1", "Phase plot of dθ<sub>1</sub>/dt against θ<sub>1</sub>", "θ<sub>1</sub>", "dθ<sub>1</sub>/dt");
}

/**
 * Generates a dtheta2 against theta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta1Dtheta2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var dtheta2 = vars[2];
    
    // Generate 2D plot
    gen2DPlot(theta1, dtheta2, "phasePlotTheta1Dtheta2", "Phase plot of dθ<sub>2</sub>/dt against θ<sub>1</sub>", "θ<sub>1</sub>", "θ<sub>2</sub>/dt");
}

/**
 * Generates a dtheta1 against theta2 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta2Dtheta1PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta2 = vars[2];
    var dtheta1 = vars[1];
    
    // Generate 2D plot
    gen2DPlot(theta2, dtheta1, "phasePlotTheta2Dtheta1", "Phase plot of dθ<sub>1</sub>/dt against θ<sub>2</sub>", "θ<sub>2</sub>", "dθ<sub>1</sub>/dt");
}

/**
 * Generates a dtheta2 against theta2 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta2Dtheta2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta2 = vars[2];
    var dtheta2 = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta2, dtheta2, "phasePlotTheta2Dtheta2", "Phase plot of dθ<sub>2</sub>/dt against θ<sub>2</sub>", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], "timePlot", "Plot of $\\theta_1$, $\\dfrac{d\\theta_1}{dt}$, $\\theta_2$ and $\\dfrac{d\\theta_2}{dt}$ against time");
}

/**
 * Generate all plots
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generatePendulumPlots(objectOfInputs, solution);
    generateTheta1Theta2PhasePlot(solution);
    generateTheta1Dtheta1PhasePlot(solution);
    generateTheta1Dtheta2PhasePlot(solution);
    generateTheta2Dtheta1PhasePlot(solution);
    generateTheta2Dtheta2PhasePlot(solution);
    generateDtheta1Dtheta2PhasePlot(solution);
    generatePendulumPlots(objectOfInputs, solution)
    generateTimePlot(solution);
}

/**
 * Remove animation
 * @return nothing
 */
function removeAnimation() {
    rmPlot("animation");
}

/**
 * Remove Theta1 Dtheta1 phase plot
 * @return nothing
 */
function removeTheta1Dtheta1PhasePlot() {
  rmPlot("phasePlotTheta1Dtheta1");
}

/**
 * Remove Theta1 Dtheta2 phase plot
 * @return nothing
 */
function removeTheta1Dtheta2PhasePlot() {
  rmPlot("phasePlotTheta1Dtheta2");
}

/**
 * Remove Theta2 Dtheta1 phase plot
 * @return nothing
 */
function removeTheta2Dtheta1PhasePlot() {
  rmPlot("phasePlotTheta2Dtheta1");
}

/**
 * Remove Theta2 Dtheta2 phase plot
 * @return nothing
 */
function removeTheta2Dtheta2PhasePlot() {
  rmPlot("phasePlotTheta2Dtheta2");
}

/**
 * Remove Dtheta1 Dtheta2 phase plot
 * @return nothing
 */
function removeDtheta1Dtheta2PhasePlot() {
  rmPlot("phasePlotDtheta1Dtheta2");
}

function generateAnimationBase(objectOfInputs, solution) {
    animatePendulum(objectOfInputs, solution, "Double pendulum");
}

/**
 * Generate animation
 * @return nothing
 */
function generateAnimation() {
  var objectOfInputs = readInputs();
  var solution = solveProblem(RKF45, objectOfInputs);
  generateAnimationBase(objectOfInputs, solution);
}

function removeAnimation() {
    rmPlot("animation");
}

function generateTheta1PhaseAnimationBase(solution) {
    animate2D(solution, {varnames: ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt"], IdSuffix: "Theta1Phase", title: "Phase plot of $\\dfrac{d\\theta_1}{dt}$ against $\\theta_1$."});
}

function generateTheta1PhaseAnimation() {
  var objectOfInputs = readInputs();
  var solution = solveProblem(RKF45, objectOfInputs);
  generateTheta1PhaseAnimationBase(solution);
}

function removeTheta1PhaseAnimation() {
    rmPlot("animationTheta1Phase");
}

function generateTheta2PhaseAnimationBase(solution) {
    animate2D(solution, {varnames: ["θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], IdSuffix: "Theta2Phase", nos: [2, 3], title: "Phase plot of $\\dfrac{d\\theta_2}{dt}$ against $\\theta_2$."});
}

function generateTheta2PhaseAnimation() {
  var objectOfInputs = readInputs();
  var solution = solveProblem(RKF45, objectOfInputs);
  generateTheta2PhaseAnimationBase(solution);
}

function removeTheta2PhaseAnimation() {
    rmPlot("animationTheta2Phase");
}

function generateTable() {
    fillTable(readInputs(), ['&theta;<sub>1</sub>', 'd&theta;<sub>1</sub>/dt', '&theta;<sub>2</sub>', 'd&theta;<sub>2</sub>/dt'])
}

function generateAnimations() {
    var objectOfInputs = readInputs();
    var solution = solveProblem(RKF45, objectOfInputs);
    generateAnimationBase(objectOfInputs, solution)
    generateTheta1PhaseAnimationBase(solution);
    generateTheta2PhaseAnimationBase(solution);
}

function removeAnimations() {
    removeAnimation();
    removeTheta1PhaseAnimation();
    removeTheta2PhaseAnimation();
}