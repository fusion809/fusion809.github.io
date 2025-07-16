/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [theta1, p1, theta2, p2]
 * @return               [dtheta1/dt, dp1/dt, dtheta2/dt, dp2/dt]
 */
// === Solve A * ddtheta = b ===
// Use simple Gauss-Jordan elimination or a numeric library (like math.js) if preferred
function solveLinearSystem(A, b) {
  const math = window.math || require("mathjs");
  return math.lusolve(math.matrix(A), math.matrix(b)).toArray().flat();
}

// === Dissipative and driving terms Qtheta ===
function Qtheta1(objectOfInputs, vars) {
    const {
    l2, l3,
    b2b, b2r, b3b, b3r,
    c2b, c2r, c3b, c3r
  } = objectOfInputs;

  const [theta1, dtheta1, theta2, dtheta2, theta3, dtheta3] = vars;
    const v1 = dtheta1;
    const v2 = dtheta2;
    const v3 = dtheta3;
    const d12 = theta1 - theta2;
  const d13 = theta1 - theta3;
  const d23 = theta2 - theta3;
    const s = Math.sin, c = Math.cos, sqrt = Math.sqrt;
    const cos12 = c(d12), cos13 = c(d13), cos23 = c(d23);

    const sqrt1 = sqrt(l1 ** 2 * v1 ** 2);
    const sqrt2 = sqrt(l1 ** 2 * v1 ** 2 + 2 * l1 * l2 * cos12 * v1 * v2 + l2 ** 2 * v2 ** 2);
    const sqrt2r = sqrt(4 * l1 ** 2 * v1 ** 2 + 4 * l1 * l2 * cos12 * v1 * v2 + l2 ** 2 * v2 ** 2);

    const sqrt3 = sqrt(
      l1 ** 2 * v1 ** 2 + 2 * l1 * l2 * cos12 * v1 * v2 + 2 * l1 * l3 * cos13 * v1 * v3
      + l2 ** 2 * v2 ** 2 + 2 * l2 * l3 * cos23 * v2 * v3 + l3 ** 2 * v3 ** 2
    );
    const sqrt3r = sqrt(
      4 * l1 ** 2 * v1 ** 2 + 8 * l1 * l2 * cos12 * v1 * v2 + 4 * l1 * l3 * cos13 * v1 * v3
      + 4 * l2 ** 2 * v2 ** 2 + 4 * l2 * l3 * cos23 * v2 * v3 + l3 ** 2 * v3 ** 2
    );

    const term1 = l1 * (
      8 * b1b * v1 + 2 * b1r * v1 + 8 * b2b * v1 + 8 * b2r * v1 + 8 * b3b * v1 + 8 * b3r * v1
      + 8 * c1b * sqrt1 * v1 + c1r * sqrt1 * v1
      + 8 * c2b * sqrt2 * v1 + 4 * c2r * sqrt2r * v1
      + 8 * c3b * sqrt3 * v1 + 4 * c3r * sqrt3r * v1
    );

    const term2 = l2 * (
      8 * b2b * cos12 * v2 + 4 * b2r * cos12 * v2 + 8 * b3b * cos12 * v2 + 8 * b3r * cos12 * v2
      + 8 * c2b * sqrt2 * cos12 * v2 + 2 * c2r * sqrt2r * cos12 * v2
      + 8 * c3b * sqrt3 * cos12 * v2 + 4 * c3r * sqrt3r * cos12 * v2
    );

    const term3 = l3 * (
      8 * b3b * cos13 * v3 + 4 * b3r * cos13 * v3
      + 8 * c3b * sqrt3 * cos13 * v3 + 2 * c3r * sqrt3r * cos13 * v3
    );

    return -(term1 + term2 + term3);
}

function Qtheta2(objectOfInputs, vars) {
  const {
    l2, l3,
    b2b, b2r, b3b, b3r,
    c2b, c2r, c3b, c3r
  } = objectOfInputs;

  const [theta1, dtheta1, theta2, dtheta2, theta3, dtheta3] = vars;

  const Delta21 = theta2 - theta1;
  const Delta32 = theta3 - theta2;

  const v2bMag = Math.abs(l2 * dtheta2) / Math.SQRT2;
  const v2rMag = Math.abs(l2 * dtheta2) / (2 * Math.SQRT2);
  const v3bMag = Math.abs(l3 * dtheta3) / Math.SQRT2;
  const v3rMag = Math.abs(l3 * dtheta3) / (2 * Math.SQRT2);

  const Qtheta2 =
    - (b2b + c2b * v2bMag) * l2 * l2 * dtheta2
    - (b2r + c2r * v2rMag) * l2 * l2 * dtheta2 / 4
    - (b3b + c3b * v3bMag) * l2 * l3 * dtheta3 * Math.cos(Delta32)
    - (b3r + c3r * v3rMag) * l2 * l3 * dtheta3 * Math.cos(Delta32) / 2;

  return Qtheta2;
}

function Qtheta3(objectOfInputs, vars) {
  const {
    l2, l3,
    b3b, b3r,
    c3b, c3r
  } = objectOfInputs;

  const [theta1, dtheta1, theta2, dtheta2, theta3, dtheta3] = vars;

  const Delta31 = theta3 - theta1;
  const Delta32 = theta3 - theta2;

  const v3bMag = Math.abs(l3 * dtheta3) / Math.SQRT2;
  const v3rMag = Math.abs(l3 * dtheta3) / (2 * Math.SQRT2);
  const v2bMag = Math.abs(l2 * dtheta2) / Math.SQRT2;
  const v2rMag = Math.abs(l2 * dtheta2) / (2 * Math.SQRT2);

  const Qtheta3 =
    - (b3b + c3b * v3bMag) * l3 * l3 * dtheta3
    - (b3r + c3r * v3rMag) * l3 * l3 * dtheta3 / 4
    - (b3b + c3b * v2bMag) * l2 * l3 * dtheta2 * Math.cos(Delta32)
    - (b3r + c3r * v2rMag) * l2 * l3 * dtheta2 * Math.cos(Delta32) / 2;

  return Qtheta3;
}


function computeDdTheta(objectOfInputs, vars) {
  const {
    l1, l2, l3, g,
    m1b, m1r, m2b, m2r, m3b, m3r,
    b1b, b1r, b2b, b2r, b3b, b3r,
    c1b, c1r, c2b, c2r, c3b, c3r
  } = objectOfInputs;

  const [theta1, dtheta1, theta2, dtheta2, theta3, dtheta3] = vars;
  
  const s = Math.sin, c = Math.cos, sqrt = Math.sqrt;

  const d12 = theta1 - theta2;
  const d13 = theta1 - theta3;
  const d23 = theta2 - theta3;

  // === Inertia Matrix A ===
  const A = [
    [
      l1 ** 2 * (m1r / 12 + 0.5 * (m1b + 0.5 * m1r + m2b + m2r + m3b + m3r)),
      l1 * l2 * c(d12) * (m2b + 0.5 * m2r + m3b + m3r),
      l1 * l3 * c(d13) * (m3b + 0.5 * m3r),
    ],
    [
      l1 * l2 * c(d12) * (m2b + 0.5 * m2r + m3b + m3r),
      l2 ** 2 * (m2r / 12 + 0.5 * (m2b + 0.5 * m2r + m3b + m3r)),
      l2 * l3 * c(d23) * (m3b + 0.5 * m3r),
    ],
    [
      l1 * l3 * c(d13) * (m3b + 0.5 * m3r),
      l2 * l3 * c(d23) * (m3b + 0.5 * m3r),
      l3 ** 2 * (m3r / 12 + 0.5 * (m3b + 0.5 * m3r))
    ]
  ];

  // === Gravitational and centrifugal force vector b ===
  const b = [
    -g * l1 * c(theta1) * (m1b + 0.5 * m1r + m2b + m2r + m3b + m3r)
    - l1 * l2 * s(d12) * dtheta2 ** 2 * (m2b + 0.5 * m2r + m3b + m3r)
    - l1 * l3 * s(d13) * dtheta3 ** 2 * (m3b + 0.5 * m3r),

    l1 * l2 * s(d12) * dtheta1 ** 2 * (m2b + 0.5 * m2r + m3b + m3r)
    - g * l2 * c(theta2) * (m2b + 0.5 * m2r + m3b + m3r)
    - l2 * l3 * s(d23) * dtheta3 ** 2 * (m3b + 0.5 * m3r),

    l1 * l3 * s(d13) * dtheta1 ** 2 * (m3b + 0.5 * m3r)
    - g * l3 * c(theta3) * (m3b + 0.5 * m3r)
    + l2 * l3 * s(d23) * dtheta2 ** 2 * (m3b + 0.5 * m3r),
  ];

  b[0] += Qtheta1(objectOfInputs, vars);
  b[1] += Qtheta2(objectOfInputs, vars);
  b[2] += Qtheta3(objectOfInputs, vars);

  // Similar functions Qtheta2 and Qtheta3 omitted for brevity
  // b[1] += Qtheta2();
  // b[2] += Qtheta3();

  // === Solve A * ddtheta = b ===
  const ddtheta = math.lusolve(math.matrix(A), math.matrix(b)).toArray().flat();

  return ddtheta; // [ddtheta1, ddtheta2, ddtheta3]
}

function f(objectOfInputs, t, vars, dt) {
    var [theta1, dtheta1, theta2, dtheta2, theta3, dtheta3] = vars;
    var d2 = computeDdTheta(objectOfInputs, vars)
    
    // Return statement
    return [dt*dtheta1, dt*d2[0], dt*dtheta2, dt*d2[1], dt*dtheta3, dt*d2[2]];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {theta10, dtheta10, theta20, dtheta20, theta30, dtheta30} = objectOfInputs;
    var vars0 = [[theta10, dtheta10, theta20, dtheta20, theta30, dtheta30]];
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
    gen2DPlot(dtheta1, dtheta2, "phasePlotDtheta2Dtheta1", "Phase plot of dθ<sub>2</sub>/dt against dθ<sub>1</sub>/dt", "dθ<sub>1</sub>/dt", "dθ<sub>2</sub>/dt");
}

/**
 * Generates a dtheta1 against theta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateDtheta1Theta1PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var dtheta1 = vars[1];
    
    // Generate 2D plot
    gen2DPlot(theta1, dtheta1, "phasePlotDtheta1Theta1", "Phase plot of dθ<sub>1</sub>/dt against θ<sub>1</sub>", "θ<sub>1</sub>", "dθ<sub>1</sub>/dt");
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
    gen2DPlot(theta1, dtheta2, "phasePlotDtheta1Theta2", "Phase plot of dθ<sub>2</sub>/dt against θ<sub>1</sub>", "θ<sub>1</sub>", "θ<sub>2</sub>/dt");
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
    gen2DPlot(theta2, dtheta1, "phasePlotDtheta2Theta1", "Phase plot of dθ<sub>1</sub>/dt against θ<sub>2</sub>", "θ<sub>2</sub>", "dθ<sub>1</sub>/dt");
}

/**
 * Generates a dtheta2 against theta2 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateDtheta2Theta2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta2 = vars[2];
    var dtheta2 = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta2, dtheta2, "phasePlotDtheta2Theta2", "Phase plot of dθ<sub>2</sub>/dt against θ<sub>2</sub>", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], "timePlot", "Plot of θ<sub>1</sub>, dθ<sub>1</sub>/dt, θ<sub>2</sub> and dθ<sub>2</sub>/dt against time");
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
    generateDtheta1Theta1PhasePlot(solution);
    generateTheta1Dtheta2PhasePlot(solution);
    generateTheta2Dtheta1PhasePlot(solution);
    generateDtheta2Theta2PhasePlot(solution);
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
function removeDtheta1Theta1PhasePlot() {
  rmPlot("phasePlotDtheta1Theta1");
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
function removeDtheta2Theta2PhasePlot() {
  rmPlot("phasePlotDtheta2Theta2");
}

/**
 * Remove Dtheta1 Dtheta2 phase plot
 * @return nothing
 */
function removeDtheta1Dtheta2PhasePlot() {
  rmPlot("phasePlotDtheta1Dtheta2");
}

function generateAnimationBase(objectOfInputs, solution) {
    animatePendulum(objectOfInputs, solution, "Triple pendulum");
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
    animate2D(solution, {varnames: ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt"], IdSuffix: "Theta1Phase", title: "Phase plot of dθ<sub>1</sub>/dt against θ<sub>1</sub>."});
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
    animate2D(solution, {timer: [0, 1.0], varnames: ["θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], IdSuffix: "Theta2Phase", nos: [2, 3], title: "Phase plot of dθ<sub>2</sub>/dt against θ<sub>2</sub>."});
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
    fillTable(readInputs(), ['&theta;<sub>1</sub>', 'd&theta;<sub>1</sub>/dt', '&theta;<sub>2</sub>', 'd&theta;<sub>2</sub>/dt', '&theta;<sub>3</sub>', 'd&theta;<sub>3</sub>/dt'])
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