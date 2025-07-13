/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [theta1, p1, theta2, p2]
 * @return               [dtheta1/dt, dp1/dt, dtheta2/dt, dp2/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {g, l1, l2, m1, m2, k1, k2, b1, c1, b2, c2} = objectOfInputs;
    var [r1, dr1, r2, dr2, theta1, dtheta1, theta2, dtheta2] = vars;

    // Define Delta = theta2 - theta1
    var Delta = theta2 - theta1;
    var cosDelta = Math.cos(Delta);
    var sinDelta = Math.sin(Delta);
    var cosTheta1 = Math.cos(theta1);
    var cosTheta2 = Math.cos(theta2);
    var sinTheta1 = Math.sin(theta1);
    var sinTheta2 = Math.sin(theta2);
    var v1 = Math.sqrt(dr1**2 + r1**2 *  dtheta1**2);
    var v2 = Math.sqrt(v1**2 + dr2**2 + r2**2*dtheta2**2 + 2*cosDelta * (dr1*dr2 + r1*r2*dtheta1*dtheta2) + 2*sinDelta * (r1*dr2*dtheta1 - dr1*r2*dtheta2));
    var Qr1 = -(b1+c1*v1)*dr1 - (b2+c2*v2)*(dr1+dr2*cosDelta - r2*dtheta2*sinDelta);
    var Qr2 = -(b2+c2*v2)*(dr1*cosDelta + r1*dtheta1 * sinDelta + dr2);
    var Qtheta1 = -(b1+c1*v1)*r1**2 *dtheta1 - (b2+c2*v2)*(r1**2*dtheta1 + r1*dr2*sinDelta + r1*r2*dtheta2*cosDelta);
    var Qtheta2 = -(b2+c2*v2)*(r2**2*dtheta2 - dr1*r2*sinDelta + r1*r2*dtheta1*cosDelta);

// varruct the 4x4 matrix A
    var A = [
        [1, m2 * cosDelta / (m1 + m2), 0, -m2 * r2 * sinDelta / (m1 + m2)],
        [cosDelta, 1, r1 * sinDelta, 0],
        [0, m2 * sinDelta / ((m1 + m2) * r1), 1, m2 * r2 * cosDelta / ((m1 + m2) * r1)],
        [-sinDelta / r2, 0, r1 * cosDelta / r2, 1]
    ];

// varruct the RHS vector B
    var b = [
        r1 * dtheta1**2 - g * sinTheta1 + m2 / (m1 + m2) * (r2 * dtheta2**2 * cosDelta + 2 * dr2 * dtheta2 * sinDelta)
            + (Qr1 - k1 * (r1 - l1)) / (m1 + m2),

        r2 * dtheta2**2 - g * sinTheta2 + r1 * dtheta1**2 * cosDelta - 2 * dr1 * dtheta1 * sinDelta
            + (Qr2 - k2 * (r2 - l2)) / m2,

        -2 * dr1 * dtheta1 / r1 - g * cosTheta1 / r1 - (m2 / ((m1 + m2) * r1)) *
            (2 * dr2 * dtheta2 * cosDelta - r2 * dtheta2**2 * sinDelta)
            + Qtheta1 / ((m1 + m2) * r1**2),

        -2 * dr2 * dtheta2 / r2 - g * cosTheta2 / r2 - 2 * dr1 * dtheta1 * cosDelta / r2
            - r1 * dtheta1**2 * sinDelta / r2 + Qtheta2 / (m2 * r2**2)
    ];
    var d2 = math.lusolve(A, b);

    // Return statement
    return [dt*dr1, dt*d2[0], dt*dr2, dt*d2[1], dt*dtheta1, dt*d2[2], dt*dtheta2, dt*d2[3]];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {r10, dr10, r20, dr20, theta10, dtheta10, theta20, dtheta20} = objectOfInputs;
    var vars0 = [[r10, dr10, r20, dr20, theta10, dtheta10, theta20, dtheta20]];
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
    var theta1 = vars[4];
    var theta2 = vars[6];

    // Generate 2D plot
    gen2DPlot(theta1, theta2, "phasePlotTheta1Theta2", "Phase plot of θ<sub>2</sub> against θ<sub>1</sub>");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["r<sub>1</sub>", "dr<sub>1</sub>/dt", "r<sub>2</sub>", "dr<sub>2</sub>/dt", "θ<sub>1</sub>", "dθ<sub>1</sub>/dt", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], "timePlot", "Plot of r<sub>1</sub>, dr<sub>1</sub>/dt, r<sub>2</sub>, dr<sub>2</sub>/dt, θ<sub>1</sub>, dθ<sub>1</sub>/dt, θ<sub>2</sub> and dθ<sub>2</sub>/dt against time");
}

function generateR1TPlot(solution) {
    var {t, vars} = solution;
    var r1 = vars[0];
    gen2DPlotXYLabs(t, r1, "R1TPlot", "Time plot of r<sub>1</sub>", "t", "r<sub>1</sub>");
}

function generateDr1TPlot(solution) {
    var {t, vars} = solution;
    var dr1 = vars[1];
    gen2DPlotXYLabs(t, dr1, "Dr1TPlot", "Time plot of dr<sub>1</sub>/dt", "t", "dr<sub>1</sub>/dt");
}

function generateDr1R1Plot(solution) {
    var {vars} = solution;
    var r1 = vars[0];
    var dr1 = vars[1];
    gen2DPlotXYLabs(r1, dr1, "Dr1R1Plot", "Phase plot of dr<sub>1</sub>/dt vs r<sub>1</sub>", "r<sub>1</sub>", "dr<sub>1</sub>/dt");
}

function generateR2TPlot(solution) {
    var {t, vars} = solution;
    var r2 = vars[2];
    gen2DPlotXYLabs(t, r2, "R2TPlot", "Time plot of r<sub>2</sub>", "t", "r<sub>2</sub>");
}

function generateDr2TPlot(solution) {
    var {t, vars} = solution;
    var dr2 = vars[3];
    gen2DPlotXYLabs(t, dr2, "Dr2TPlot", "Time plot of dr<sub>2</sub>/dt", "t", "dr<sub>2</sub>/dt");
}

function generateDr2R2Plot(solution) {
    var {vars} = solution;
    var r2 = vars[2];
    var dr2 = vars[3];
    gen2DPlotXYLabs(r2, dr2, "Dr2R2Plot", "Phase plot of dr<sub>2</sub>/dt vs r<sub>2</sub>", "r<sub>2</sub>", "dr<sub>2</sub>/dt");
}

function generateTheta1TPlot(solution) {
    var {t, vars} = solution;
    var theta1 = vars[4];
    gen2DPlotXYLabs(t, theta1, "Theta1TPlot", "Time plot of θ<sub>1</sub>", "t", "θ<sub>1</sub>");
}

function generateDtheta1TPlot(solution) {
    var {t, vars} = solution;
    var dtheta1 = vars[5];
    gen2DPlotXYLabs(t, dtheta1, "Dtheta1TPlot", "Time plot of dθ<sub>1</sub>/dt", "t", "dθ<sub>1</sub>/dt");
}

function generateDtheta1Theta1Plot(solution) {
    var {vars} = solution;
    var theta1 = vars[4];
    var dtheta1 = vars[5];
    gen2DPlotXYLabs(theta1, dtheta1, "Dtheta1Theta1Plot", "Phase plot of dθ<sub>1</sub>/dt vs θ<sub>1</sub>", "θ<sub>1</sub>", "dθ<sub>1</sub>/dt");
}

function generateTheta2TPlot(solution) {
    var {t, vars} = solution;
    var theta2 = vars[6];
    gen2DPlotXYLabs(t, theta2, "Theta2TPlot", "Time plot of θ<sub>2</sub>", "t", "θ<sub>2</sub>");
}

function generateDtheta2TPlot(solution) {
    var {t, vars} = solution;
    var dtheta2 = vars[7];
    gen2DPlotXYLabs(t, dtheta2, "Dtheta2TPlot", "Time plot of dθ<sub>2</sub>/dt", "t", "dθ<sub>2</sub>/dt");
}

function generateDtheta2Theta2Plot(solution) {
    var {vars} = solution;
    var theta2 = vars[6];
    var dtheta2 = vars[7];
    gen2DPlotXYLabs(theta2, dtheta2, "Dtheta2Theta2Plot", "Phase plot of dθ<sub>2</sub>/dt vs θ<sub>2</sub>", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt");
}

function generateR1R2PhasePlot(solution) {
    var {vars} = solution;
    var r1 = vars[0];
    var r2 = vars[2];
    gen2DPlotXYLabs(r1, r2, "R2R1Plot", "Phase plot of r<sub>2</sub> vs r<sub>1</sub>", "r<sub>1</sub>", "r<sub>2</sub>");
}


function removePendulumPlots() {
    rmPlot("pendulumPlot");
    rmPlot("pendulumTimePlot");
}

function removeDr1TPlot() {
    rmPlot("Dr1TPlot");
}

function removeDr1R1Plot() {
    rmPlot("Dr1R1Plot");
}

function removeR2TPlot() {
    rmPlot("R2TPlot");
}

function removeDr2TPlot() {
    rmPlot("Dr2TPlot");
}

function removeDr2R2Plot() {
    rmPlot("Dr2R2Plot");
}

function removeTheta1TPlot() {
    rmPlot("Theta1TPlot");
}

function removeDtheta1TPlot() {
    rmPlot("Dtheta1TPlot");
}

function removeDtheta1Theta1Plot() {
    rmPlot("Dtheta1Theta1Plot");
}

function removeTheta2TPlot() {
    rmPlot("Theta2TPlot");
}

function removeDtheta2TPlot() {
    rmPlot("Dtheta2TPlot");
}

function removeDtheta2Theta2Plot() {
    rmPlot("Dtheta2Theta2Plot");
}

function removeR1R2PhasePlot() {
    rmPlot("R2R1Plot");
}

function removeTimePlot() {
    rmPlot("timePlot");
}

function rmPlots() {
    removePendulumPlots();
    removeTimePlot();
    removeTheta1Theta2PhasePlot();
    removeR1TPlot();
    removeDr1TPlot();
    removeDr1R1Plot()
    removeR2TPlot()
    removeDr2TPlot()
    removeDr2R2Plot()
    removeTheta1TPlot()
    removeDtheta1TPlot()
    removeDtheta1Theta1Plot()
    removeTheta2TPlot()
    removeDtheta2TPlot()
    removeDtheta2Theta2Plot()
    removeR1R2PhasePlot();
}

/**
 * Generate cartesian coordinates 
 * @param func           Function being used to integrate problem
 * @param objectOfInputs Problem parameters.
 */
function generatePendulumCoords(solution) {
    // Extract solution values and pendulum lengths
    var {t, vars} = solution;
    var [r1, dr1, r2, dr2, theta1, dtheta1, theta2, dtheta2] = vars;
    var N = theta1.length;

    // Initialize arrays that will store x and y coords
    var x1 = new Array(N);
    var x2 = new Array(N);
    var y1 = new Array(N);
    var y2 = new Array(N);
    for (let i = 0; i < N; i++) {
        x1[i] = r1[i]*Math.cos(theta1[i]);
        y1[i] = r1[i]*Math.sin(theta1[i]);
        x2[i] = x1[i] + r2[i]*Math.cos(theta2[i]);
        y2[i] = y1[i] + r2[i]*Math.sin(theta2[i]);
    }

    // Return t and Cartesian coordinates of the pendulum bobs
    return [t, x1, y1, x2, y2];
}

/**
 * Generates two plots pertaining to the location of the bobs
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generatePendulumPlots(solution) {
    var [t, x1, y1, x2, y2] = generatePendulumCoords(solution);
    adjustPlotHeight("pendulumPlot");
    adjustPlotHeight("pendulumTimePlot");
    
    // Show two pendulum bob locations on the same plot
    var plotPen1 = {
        x: x1,
        y: y1,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: "Pendulum 1 bob"
    }
    var plotPen2 = {
        x: x2,
        y: y2,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: "Pendulum 2 bob"
    }
    var dataPen = [plotPen1, plotPen2];
    var layoutPen = {
        title: "Pendulum coordinate plots"
    };
    Plotly.newPlot("pendulumPlot", dataPen, layoutPen);
    
    // Plot pendulum bob location against time plot
    var plotPen1Time = {
        x: t,
        y: x1,
        z: y1,
        type: 'scatter3d',
        mode: 'lines',
        opacity: 1,
        line: {
            width: 6,
            reversescale: false
        },
        name: "Pendulum 1 bob"
    };
    var plotPen2Time = {
        x: t,
        y: x2,
        z: y2,
        type: 'scatter3d',
        mode: 'lines',
        opacity: 1,
        line: {
            width: 6,
            reversescale: false
        },
        name: 'Pendulum 2 bob'
    };
    var dataPen = [plotPen1Time, plotPen2Time];
    var layoutPenTime = {
        title: "Pendulum bob position against time plot",
        xaxis: {
            title: {text: "x"}
        },
        yaxis: {
            title: {text: "y"}
        }
    };
    Plotly.newPlot("pendulumTimePlot", dataPen, layoutPenTime);
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
    generatePendulumPlots(solution);
    generateTimePlot(solution);
    generateTheta1Theta2PhasePlot(solution);
    generateR1R2PhasePlot(solution);
    generateR1TPlot(solution);
    generateDr1TPlot(solution);
    generateDr1R1Plot(solution);
    generateR2TPlot(solution);
    generateDr2TPlot(solution);
    generateDr2R2Plot(solution);
    generateTheta1TPlot(solution);
    generateDtheta1TPlot(solution);
    generateDtheta1Theta1Plot(solution);
    generateTheta2TPlot(solution);
    generateDtheta2TPlot(solution);
    generateDtheta2Theta2Plot(solution);
}

function generateAnimationBase(objectOfInputs, solution) {
    animatePendulum(objectOfInputs, solution, "Double elastic pendulum");
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

/**
 * Remove animation
 * @return nothing
 */
function removeAnimation() {
    rmPlot("animation");
}

function generateR1PhaseAnimationBase(solution) {
    animate2D(solution, {varnames: ["r<sub>1</sub>", "dr<sub>1</sub>/dt"],IdSuffix: "R1Phase", title: "Phase plot of $\\dfrac{dr_1}{dt}$ against $r_1$."});
}

function generateR1PhaseAnimation() {
  var objectOfInputs = readInputs();
  var solution = solveProblem(RKF45, objectOfInputs);
  generateR1PhaseAnimationBase(solution);
}

function removeR1PhaseAnimation() {
    rmPlot("animationR1Phase");
}

function generateR2PhaseAnimationBase(solution) {
    animate2D(solution, {varnames: ["r<sub>2</sub>", "dr<sub>2</sub>/dt"], IdSuffix: "R2Phase", nos: [2, 3], title: "Phase plot of $\\dfrac{dr_2}{dt}$ against $r_2$."});
}

function generateR2PhaseAnimation() {
  var objectOfInputs = readInputs();
  var solution = solveProblem(RKF45, objectOfInputs);
  generateR2PhaseAnimationBase(solution);
}

function removeR2PhaseAnimation() {
    rmPlot("animationR2Phase");
}

function generateTheta1PhaseAnimationBase(solution) {
    animate2D(solution, {varnames: ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt"],IdSuffix: "Theta1Phase", nos: [4, 5], title: "Phase plot of $\\dfrac{d\\theta_1}{dt}$ against $\\theta_1$."});
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
    animate2D(solution, {varnames: ["θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], timer: [0.98, 0.02], IdSuffix: "Theta2Phase", nos: [6, 7], title: "Phase plot of $\\dfrac{d\\theta_2}{dt}$ against $\\theta_2$."});
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
    fillTable(readInputs(), ['r<sub>1</sub>', 'dr<sub>1</sub>/dt', 'r<sub>2</sub>', 'dr<sub>2</sub>/dt', '&theta;<sub>1</sub>', 'd&theta;<sub>1</sub>/dt', '&theta;<sub>2</sub>', 'd&theta;<sub>2</sub>/dt'])
}

function generateAnimations() {
    var objectOfInputs = readInputs();
    var solution = solveProblem(RKF45, objectOfInputs);
    generateAnimationBase(objectOfInputs, solution)
    generateR1PhaseAnimationBase(solution);
    generateR2PhaseAnimationBase(solution);
    generateTheta1PhaseAnimationBase(solution);
    generateTheta2PhaseAnimationBase(solution);
}

function removeAnimations() {
    removeAnimation();
    removeR1PhaseAnimation();
    removeR2PhaseAnimation();
    removeTheta1PhaseAnimation();
    removeTheta2PhaseAnimation();
}