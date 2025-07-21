/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [x, xDot, theta, dtheta]
 * @return               [dx/dt, d2x/dt2, dtheta/dt, d2theta/dt2]
 */
function f(objectOfInputs, t, vars, dt) {
    var {g, l0, k, mb, mr, bb, br, cb, cr} = objectOfInputs;
    var [r, dr, theta, dtheta] = vars;
    var M1 = mb + mr/3;
    var mu1 = mb + mr/2;
    var vb = Math.sqrt(dr**2+r**2*dtheta**2)
    var vr = 1/2*vb;
    var Qr = -(bb+cb*vb)*dr - (br+cr*vr)*dr/4;
    var Qtheta = -(bb+cb*vb)*r**2*dtheta - (br+cr*vr)*r**2*dtheta/4;
    var d2r = r*dtheta**2 - mu1/M1 *g*Math.sin(theta) - k*(r-l0)/M1 + Qr/M1;
    var d2theta = -2*dr*dtheta/r - mu1*g*Math.cos(theta)/(M1*r) + Qtheta/(M1*r**2);
    return [dt*dr, dt*d2r, dt*dtheta, dt*d2theta];
}

/**
 * Generates a 2D phase plot of r against theta
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateRThetaPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var r = vars[0];
    var theta = vars[2];

    // Generate 2D plot
    gen2DPlot(r, theta, "phasePlotRTheta", "Phase plot of θ against r", "r", "θ");
}

/**
 * Generates a rdot against r phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateDrRPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var r = vars[0];
    var dr = vars[1];

    // Generate 2D plot
    gen2DPlot(r, dr, "phasePlotDrR", "Phase plot of dr/dt against r", "r", "dr/dt");
}

/**
 * Generates a theta dot against theta phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateDThetaThetaPhasePlot(solution) {    
    // Extract solution data from solution object
    var {vars} = solution;
    var theta = vars[2];
    var dtheta = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta, dtheta, "phasePlotDthetaTheta", "Phase plot of dθ/dt against θ", "θ", "dθ/dt");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["r", "dr/dt", "θ", "dθ/dt"], "timePlot", "Plot of r, dr/dt, θ and dθ/dt against time");
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
    generateRThetaPhasePlot(solution);
    generateDrRPhasePlot(solution);
    generateDThetaThetaPhasePlot(solution);
    generateTimePlot(solution);
    generatePendulumPlots(objectOfInputs, solution);
}

/**
 * Generate animation
 * @return nothing.
 */
function generateAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animatePendulum(objectOfInputs, solution, "Elastic pendulum");
}

function generateThetaPhaseAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["θ", "dθ/dt"], timer: [0.9, 0.98], IdSuffix: "ThetaPhase", nos: [2, 3], title: "Elastic pendulum: phase plot of dθ/dt against θ."});
}

function generateRPhaseAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["r", "dr/dt"], timer: [0.9, 0.98], IdSuffix: "RPhase", title: "Elastic pendulum: phase plot of dr/dt against r."});  
}

function generateAnimations(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    generateAnimation(objectOfInputs, solution);
    generateThetaPhaseAnimation(objectOfInputs, solution);
    generateRPhaseAnimation(objectOfInputs, solution); 
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    fillTable(objectOfInputs, ['r', 'dr/dt', '&theta;', 'd&theta;/dt'], solution)
}