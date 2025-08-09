/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [theta1, p1, theta2, p2]
 * @return               [dtheta1/dt, dp1/dt, dtheta2/dt, dp2/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    //var {g, l1, l2, m1, m2, k1, k2, b1, c1, b2, c2} = objectOfInputs;
    var {g, l1, l2, m1b, m1r, m2b, m2r, k1, k2, b1b, c1b, b1r, c1r, b2b, b2r, c2b, c2r} = objectOfInputs;
    var M1 = m1b + m1r/3 + m2b + m2r;
    var M2 = m2b + m2r/3;
    var mu1 = m1b + m1r/2 + m2b + m2r;
    var mu2 = m2b + m2r/2;
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
    // var Qr1 = -(b1+c1*v1)*dr1 - (b2+c2*v2)*(dr1+dr2*cosDelta - r2*dtheta2*sinDelta);
    // var Qr2 = -(b2+c2*v2)*(dr1*cosDelta + r1*dtheta1 * sinDelta + dr2);
    // var Qtheta1 = -(b1+c1*v1)*r1**2 *dtheta1 - (b2+c2*v2)*(r1**2*dtheta1 + r1*dr2*sinDelta + r1*r2*dtheta2*cosDelta);
    // var Qtheta2 = -(b2+c2*v2)*(r2**2*dtheta2 - dr1*r2*sinDelta + r1*r2*dtheta1*cosDelta);
    var v1b = v1;
    var v2b = v2;
    var v1r = v1b/2;
    var v2r = Math.sqrt(r1**2*dtheta1**2 + r1*r2*Math.cos(Delta)*dtheta1*dtheta2 - r1*-Math.sin(Delta)*dr2*dtheta1 + r2**2*dtheta2**2/4 + r2*-Math.sin(Delta)*dr1*dtheta2 + Math.cos(Delta)*dr1*dr2 + dr1**2 + dr2**2/4)

// varruct the 4x4 matrix A
    // var A = [
    //     [1, m2 * cosDelta / (m1 + m2), 0, -m2 * r2 * sinDelta / (m1 + m2)],
    //     [cosDelta, 1, r1 * sinDelta, 0],
    //     [0, m2 * sinDelta / ((m1 + m2) * r1), 1, m2 * r2 * cosDelta / ((m1 + m2) * r1)],
    //     [-sinDelta / r2, 0, r1 * cosDelta / r2, 1]
    // ];
// SymPy
    // var A = [
    //     [M1, 2*mu2*Math.cos(Delta), 0, -2*mu2*r2*Math.sin(Delta)],
    //     [2*mu2*Math.cos(Delta), M2, 2*mu2*r1*Math.sin(Delta), 0],
    //     [0, 2*mu2*r1*Math.sin(Delta), M1*r1**2, 2*mu2*r1*r2*Math.cos(Delta)],
    //     [-2*mu2*r2*Math.sin(Delta), 0, 2*mu2*r1*r2*Math.cos(Delta), M2*r2**2]
    // ]
    var A = [
        [M1, mu2*Math.cos(Delta), 0, -mu2*r2*Math.sin(Delta)],
        [mu2*Math.cos(Delta), M2, mu2*r1*Math.sin(Delta), 0],
        [0, mu2*r1*Math.sin(Delta), M1*r1**2, mu2*r1*r2*Math.cos(Delta)],
        [-mu2*r2*Math.sin(Delta), 0, mu2*r1*r2*Math.cos(Delta), M2*r2**2]
    ]
    // var A = [
    //     [M1, mu2*Math.cos(Delta), mu2*r1*r2*Math.cos(Delta), -mu2*r2*Math.sin(Delta)],
    //     [mu2*Math.cos(Delta), M2, -mu2*r1*Math.sin(Delta), mu2*r1*r2*Math.cos(Delta)],
    //     [mu2*r1*r2*Math.cos(Delta), -mu2*r1*Math.sin(Delta), M1*r1**2, 0],
    //     [-mu2*r2*Math.sin(Delta), mu2*r1*r2*Math.cos(Delta), 0, M2*r2**2]
    // ]

// varruct the RHS vector B
    // var b = [
    //     r1 * dtheta1**2 - g * sinTheta1 + m2 / (m1 + m2) * (r2 * dtheta2**2 * cosDelta + 2 * dr2 * dtheta2 * sinDelta)
    //         + (Qr1 - k1 * (r1 - l1)) / (m1 + m2),

    //     r2 * dtheta2**2 - g * sinTheta2 + r1 * dtheta1**2 * cosDelta - 2 * dr1 * dtheta1 * sinDelta
    //         + (Qr2 - k2 * (r2 - l2)) / m2,

    //     -2 * dr1 * dtheta1 / r1 - g * cosTheta1 / r1 - (m2 / ((m1 + m2) * r1)) *
    //         (2 * dr2 * dtheta2 * cosDelta - r2 * dtheta2**2 * sinDelta)
    //         + Qtheta1 / ((m1 + m2) * r1**2),

    //     -2 * dr2 * dtheta2 / r2 - g * cosTheta2 / r2 - 2 * dr1 * dtheta1 * cosDelta / r2
    //         - r1 * dtheta1**2 * sinDelta / r2 + Qtheta2 / (m2 * r2**2)
    // ];
    // From SymPy
    // var b = [
    //     M1*r1*dtheta1**2 - k1*(r1-l1) + 2*mu2*r2*dtheta2**2*Math.cos(Delta) - mu1*g*Math.sin(theta1)+4*mu2*dr2*dtheta2*Math.sin(Delta) + Q[0],
    //     M2*r2*dtheta2**2 - k2*(r2-l2) + 2*mu2*r1*dtheta1**2*Math.cos(Delta) - mu2*g*Math.sin(theta2)-4*mu2*dr1*dtheta1*Math.sin(Delta) + Q[1],
    //     -2*M1*r1*dr1*dtheta1 - mu1*g*r1*Math.cos(theta1) + 2*mu2*r1*r2*dtheta2**2*Math.sin(Delta) - 4*mu2*r1*dr2*dtheta2*Math.cos(Delta) + Q[2],
    //     -2*M2*r2*dr2*dtheta2 - mu2*g*r2*Math.cos(theta2) - 2*mu2*r1*r2*dtheta1**2*Math.sin(Delta) - 4*mu2*dr1*r2*dtheta1*Math.cos(Delta)+Q[3]
    // ]
    var Qr1 = -(b1b+c1b*v1b)*dr1 - dr1/4*(b1r+c1r*v1r) - (b2b+c2b*v2b)*(dr1+dr2*Math.cos(Delta)-r2*dtheta2*Math.sin(Delta)) - (b2r+c2r*v2r)*(dr1+(dr2*Math.cos(Delta)-r2*dtheta2*Math.sin(Delta))/2);
    var Qr2 = -(b2b+c2b*v2b)*(dr1*Math.cos(Delta) + r1*dtheta1*Math.sin(Delta)+dr2) -1/2*(b2r+c2r*v2r)*(dr1*Math.cos(Delta) + r1*dtheta1*Math.sin(Delta)+dr2/2)
    var Qtheta1 = -(b1b+c1b*v1b)*r1**2*dtheta1 - 1/4*(b1r+c1r*v1r)*(r1**2*dtheta1) - (b2b+c2b*v2b)*(r1**2*dtheta1+r1*dr2*Math.sin(Delta)+r1*r2*dtheta2*Math.cos(Delta)) - (b2r+c2r*v2r)*(r1**2*dtheta1+(r1*dr2*Math.sin(Delta)+r1*r2*dtheta2*Math.cos(Delta))/2);
    var Qtheta2 = -(b2b+c2b*v2b)*(r2**2*dtheta2-dr1*r2*Math.sin(Delta)+r1*r2*dtheta1*Math.cos(Delta)) - 1/2*(b2r+c2r*v2r)*(r2**2*dtheta2/2 - dr1*r2*Math.sin(Delta) + r1*r2*dtheta1*Math.cos(Delta))
    var b = [
        M1*r1*dtheta1**2 - g*mu1*Math.sin(theta1) + k1*(l1 - r1) + mu2*(r2*Math.cos(Delta)*dtheta2**2 + 2*Math.sin(Delta)*dr2*dtheta2) + Qr1,
        M2*r2*dtheta2**2 + k2*(l2 - r2) + mu2*(-g*Math.sin(theta2) + r1*Math.cos(Delta)*dtheta1**2 - 2*Math.sin(Delta)*dr1*dtheta1) + Qr2,
        -2*M1*r1*dr1*dtheta1 - g*mu1*r1*Math.cos(theta1) + mu2*(r1*r2*Math.sin(Delta)*dtheta2**2 - 2*r1*Math.cos(Delta)*dr2*dtheta2) + Qtheta1,
        -2*M2*r2*dr2*dtheta2 + mu2*(-g*r2*Math.cos(theta2) - r1*r2*Math.sin(Delta)*dtheta1**2 - 2*r2*Math.cos(Delta)*dr1*dtheta1) + Qtheta2
    ]
    // var b = [
    //     mu2*(Math.cos(Delta)*(dtheta1*dtheta2*r2+dr2*dtheta1*dtheta2)+Math.sin(Delta)*(dr2*dtheta1-r2*dtheta1*dtheta2))-mu1*g*Math.sin(theta1)-k1*(r1-l1) + Qr1,
    //     mu2*(Math.cos(Delta)*(-dr1*dtheta1*dtheta2+r1*dtheta1**2)-Math.sin(Delta)*(dr1*dtheta2))
    // ]
    var d2 = math.lusolve(A, b);

    // Return statement
    return [dt*dr1, dt*d2[0], dt*dr2, dt*d2[1], dt*dtheta1, dt*d2[2], dt*dtheta2, dt*d2[3]];
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
    gen2DPlot(theta1, theta2, "plotTheta1Theta2", "Phase plot of θ<sub>2</sub> against θ<sub>1</sub>.", "θ<sub>1</sub>", "θ<sub>2</sub>");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["r<sub>1</sub>", "dr<sub>1</sub>/dt", "r<sub>2</sub>", "dr<sub>2</sub>/dt", "θ<sub>1</sub>", "dθ<sub>1</sub>/dt", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], "timePlot", "Plot of r<sub>1</sub>, dr<sub>1</sub>/dt, r<sub>2</sub>, dr<sub>2</sub>/dt, θ<sub>1</sub>, dθ<sub>1</sub>/dt, θ<sub>2</sub> and dθ<sub>2</sub>/dt against time.");
}

function generateR1TPlot(solution) {
    var {t, vars} = solution;
    var r1 = vars[0];
    gen2DPlot(t, r1, "plotR1T", "Time plot of r<sub>1</sub>.", "t", "r<sub>1</sub>");
}

function generateDr1TPlot(solution) {
    var {t, vars} = solution;
    var dr1 = vars[1];
    gen2DPlot(t, dr1, "plotDr1T", "Time plot of dr<sub>1</sub>/dt.", "t", "dr<sub>1</sub>/dt");
}

function generateDr1R1Plot(solution) {
    var {vars} = solution;
    var r1 = vars[0];
    var dr1 = vars[1];
    gen2DPlot(r1, dr1, "plotDr1R1", "Phase plot of dr<sub>1</sub>/dt vs r<sub>1</sub>.", "r<sub>1</sub>", "dr<sub>1</sub>/dt");
}

function generateR2TPlot(solution) {
    var {t, vars} = solution;
    var r2 = vars[2];
    gen2DPlot(t, r2, "plotR2T", "Time plot of r<sub>2</sub>.", "t", "r<sub>2</sub>");
}

function generateDr2TPlot(solution) {
    var {t, vars} = solution;
    var dr2 = vars[3];
    gen2DPlot(t, dr2, "plotDr2T", "Time plot of dr<sub>2</sub>/dt.", "t", "dr<sub>2</sub>/dt");
}

function generateDr2R2Plot(solution) {
    var {vars} = solution;
    var r2 = vars[2];
    var dr2 = vars[3];
    gen2DPlot(r2, dr2, "plotDr2R2", "Phase plot of dr<sub>2</sub>/dt vs r<sub>2</sub>.", "r<sub>2</sub>", "dr<sub>2</sub>/dt");
}

function generateTheta1TPlot(solution) {
    var {t, vars} = solution;
    var theta1 = vars[4];
    gen2DPlot(t, theta1, "plotTheta1T", "Time plot of θ<sub>1</sub>.", "t", "θ<sub>1</sub>");
}

function generateDtheta1TPlot(solution) {
    var {t, vars} = solution;
    var dtheta1 = vars[5];
    gen2DPlot(t, dtheta1, "plotDtheta1T", "Time plot of dθ<sub>1</sub>/dt.", "t", "dθ<sub>1</sub>/dt");
}

function generateDtheta1Theta1Plot(solution) {
    var {vars} = solution;
    var theta1 = vars[4];
    var dtheta1 = vars[5];
    gen2DPlot(theta1, dtheta1, "plotDtheta1Theta1", "Phase plot of dθ<sub>1</sub>/dt vs θ<sub>1</sub>.", "θ<sub>1</sub>", "dθ<sub>1</sub>/dt");
}

function generateTheta2TPlot(solution) {
    var {t, vars} = solution;
    var theta2 = vars[6];
    gen2DPlot(t, theta2, "plotTheta2T", "Time plot of θ<sub>2</sub>.", "t", "θ<sub>2</sub>");
}

function generateDtheta2TPlot(solution) {
    var {t, vars} = solution;
    var dtheta2 = vars[7];
    gen2DPlot(t, dtheta2, "plotDtheta2T", "Time plot of dθ<sub>2</sub>/dt.", "t", "dθ<sub>2</sub>/dt");
}

function generateDtheta2Theta2Plot(solution) {
    var {vars} = solution;
    var theta2 = vars[6];
    var dtheta2 = vars[7];
    gen2DPlot(theta2, dtheta2, "plotDtheta2Theta2", "Phase plot of dθ<sub>2</sub>/dt vs θ<sub>2</sub>.", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt");
}

function generateR1R2PhasePlot(solution) {
    var {vars} = solution;
    var r1 = vars[0];
    var r2 = vars[2];
    gen2DPlot(r1, r2, "plotR2R1", "Phase plot of r<sub>2</sub> vs r<sub>1</sub>.", "r<sub>1</sub>", "r<sub>2</sub>");
}

/**
 * Generate all plots
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs, solution=undefined) {
    // Solve problem
    if (solution == undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    // Generate plots
    generatePendulumPlots(objectOfInputs, solution);
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

/**
 * Generate animation
 * @return nothing
 */
function generateAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animatePendulum(objectOfInputs, solution, "Double elastic pendulum");
}

function generateR1PhaseAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["r<sub>1</sub>", "dr<sub>1</sub>/dt"],IdSuffix: "R1Phase", title: "Phase plot of dr<sub>1</sub>/dt vs r<sub>1</sub>."});
}

function generateR2PhaseAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {timer: [1.0, 1.0], varnames: ["r<sub>2</sub>", "dr<sub>2</sub>/dt"], IdSuffix: "R2Phase", nos: [2, 3], title: "Phase plot of dr<sub>2</sub>/dt against r<sub>2</sub>."});
}

function generateTheta1PhaseAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt"],IdSuffix: "Theta1Phase", nos: [4, 5], title: "Phase plot of dθ<sub>1</sub>/dt against θ<sub>1</sub> ."});
}

function generateTheta2PhaseAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], timer: [0.0, 0.0], IdSuffix: "Theta2Phase", nos: [6, 7], title: "Phase plot of dθ<sub>2</sub>/dt against θ<sub>2</sub>."});
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    fillTable(objectOfInputs, ['r<sub>1</sub>', 'dr<sub>1</sub>/dt', 'r<sub>2</sub>', 'dr<sub>2</sub>/dt', '&theta;<sub>1</sub>', 'd&theta;<sub>1</sub>/dt', '&theta;<sub>2</sub>', 'd&theta;<sub>2</sub>/dt'], solution)
}

function generateAnimations(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    generateAnimation(objectOfInputs, solution)
    generateR1PhaseAnimation(objectOfInputs, solution);
    generateR2PhaseAnimation(objectOfInputs, solution);
    generateTheta1PhaseAnimation(objectOfInputs, solution);
    generateTheta2PhaseAnimation(objectOfInputs, solution);
}