function f(objectOfInputs, t, vars, dt) {
    var {g, l, mr, mb, br, cr, bb, cb} = objectOfInputs;
    var [theta, dtheta] = vars;
    // Calculate derivatives
    var mT = mb+mr/3;
    var mV = mb+mr/2;
    var d2theta = -g*mV/(mT*l)*Math.cos(theta) - (br+cr*l*Math.abs(dtheta)/2)/mT*dtheta/4 - (bb + cb*l*Math.abs(dtheta))/mT * dtheta; 

    // Put into return value
    return [dt*dtheta, dt*d2theta];
}

function prntObj(obj) {
    var keys = Object.keys(obj);
    var N = keys.length;
    var str = " ";
    for (let i =0; i<N; i++) {
        str += keys[i] + " = " + obj[keys[i]] + ", "
    }
    return str;
}

/**
 * Generate phase plot of theta dot against theta
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(solution) {
    // Extract solution data
    var {vars} = solution;
    var [theta, dtheta] = vars;

    // Generate 2D plot
    gen2DPlot(theta, dtheta, "phasePlot", "Phase plot of dθ/dt and θ", "θ", "dθ/dt");
}

/**
 * Generate plot of theta and theta dot against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(solution) {
    // Extract solution values
    var {t, vars} = solution;
    var theta = vars[0];
    var dtheta = vars[1];
    var sol = {
        t: t,
        vars: [theta, dtheta]
    }

    // Generate time plot
    genMultPlot(sol, ["θ", "dθ/dt"], "timePlot", "Plot of θ and dθ/dt against time");
}

/**
 * Generate two plots:
 * - one of dtheta and theta against t;
 * - a phase plot of dtheta against theta.
 * 
 * @param objectOfInputs An object that contains all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate solution plots
    generateTimePlot(solution);
    generatePhasePlot(solution);
};

/**
 * Generate x and y coordinates of pendulum.
 * @param objectOfInputs An object that contains all the problem parameters.
 * @param solution       solution object; solution.t is a vector of time
 *                       values, solution.vars contains theta and theta dot. 
 * @returns              [x, y] coordinates.
 */
function generatePendulumCoords(objectOfInputs, solution) {
    var theta = solution.vars[0];
    var l = objectOfInputs.l;
    var N = theta.length;
    var x = new Array(N);
    var y = new Array(N);
    for (let i = 0; i < N; i++) {
        x[i] = l*Math.cos(theta[i]);
        y[i] = l*Math.sin(theta[i]);
    }
    return [x, y];
}

/**
 * Generates animation.
 * @return nothing. 
 */
function generateAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animatePendulum(objectOfInputs, solution, "Single pendulum");
}

function generatePhaseAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    animate2D(solution, {varnames: ["θ", "dθ/dt"], IdSuffix: "Phase", title: "Single pendulum: phase plot of dθ/dt against θ."});
}

function generateAnimations(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    generateAnimation(objectOfInputs, solution);
    generatePhaseAnimation(solution);
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    fillTable(objectOfInputs, ['θ', 'dθ/dt'], solution)
}