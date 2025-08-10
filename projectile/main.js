/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param t          Time (seconds).
 */
function f(objectOfInputs, t, vars, dt) {
    var {g, b, c} = objectOfInputs;
    var [x, dx, y, dy] = vars;
    var v = Math.sqrt(dx**2+dy**2);
    if ((y <= 0) && (dy <= 0)) {
        dy = 0
        var d2y = 0
    } else {
        var d2y = -(g+dy*(b+c*v));
    }
    
    return [dt*dx, -dt*dx*(b+c*v), dt*dy, dt*d2y];
}

function trimAndSol(t, x, dx, y, dy) {
    var j;
    for (let i = 0; i < y.length; i++) {
        if (y[i] < 0 && dy[i] < 0) {
            if (j == undefined) {
                j = i;
                break;
            }
        }
    }
    x = x.slice(0, j);
    dx = dx.slice(0, j);
    y = y.slice(0, j);
    dy = dy.slice(0, j);
    t = t.slice(0, j);
    var vars = [x, dx, y, dy]
    var solution = new SolClass(t, vars);
    return solution;
}

function solve(objectOfInputs=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    var N = objectOfInputs.N;
    if (objectOfInputs.b == 0 && objectOfInputs.c == 0) {
        var t = linspace(0, objectOfInputs.tf, N);
        var x = t.map(val => objectOfInputs.dx0*val + objectOfInputs.x0);
        var dx = Array(N).fill(objectOfInputs.dx0);
        var y = t.map(val => -objectOfInputs.g/2 * val**2 + objectOfInputs.y0 + objectOfInputs.dy0 * val)
        var dy = t.map(val => -objectOfInputs.g*val + objectOfInputs.dy0)
        var solution = trimAndSol(t, x, dx, y, dy);
    } else if (objectOfInputs.c == 0) {
        var t = linspace(0, objectOfInputs.tf, N);
        var x = t.map(val => objectOfInputs.dx0/objectOfInputs.b * (1-Math.exp(-objectOfInputs.b*val)) + objectOfInputs.x0);
        var dx = t.map(val => Math.exp(-objectOfInputs.b*val) + objectOfInputs.dx0);
        var dy = t.map(val => -objectOfInputs.g*val + objectOfInputs.dy0*Math.exp(-objectOfInputs.b*val))
        var y = t.map(val => -objectOfInputs.g/2*val**2 + objectOfInputs.y0 + (objectOfInputs.dy0/objectOfInputs.b)*(1-Math.exp(-objectOfInputs.b*val)))
        var solution = trimAndSol(t, x, dx, y, dy);
    } else {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    return solution;
}

/**
 * Generate phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the phase plot.
 */
function generateXPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, dx] = vars;

    // Generate 2D phase plot
    gen2DPlot(x, dx, "xPhasePlot", "Phase plot of dx/dt vs x.", "x", "dx/dt")
}

/**
 * Generate phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the phase plot.
 */
function generateYPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, dx, y, dy] = vars;

    // Generate 2D phase plot
    gen2DPlot(y, dy, "yPhasePlot", "Phase plot of dy/dt vs y.", "y", "dy/dt")
}

/**
 * Generate phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the phase plot.
 */
function generateXYPlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, dx, y, dy] = vars;

    // Generate 2D phase plot
    gen2DPlot(x, y, "xYPlot", "Plot of projectile motion.", "x", "y")
}

/**
 * Generate a plot of x and y against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates a plot of x and y against time.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "y"], "timePlot", "Plot of prey and predator animal numbers against time");
}

/**
 * Generate two plots:
 * - one of y and x against t; and
 * - a phase plot of y against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solve(objectOfInputs);
    }

    // Generate plots
    generateXPhasePlot(solution);
    generateYPhasePlot(solution);
    generateXYPlot(solution);
    generateTimePlot(solution);
}

/**
 * Generate animation
 * @return nothing.
 */
function generateXAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solve(objectOfInputs);
    }
    animate2D(solution, {timer: [0.98, 0.98], varnames: ["x", "dx/dt"], title: "Phase animation of dx/dt vs x.", IdSuffix: "X"});
}

/**
 * Generate animation
 * @return nothing.
 */
function generateYAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solve(objectOfInputs);
    }
    animate2D(solution, {timer: [0.98, 0.98], nos: [2, 3], varnames: ["y", "dy/dt"], title: "Phase animation of dy/dt vs y.", IdSuffix: "Y"});
}

/**
 * Generate animation
 * @return nothing.
 */
function generateXYAnimation(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solve(objectOfInputs);
    }
    animate2D(solution, {timer: [0.98, 0.98], nos: [0, 2], varnames: ["x", "y"], title: "Animation of projectile motion.", IdSuffix: "XY"});
}

function generateAnimations(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solve(objectOfInputs);
    }
    generateXAnimation(objectOfInputs, solution);
    generateYAnimation(objectOfInputs, solution);
    generateXYAnimation(objectOfInputs, solution);
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solve(objectOfInputs);
    }
    fillTable(objectOfInputs, ['x', 'dx/dt', 'y', 'dy/dt'], solution)
}