/** 
* Solve the problem using RKF45
*
* @param objectOfInputs An object containing all the problem parameters.
* @return               [t, vars]
*/
RKF45 = function(objectOfInputs) {
    // Extract initial conditions from object and enter it into RKF45Body
    var {x0, y0, z0} = objectOfInputs;
    var vars0 = [[x0, y0, z0]];
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
    // Extract solution data from solution object
    var vars = solution.vars;
    var [x, y, z] = vars;

    gen3DPlot(x, y, z, phasePlotXYZ, Phase plot of the solution to the Lorenz system, {view: [2, 0, 0]})
}

/**
* Generates a XY phase plot
* 
* @param solution       An object containing solution data.
* @return               Nothing.
*/
function generateXYPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var y = vars[1];

    // Generate plot
    gen2DPlot(x, y, phasePlotXY, y against x phase plot, x, y);
}

/**
* Generates a XZ phase plot
* 
* @param solution       An object containing solution data.
* @return               Nothing.
*/
function generateXZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var z = vars[2];

    // Generate plot
    gen2DPlot(x, z, phasePlotXZ, z against x phase plot, x, z);
}

/**
* Generates a YZ phase plot
* 
* @param solution       An object containing solution data.
* @return               Nothing.
*/
function generateYZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var y = vars[1];
    var z = vars[2];

    // Generate plot
    gen2DPlot(y, z, phasePlotYZ, z against y phase plot, y, z);
}

/**
* Generates a time plot
* 
* @param solution       An object containing solution data.
* @return               Nothing.
*/
function generateTimePlot(solution) {
    // Extract solution data from solution object
    genMultPlot(solution, [x, y, z], timePlot, Time plots of the solution to the problem)
}

/**
* Remove animation
* @return nothing
*/
function removeAnimation() {
    rmPlot(animation);
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
    animate3D(solution, {view: [2, 0, 0], title: Lorenz system: x, y and z phase plot.});
}

/**
* Generate five plots:
* - The first is a 3D phase plot of x, y and z.
* - The second is a 2D phase plot of y against x.
* - The third is a 2D phase plot of z against x.
* - The fourth is a 2D phase plot of z against y.
* - The fifth is a plot of x, y and z against time.
* 
* @param objectOfInputs An object containing all the problem parameters.
* @return               Nothing. Just generates the plots.
*/
function generatePlots(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    // Generate plots
    generate3DPhasePlot(solution);
    generateXYPhasePlot(solution);
    generateXZPhasePlot(solution);
    generateYZPhasePlot(solution);
    generateTimePlot(solution);
}

function generateAllOutputs(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    fillTable(objectOfInputs, ['x', 'y', 'z'], solution)
    generatePlots(objectOfInputs, solution);
    generateAnimation(objectOfInputs, solution);
}

function removeAllOutputs() {
    removeTable();
    removePlots();
    removeAnimation();
}