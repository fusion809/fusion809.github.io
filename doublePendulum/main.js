/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [theta1, p1, theta2, p2]
 * @return               [dtheta1/dt, dp1/dt, dtheta2/dt, dp2/dt]
 */
function f(objectOfInputs, t, vars) {
    var {g, l1, l2, m1, m2} = objectOfInputs;
    var [theta1, ptheta1, theta2, ptheta2] = vars;

    // Simplifying terms
    var C1 = ptheta1*ptheta2*Math.sin(theta1-theta2)/(l1*l2*(m1+m2*(Math.sin(theta1-theta2)**2)));
    var C2 = ((l2**2)*m2*ptheta1**2 + (l1**2)*(m1+m2)*ptheta2**2 - l1*l2*m2*ptheta1*ptheta2*Math.cos(theta1-theta2))/(2*(l1**2)*(l2**2)*(m1+m2*Math.sin(theta1-theta2)**2))*Math.sin(2*(theta1-theta2));
    
    // Derivatives
    var thetaDot1 = (l2 * ptheta1 - l1*ptheta2 * Math.cos(theta1))/((l1**2)*l2*(m1+m2*((Math.sin(theta1-theta2))**2)));
    var pthetaDot1 = -(m1+m2)*g*l1*Math.sin(theta1) - C1 + C2;
    var thetaDot2 = (l1*(m1+m2)*ptheta2-l2*m2*ptheta1*Math.cos(theta1-theta2))/(l1*(l2**2)*m2*(m1+m2*(Math.sin(theta1-theta2)**2)));
    var pthetaDot2 = -m2*g*l2*Math.sin(theta2)+C1-C2;

    // Return statement
    return [thetaDot1, pthetaDot1, thetaDot2, pthetaDot2];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {theta10, p10, theta20, p20} = objectOfInputs;
    var vars0 = [[theta10, p10, theta20, p20]];
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
    gen2DPlot(theta1, theta2, "phasePlotTheta1Theta2", "Phase plot of theta2 against theta1");
}

/**
 * Generates a ptheta2 against ptheta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateP1P2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var ptheta1 = vars[1];
    var ptheta2 = vars[3];

    // Generate 2D plot
    gen2DPlot(ptheta1, ptheta2, "phasePlotP1P2", "Phase plot of ptheta2 against ptheta1");
}

/**
 * Generates a ptheta1 against theta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta1P1PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var ptheta1 = vars[1];
    
    // Generate 2D plot
    gen2DPlot(theta1, ptheta1, "phasePlotTheta1P1", "Phase plot of ptheta1 against theta1");
}

/**
 * Generates a ptheta2 against theta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta1P2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var ptheta2 = vars[2];
    
    // Generate 2D plot
    gen2DPlot(theta1, ptheta2, "phasePlotTheta1P2", "Phase plot of ptheta2 against theta1");
}

/**
 * Generates a ptheta1 against theta2 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta2P1PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta2 = vars[2];
    var ptheta1 = vars[1];
    
    // Generate 2D plot
    gen2DPlot(theta2, ptheta1, "phasePlotTheta2P1", "Phase plot of ptheta1 against theta2");
}

/**
 * Generates a ptheta2 against theta2 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta2P2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta2 = vars[2];
    var ptheta2 = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta2, ptheta2, "phasePlotTheta2P2", "Phase plot of ptheta2 against theta2");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["theta1", "ptheta1", "theta2", "ptheta2"], "timePlot", "Plot of theta1, ptheta1, theta2 and ptheta2 against time");
}

/**
 * Generate cartesian coordinates 
 * @param func           Function being used to integrate problem
 * @param objectOfInputs Problem parameters.
 */
function generatePendulumCoords(solution) {
    // Extract solution values and pendulum lengths
    var {t, vars} = solution;
    var [theta1, ptheta1, theta2, ptheta2] = vars;
    var {l1, l2} = objectOfInputs;
    var N = theta1.length;

    // Initialize arrays that will store x and y coords
    var x1 = new Array(N);
    var x2 = new Array(N);
    var y1 = new Array(N);
    var y2 = new Array(N);
    for (let i = 0; i < N; i++) {
        x1[i] = l1*Math.sin(theta1[i]);
        y1[i] = -l1*Math.cos(theta1[i]);
        x2[i] = x1[i] + l2*Math.sin(theta2[i]);
        y2[i] = y1[i] - l2*Math.cos(theta2[i]);
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
        name: "Pendulum 1 plot"
    }
    var plotPen2 = {
        x: x2,
        y: y2,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: "Pendulum 2 plot"
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
        name: "Pendulum 1 plot against time"
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
        name: 'Pendulum 2 plot against time'
    };
    var dataPen = [plotPen1Time, plotPen2Time];
    var layoutPenTime = {
        title: "Pendulum position against time plot"
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
    generateTheta1Theta2PhasePlot(solution);
    generateTheta1P1PhasePlot(solution);
    generateTheta1P2PhasePlot(solution);
    generateTheta2P1PhasePlot(solution);
    generateTheta2P2PhasePlot(solution);
    generateP1P2PhasePlot(solution);
    generatePendulumPlots(solution)
    generateTimePlot(solution);
}