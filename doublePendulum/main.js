/**
 * Right-hand side of our second-order ODE written as a simple of first-order ODEs.
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
    var thetaDot2 = (l1*(m1+m2)*ptheta2-l2*m2*ptheta1*Math.cos(theta1-theta2))/(l1*(l2**2)*m2*(m1+m2*Math.sin(theta1-theta2)**2));
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
    var [t, vars] = RKF45Body(objectOfInputs, vars0); 
    return [t, vars];
}

/**
 * Generates a 2D phase plot of theta2 against theta1
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateTheta1Theta2PhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);

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
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateP1P2PhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);

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
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateTheta1P1PhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);
    
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var ptheta1 = vars[1];
    
    // Generate 2D plot
    gen2DPlot(theta1, ptheta1, "phasePlotTheta1P1", "Phase plot of ptheta1 against theta1");
}

/**
 * Generates a ptheta2 against theta2 phase plot
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateTheta2P2PhasePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);
    
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
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing.
 */
function generateTimePlot(objectOfInputs) {
    // Run solveProblem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate time plot
    genMultPlot(solution, ["theta1", "ptheta1", "theta2", "ptheta2"], "timePlot", "Plot of theta1, ptheta1, theta2 and ptheta2 against time");
}

/**
 * Generate all plots
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generateTheta1Theta2PhasePlot(objectOfInputs);
    generateTheta1P1PhasePlot(objectOfInputs);
    generateP1P2PhasePlot(objectOfInputs);
    generateTheta2P2PhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
}

/**
 * Remove theta1 theta2 plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeTheta1Theta2PhasePlot() {
    rmPlot("phasePlotTheta1Theta2");
}

/**
 * Remove theta1 ptheta1 plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeTheta1P1PhasePlot() {
    rmPlot("phasePlotTheta1P1");
}

/**
 * Remove ptheta1 ptheta2 plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeP1P2PhasePlot() {
    rmPlot("phasePlotP1P2");
}

/**
 * Remove theta1 ptheta1 plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeTheta2P2PhasePlot() {
    rmPlot("phasePlotTheta2P2");
}