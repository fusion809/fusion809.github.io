function f(objectOfInputs, t, vars, dt) {
    var [theta1, dtheta1, theta2, dtheta2, theta3, dtheta3] = vars;
    var {b1, c1, b2, c2, b3, c3, l1, l2, l3, m1, m2, m3, g} = objectOfInputs;
    var Delta21 = theta2-theta1;
    var Delta32 = theta3-theta2;
    var Delta31 = theta3 - theta1;
    var v1 = l1*Math.abs(dtheta1);
    var v2 = Math.sqrt(l1**2*dtheta1**2 + 2*l1*l2*dtheta1*dtheta2*Math.cos(Delta21)+l2**2*dtheta2**2);
    var v3 = Math.sqrt(l1**2*dtheta1**2 + l2**2*dtheta2**2 + l3**2*dtheta3**2 + 2*l1*l2*dtheta1*dtheta2*Math.cos(Delta21) + 2*l1*l3*dtheta1*dtheta3*Math.cos(Delta31)+2*l2*l3*dtheta2*dtheta3*Math.cos(Delta32))
    var m123 = m1 + m2 + m3;
    var m23 = m2 + m3;
    var Qtheta1 = -(b1+c1*v1)*(l1**2*dtheta1) - (b2+c2*v2)*(l1**2*dtheta1 + l1*l2*dtheta2*Math.cos(Delta21))-(b3+c3*v3)*(l1**2*dtheta1 + l1*l2*dtheta2*Math.cos(Delta21) + l1*l3*dtheta3*Math.cos(Delta31));
    var Qtheta2 = -(b2+c2*v2)*(l2**2*dtheta2 + l1*l2*dtheta1*Math.cos(Delta21)) - (b3+c3*v3)*(l2**2*dtheta2 + l1*l2*dtheta1*Math.cos(Delta21) + l2*l3*dtheta3*Math.cos(Delta32))
    var Qtheta3 = -(b3+c3*v3)*(l3**2*dtheta3 + l1*l3*dtheta1*Math.cos(Delta31) + l2*l3*dtheta2*Math.cos(Delta32));
    // var A = [
    //   [m123*l1**2, m23*l1*l2*Math.cos(Delta21), m3*l1*l3*Math.cos(Delta31)],
    //   [m23*l1*l2*Math.cos(Delta21), m23*l2**2, m3*l2*l3*Math.cos(Delta32)],
    //   [m3*l1*l3*Math.cos(Delta31), m3*l2*l3*Math.cos(Delta32), m3*l3**2]
    // ]
    // var b = [
    //   Qtheta1 - m123*g*l1*Math.cos(theta1) + m23*l1*l2*dtheta2**2*Math.sin(Delta21)+m3*l1*l3*dtheta3**2*Math.sin(Delta31),
    //   Qtheta2 - m23*l1*l2*dtheta1**2*Math.sin(Delta21) - m23*l2*g*Math.cos(theta2) + m3*l2*l3*dtheta3**2*Math.sin(Delta32),
    //   Qtheta3 - m3*l1*3*dtheta1**2*Math.sin(Delta31) - m3*l2*l3*dtheta2**2*Math.sin(Delta32) - m3*l3*g*Math.cos(theta3)
    // ]
    var A = [
      [m1*(2*l1**2*Math.sin(theta1)**2 + 2*l1**2*Math.cos(theta1)**2)/2 + m2*(2*l1**2*Math.sin(theta1)**2 + 2*l1**2*Math.cos(theta1)**2)/2 + m3*(2*l1**2*Math.sin(theta1)**2 + 2*l1**2*Math.cos(theta1)**2)/2, m2*(2*l1*l2*Math.sin(theta1)*Math.sin(theta2) + 2*l1*l2*Math.cos(theta1)*Math.cos(theta2))/2 + m3*(2*l1*l2*Math.sin(theta1)*Math.sin(theta2) + 2*l1*l2*Math.cos(theta1)*Math.cos(theta2))/2, m3*(2*l1*l3*Math.sin(theta1)*Math.sin(theta3) + 2*l1*l3*Math.cos(theta1)*Math.cos(theta3))/2],
      [m2*(2*l1*l2*Math.sin(theta1)*Math.sin(theta2) + 2*l1*l2*Math.cos(theta1)*Math.cos(theta2))/2 + m3*(2*l1*l2*Math.sin(theta1)*Math.sin(theta2) + 2*l1*l2*Math.cos(theta1)*Math.cos(theta2))/2,                                                 m2*(2*l2**2*Math.sin(theta2)**2 + 2*l2**2*Math.cos(theta2)**2)/2 + m3*(2*l2**2*Math.sin(theta2)**2 + 2*l2**2*Math.cos(theta2)**2)/2, m3*(2*l2*l3*Math.sin(theta2)*Math.sin(theta3) + 2*l2*l3*Math.cos(theta2)*Math.cos(theta3))/2],
      [m3*(2*l1*l3*Math.sin(theta1)*Math.sin(theta3) + 2*l1*l3*Math.cos(theta1)*Math.cos(theta3))/2, m3*(2*l2*l3*Math.sin(theta2)*Math.sin(theta3) + 2*l2*l3*Math.cos(theta2)*Math.cos(theta3))/2, m3*(2*l3**2*Math.sin(theta3)**2 + 2*l3**2*Math.cos(theta3)**2)/2]
    ]
    var b = [
      -l1*g*m1*Math.cos(theta1) - l1*g*m2*Math.cos(theta1) - l1*g*m3*Math.cos(theta1) + m2*(-2*l1*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2)*Math.cos(theta1)*dtheta1 - 2*l1*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2)*Math.sin(theta1)*dtheta1)/2 - m2*(-2*l1*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2)*Math.cos(theta1)*dtheta1 + 2*l1*(-l1*Math.sin(theta1)*dtheta1**2 - l2*Math.sin(theta2)*dtheta2**2)*Math.cos(theta1) - 2*l1*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2)*Math.sin(theta1)*dtheta1 - 2*l1*(-l1*Math.cos(theta1)*dtheta1**2 - l2*Math.cos(theta2)*dtheta2**2)*Math.sin(theta1))/2 + m3*(-2*l1*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2 - l3*Math.sin(theta3)*dtheta3)*Math.cos(theta1)*dtheta1 - 2*l1*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2 + l3*Math.cos(theta3)*dtheta3)*Math.sin(theta1)*dtheta1)/2 - m3*(-2*l1*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2 - l3*Math.sin(theta3)*dtheta3)*Math.cos(theta1)*dtheta1 + 2*l1*(-l1*Math.sin(theta1)*dtheta1**2 - l2*Math.sin(theta2)*dtheta2**2 - l3*Math.sin(theta3)*dtheta3**2)*Math.cos(theta1) - 2*l1*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2 + l3*Math.cos(theta3)*dtheta3)*Math.sin(theta1)*dtheta1 - 2*l1*(-l1*Math.cos(theta1)*dtheta1**2 - l2*Math.cos(theta2)*dtheta2**2 - l3*Math.cos(theta3)*dtheta3**2)*Math.sin(theta1))/2 + Qtheta1, m2*(-2*l2*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2)*Math.cos(theta2)*dtheta2 - 2*l2*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2)*Math.sin(theta2)*dtheta2)/2 - m2*(-2*l2*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2)*Math.cos(theta2)*dtheta2 + 2*l2*(-l1*Math.sin(theta1)*dtheta1**2 - l2*Math.sin(theta2)*dtheta2**2)*Math.cos(theta2) - 2*l2*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2)*Math.sin(theta2)*dtheta2 - 2*l2*(-l1*Math.cos(theta1)*dtheta1**2 - l2*Math.cos(theta2)*dtheta2**2)*Math.sin(theta2))/2 + m3*(-2*l2*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2 - l3*Math.sin(theta3)*dtheta3)*Math.cos(theta2)*dtheta2 - 2*l2*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2 + l3*Math.cos(theta3)*dtheta3)*Math.sin(theta2)*dtheta2)/2 - m3*(-2*l2*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2 - l3*Math.sin(theta3)*dtheta3)*Math.cos(theta2)*dtheta2 + 2*l2*(-l1*Math.sin(theta1)*dtheta1**2 - l2*Math.sin(theta2)*dtheta2**2 - l3*Math.sin(theta3)*dtheta3**2)*Math.cos(theta2) - 2*l2*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2 + l3*Math.cos(theta3)*dtheta3)*Math.sin(theta2)*dtheta2 - 2*l2*(-l1*Math.cos(theta1)*dtheta1**2 - l2*Math.cos(theta2)*dtheta2**2 - l3*Math.cos(theta3)*dtheta3**2)*Math.sin(theta2))/2 + Qtheta2, m3*(-2*l3*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2 - l3*Math.sin(theta3)*dtheta3)*Math.cos(theta3)*dtheta3 - 2*l3*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2 + l3*Math.cos(theta3)*dtheta3)*Math.sin(theta3)*dtheta3)/2 - m3*(-2*l3*(-l1*Math.sin(theta1)*dtheta1 - l2*Math.sin(theta2)*dtheta2 - l3*Math.sin(theta3)*dtheta3)*Math.cos(theta3)*dtheta3 + 2*l3*(-l1*Math.sin(theta1)*dtheta1**2 - l2*Math.sin(theta2)*dtheta2**2 - l3*Math.sin(theta3)*dtheta3**2)*Math.cos(theta3) - 2*l3*(l1*Math.cos(theta1)*dtheta1 + l2*Math.cos(theta2)*dtheta2 + l3*Math.cos(theta3)*dtheta3)*Math.sin(theta3)*dtheta3 - 2*l3*(-l1*Math.cos(theta1)*dtheta1**2 - l2*Math.cos(theta2)*dtheta2**2 - l3*Math.cos(theta3)*dtheta3**2)*Math.sin(theta3))/2 + Qtheta3
    ]
    var d2 = math.lusolve(A, b);
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
 * Generates a dtheta2 against theta2 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateDtheta3Theta3PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta3 = vars[4];
    var dtheta3 = vars[5];
    
    // Generate 2D plot
    gen2DPlot(theta3, dtheta3, "phasePlotDtheta3Theta3", "Phase plot of dθ<sub>3</sub>/dt against θ<sub>3</sub>", "θ<sub>3</sub>", "dθ<sub>3</sub>/dt");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt", "θ<sub>3</sub>", "dθ<sub>3</sub>/dt"], "timePlot", "Plot of θ<sub>1</sub>, dθ<sub>1</sub>/dt, θ<sub>2</sub>, dθ<sub>2</sub>/dt, θ<sub>3</sub> and dθ<sub>3</sub>/dt against time");
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
    generateDtheta3Theta3PhasePlot(solution);
    generatePendulumPlots(objectOfInputs, solution)
    generateTimePlot(solution);
}

/**
 * Generate animation
 * @return nothing
 */
function generateAnimation(objectOfInputs=undefined, solution=undefined) {
  if (objectOfInputs == undefined) {
    var objectOfInputs = readInputs();
  } 
  if (solution == undefined) {
    var solution = solveProblem(RKF45, objectOfInputs);
  }
  animatePendulum(objectOfInputs, solution, "Triple pendulum");
}

function generateTheta1PhaseAnimation(objectOfInputs=undefined, solution=undefined) {
  if (objectOfInputs == undefined) {
    var objectOfInputs = readInputs();
  } 
  if (solution == undefined) {
    var solution = solveProblem(RKF45, objectOfInputs);
  }
  animate2D(solution, {varnames: ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt"], IdSuffix: "Theta1Phase", title: "Phase plot of dθ<sub>1</sub>/dt against θ<sub>1</sub>."});
}

function generateTheta2PhaseAnimation(objectOfInputs=undefined, solution=undefined) {
  if (objectOfInputs == undefined) {
    var objectOfInputs = readInputs();
  } 
  if (solution == undefined) {
    var solution = solveProblem(RKF45, objectOfInputs);
  }
  animate2D(solution, {timer: [0, 1.0], varnames: ["θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], IdSuffix: "Theta2Phase", nos: [2, 3], title: "Phase plot of dθ<sub>2</sub>/dt against θ<sub>2</sub>."});
}

function generateTheta3PhaseAnimation(objectOfInputs=undefined, solution=undefined) {
  if (objectOfInputs == undefined) {
    var objectOfInputs = readInputs();
  } 
  if (solution == undefined) {
    var solution = solveProblem(RKF45, objectOfInputs);
  }
  animate2D(solution, {timer: [0, 1.0], varnames: ["θ<sub>3</sub>", "dθ<sub>3</sub>/dt"], IdSuffix: "Theta3Phase", nos: [4, 5], title: "Phase plot of dθ<sub>3</sub>/dt against θ<sub>3</sub>."});
}

function generateTable(objectOfInputs=undefined, solution=undefined) {
  if (objectOfInputs == undefined) {
    var objectOfInputs = readInputs();
  } 
  if (solution == undefined) {
    var solution = solveProblem(RKF45, objectOfInputs);
  }
  fillTable(objectOfInputs, ['&theta;<sub>1</sub>', 'd&theta;<sub>1</sub>/dt', '&theta;<sub>2</sub>', 'd&theta;<sub>2</sub>/dt', '&theta;<sub>3</sub>', 'd&theta;<sub>3</sub>/dt'], solution)
}

function generateAnimations(objectOfInputs=undefined, solution=undefined) {
  if (objectOfInputs == undefined) {
    var objectOfInputs = readInputs();
  } 
  if (solution == undefined) {
    var solution = solveProblem(RKF45, objectOfInputs);
  }
  generateAnimation(objectOfInputs, solution)
  generateTheta1PhaseAnimation(objectOfInputs, solution);
  generateTheta2PhaseAnimation(objectOfInputs, solution);
  generateTheta3PhaseAnimation(objectOfInputs, solution);
}