function f(objectOfInputs, t, vars, dt) {
    var [theta1, dtheta1, theta2, dtheta2, theta3, dtheta3] = vars;
    var {l1, l2, l3, g,
    m1b, m1r, m2b, m2r, m3b, m3r,
    b1b, b1r, b2b, b2r, b3b, b3r,
    c1b, c1r, c2b, c2r, c3b, c3r} = objectOfInputs;
    var Delta21 = theta2-theta1;
    var Delta32 = theta3-theta2;
    var Delta31 = theta3 - theta1;
    var v2b = Math.sqrt(l1**2*dtheta1**2 + 2*l1*l2*dtheta1*dtheta2*Math.cos(Delta21)+l2**2*dtheta2**2);
    var v3b = Math.sqrt(l1**2*dtheta1**2 + l2**2*dtheta2**2 + l3**2*dtheta3**2 + 2*l1*l2*dtheta1*dtheta2*Math.cos(Delta21) + 2*l1*l3*dtheta1*dtheta3*Math.cos(Delta31)+2*l2*l3*dtheta2*dtheta3*Math.cos(Delta32))
    var v2r = Math.sqrt(l1**2*dtheta1**2 + l1*l2*dtheta1*dtheta2*Math.cos(Delta21)+l2**2*dtheta2**2/4);
    var v3r = Math.sqrt(l1**2*dtheta1**2 + l2**2*dtheta2**2 + l3**2*dtheta3**2/4 + 2*l1*l2*dtheta1*dtheta2*Math.cos(Delta21) + l1*l3*dtheta1*dtheta3*Math.cos(Delta31)+l2*l3*dtheta2*dtheta3*Math.cos(Delta32))
    var Qtheta1 = -(b1b + c1b*l1*Math.abs(dtheta1))*(l1**2*dtheta1) - (b1r+c1r*l1*Math.abs(dtheta1)/2)*(l1**2*dtheta1/4) - (b2b+c2b*v2b)*(l1**2*dtheta1 + l1*l2*dtheta2*Math.cos(Delta21)) - (b2r + c2r*v2r)*(l1**2*dtheta1+l1*l2*dtheta2*Math.cos(Delta21)/2) - (b3b+c3b*v3b) * (l1**2*dtheta1 + l1*l2*dtheta2*Math.cos(Delta21) + l1*l3*dtheta3*Math.cos(Delta31)) - (b3r+c3r*v3r)*(l1**2*dtheta1 + l1*l2*dtheta2*Math.cos(Delta21) + l1*l3*dtheta3*Math.cos(Delta31)/2)
    var Qtheta2 = -(b2b + c2b*v2b)*(l2**2*dtheta2 + l1*l2*dtheta1*Math.cos(Delta21)) - (b2r+c2r*v2r)*(l2**2*dtheta2/4 + l1*l2*dtheta1*Math.cos(Delta21)/2) - (b3b + c3b*v3b)*(l1*l2*dtheta1*Math.cos(Delta21)+l2**2*dtheta2+l2*l3*dtheta3*Math.cos(Delta32)) - (b3r+c3r*v3r)*(l1*l2*dtheta1*Math.cos(Delta21) + l2**2*dtheta2 + l2*l3*dtheta3*Math.cos(Delta32)/2);
    var Qtheta3 = -(b3b+c3b*v3b)*(l1*l3*dtheta1*Math.cos(Delta31) + l2*l3*dtheta2*Math.cos(Delta32) + l3**2*dtheta3) - (b3r+c3r*v3r)*(l1*l3*dtheta1*Math.cos(Delta31)/2 + l2*l3*dtheta2*Math.cos(Delta32)/2+l3**2*dtheta3/4)
    // Massive rods, self calc
    var M1 = m1b + m1r/3 + m2b + m2r + m3b + m3r;
    var M2 = m2b + m2r/3 + m3b + m3r;
    var M3 = m3b + m3r/3;
    var mu1 = m1b + m1r/2 + m2b + m2r + m3b + m3r;
    var mu2 = m2b + m2r/2 + m3b + m3r;
    var mu3 = m3b + m3r/2;
    var A = [
      [M1*l1**2, mu2*l1*l2*Math.cos(Delta21), mu3*l1*l3*Math.cos(Delta31)],
      [mu2*l2*l1*Math.cos(Delta21), M2*l2**2, mu3*l2*l3*Math.cos(Delta32)],
      [mu3*l3*l1*Math.cos(Delta31), mu3*l3*l2*Math.cos(Delta32), M3*l3**2]
    ];
    var b = [
      Qtheta1 - mu1*g*l1*Math.cos(theta1) + mu2*l1*l2*dtheta2**2*Math.sin(Delta21) + mu3*l1*l3*dtheta3**2*Math.sin(Delta31),
      Qtheta2 - mu2*l2*(l1*dtheta1**2*Math.sin(Delta21)+g*Math.cos(theta2)) + mu3*l2*l3*dtheta3**2*Math.sin(Delta32),
      Qtheta3 - mu3*l3*(l1*dtheta1**2*Math.sin(Delta31) + l2*dtheta2**2*Math.sin(Delta32)) - mu3*g*l3*Math.cos(theta3)
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
 * Generate a theta1 vs time plot
 * 
 */
function generateTheta1TPlot(solution) {
  var {t, vars} = solution;
  var theta1 = vars[0];
  gen2DPlot(t, theta1, "theta1TPlot", "θ<sub>1</sub> against time", "t", "θ<sub>1</sub>");
}

/**
 * Generate a dtheta1 vs time plot
 * 
 */
function generateDtheta1TPlot(solution) {
  var {t, vars} = solution;
  var dtheta1 = vars[1];
  gen2DPlot(t, dtheta1, "dtheta1TPlot", "dθ<sub>1</sub>/dt against time", "t", "dθ<sub>1</sub>/dt");
}

/**
 * Generate a theta2 vs time plot
 * 
 */
function generateTheta2TPlot(solution) {
  var {t, vars} = solution;
  var theta2 = vars[2];
  gen2DPlot(t, theta2, "theta2TPlot", "θ<sub>2</sub> against time", "t", "θ<sub>2</sub>");
}

/**
 * Generate a dtheta2 vs time plot
 * 
 */
function generateDtheta2TPlot(solution) {
  var {t, vars} = solution;
  var dtheta2 = vars[3];
  gen2DPlot(t, dtheta2, "dtheta2TPlot", "dθ<sub>2</sub>/dt against time", "t", "dθ<sub>2</sub>/dt");
}

/**
 * Generate a dtheta3 vs time plot
 * 
 */
function generateDtheta3TPlot(solution) {
  var {t, vars} = solution;
  var dtheta3 = vars[5];
  gen2DPlot(t, dtheta3, "dtheta3TPlot", "dθ<sub>3</sub>/dt against time", "t", "dθ<sub>3</sub>/dt");
}

/**
 * Generate a theta3 vs time plot
 * 
 */
function generateTheta3TPlot(solution) {
  var {t, vars} = solution;
  var theta3 = vars[4];
  gen2DPlot(t, theta3, "theta3TPlot", "θ<sub>3</sub> against time", "t", "θ<sub>3</sub>");
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
    generateTheta1TPlot(solution);
    generateTheta2TPlot(solution);
    generateTheta3TPlot(solution);
    generateDtheta1TPlot(solution);
    generateDtheta2TPlot(solution);
    generateDtheta3TPlot(solution);
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

function generateThetasPhaseAnimation(objectOfInputs=undefined, solution=undefined) {
  if (objectOfInputs == undefined) {
    var objectOfInputs = readInputs();
  } 
  if (solution == undefined) {
    var solution = solveProblem(RKF45, objectOfInputs);
  }
  animate3D(solution, {view: [2, 0, 0], IdSuffix: "ThetasPhase", varnames: ['θ<sub>1</sub>', 'θ<sub>2</sub>', 'θ<sub>3</sub>'], nos: [0, 2, 4], title: "Triple pendulum: 3D phase animation."});
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
  generateThetasPhaseAnimation(objectOfInputs, solution);
}