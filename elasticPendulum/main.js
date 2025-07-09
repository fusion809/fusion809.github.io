/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [x, xDot, theta, thetaDot]
 * @return               [dx/dt, d2x/dt2, dtheta/dt, d2theta/dt2]
 */
function f(objectOfInputs, t, vars, dt) {
    var {g, l0, k, m} = objectOfInputs;
    var [x, xDot, theta, thetaDot] = vars;
    var xDDot = (l0+x)*thetaDot**2 - k*x/m + g*Math.sin(theta);
    var thetaDDot = -g*Math.cos(theta)/(l0+x)-2*xDot*thetaDot/(l0+x);
    return [dt*xDot, dt*xDDot, dt*thetaDot, dt*thetaDDot];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {x0, xDot0, theta0, thetaDot0} = objectOfInputs;
    var vars0 = [[x0, xDot0, theta0, thetaDot0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0); 
    return [t, vars];
}

/**
 * Generates a 2D phase plot of x against theta
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXThetaPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var theta = vars[2];

    // Generate 2D plot
    gen2DPlot(x, theta, "phasePlotXTheta", "Phase plot of theta against x");
}

/**
 * Generates a xdot against x phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateXXDotPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var xDot = vars[1];

    // Generate 2D plot
    gen2DPlot(x, xDot, "phasePlotXXDot", "Phase plot of x dot against x");
}

/**
 * Generates a theta dot against theta phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateThetaThetaDotPhasePlot(solution) {    
    // Extract solution data from solution object
    var {vars} = solution;
    var theta = vars[2];
    var thetaDot = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta, thetaDot, "phasePlotThetaThetaDot", "Phase plot of theta dot against theta");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "x dot", "theta", "theta dot"], "timePlot", "Plot of x, x dot, theta and theta dot against time");
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
    generateXThetaPhasePlot(solution);
    generateXXDotPhasePlot(solution);
    generateThetaThetaDotPhasePlot(solution);
    generateTimePlot(solution);
}

function generatePendulumCoords(r, th) {
    var N = th.length;
    var x = new Array(N);
    var y = new Array(N);
    for (let i = 0; i < N; i++) {
        x[i] = r[i]*Math.cos(th[i]);
        y[i] = r[i]*Math.sin(th[i]);
    }
    return [x, y];
}

function animatePendulum(solution) {
  const t = solution.t;
  const vars = solution.vars;
  const th = vars[2];
  const z = vars[0];
  var [x, y] = generatePendulumCoords(z, th);
  var xmin = Math.min(...x);
  var ymin = Math.min(...y);
  var xmax = Math.max(...x);
  var ymax = Math.max(...y);
  
  const trace1 = {
    x: [], y: [],
    mode: "lines+markers",
    marker: { size: 8 },
    line: { color: "blue" },
    name: "Rod"
  };
  const data = [trace1];
  const padding = 1;
    if (isFinite(xmin)) {
    xmin = xmin - padding;
  }
  if (isFinite(xmax)) {
    xmax = xmax + padding;
  }

  if (isFinite(ymin)) {
    ymin = ymin - padding;
  }

  if (isFinite(ymax)) {
    ymax = ymax + padding;
  }

  const layout = {
    title: "Elastic pendulum",
    xaxis: { range: [xmin, xmax], title: "x" },
    yaxis: { range: [ymin, ymax], title: "y", scaleanchor: "x" },
    showlegend: false,
    annotations: [{
    x: 0,
    y: 1.1,  // slightly above plot
    xref: 'paper',
    yref: 'paper',
    text: 'Time: 0.00 s',
    showarrow: false,
    font: { size: 16 }
  }]
  };
  Plotly.newPlot("animation", data, layout).then(() => {
    let frame = 0;
    let startTime = null;
    let paused = false;
    let pauseStart = 0;
    let totalPausedTime = 0;

    const button = document.getElementById("toggleButton");
    button.addEventListener("click", () => {
      paused = !paused;
      button.textContent = paused ? "Play" : "Pause";
      if (paused) {
        pauseStart = Date.now();
      } else {
        totalPausedTime += Date.now() - pauseStart;
        if (!paused) requestAnimationFrame(animateFrame);
      }
    });

    let cycle = 0;
    function animateFrame(timestamp) {
      if (frame == t.length - 1 || !startTime) {
        frame = 0;
        startTime = timestamp; 
      }
      if (paused) return;

      const elapsedSec = (timestamp - startTime - totalPausedTime) / 1000;

      while (frame < t.length - 1 && t[frame] < elapsedSec) {
        frame++;
      }
      if (frame >= t.length) {
        frame = t.length - 1;
        cycle++;
      }

    Plotly.animate("animation", {
      data: [
        { x: [0, x[frame]], y: [0, y[frame]] }
      ]
    }, {
      transition: { duration: 0 },
      frame: { duration: 0, redraw: true }
    });
    layout.annotations[0].text = `Time: ${t[frame].toFixed(2)} s`;
    Plotly.relayout('animation', layout);

    requestAnimationFrame(animateFrame);
  }

  requestAnimationFrame(animateFrame);
});
}

function removeAnimation() {
    rmPlot("animation");
}

function animSim() {
    var solution = solveProblem(RKF45, readInputs());
    animatePendulum(solution);
}