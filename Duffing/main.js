/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param objectOfInputs Problem parameter.
 * @param t              Time (seconds).
 * @param x              x coordinate.
 * @param xDot           dx/dt
 * @return               [dx/dt, d2x/dt2]
 */
function f(objectOfInputs, t, vars, dt) {
    var {alpha, beta, gamma, delta, omega} = objectOfInputs;
    var [x, xDot] = vars;
    return [dt*xDot, dt*(- delta*xDot - alpha*x - beta*x**3 + gamma * Math.cos(omega*t))];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions and enter into 2d array
    var {x0, xDot0} = objectOfInputs;
    var vars0 = [[x0, xDot0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generate phase plot of x dot against x
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, xDot] = vars;

    // Generate plot
    gen2DPlot(x, xDot, "phasePlot", "Phase plot of x dot against x");
}

/**
 * Generate plot of x and x dot against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(solution) {
    // Plot
    genMultPlot(solution, ["x", "x dot"], "timePlot", "Plot of x dot and x against time");
}

/**
 * Generate two plots:
 * - one of xDot and x against t; and
 * - a phase plot of xDot against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    
    // Generate plots
    generateTimePlot(solution);
    generatePhasePlot(solution);
}

function animate(solution) {
  const x = solution.vars[0];
  const xdot = solution.vars[1];
  const t = solution.t;

  const tracePath = {
    x: x,
    y: xdot,
    mode: 'lines',
    type: 'scatter',
    line: { color: 'blue', width: 4 },
    name: 'Path',
  };

  const traceMarker = {
    x: [x[0]],
    y: [xdot[0]],
    mode: 'markers',
    type: 'scatter',
    marker: { color: 'red', size: 5 },
    name: 'Object'
  };

  const layout = {
    margin: { l: 0, r: 0, b: 0, t: 0 },
    scene: {
      xaxis: { title: 'x' },
      yaxis: { title: 'dx/dt' },
    },
    annotations: [{
      text: `Time: 0.00 s`,
      xref: 'paper',
      yref: 'paper',
      x: 0.05,
      y: 0.95,
      showarrow: false,
      font: { size: 16 }
    }],
    showlegend: true
  };

  Plotly.newPlot('animation', [tracePath, traceMarker], layout).then(() => {
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
      Plotly.animate('animation', {
        data: [
          tracePath,
          {
            x: [x[frame]],
            y: [xdot[frame]],
            mode: 'markers',
            type: 'scatter',
            marker: { color: 'red', size: 5 },
            name: 'Object'
          }
        ],
        layout: {
          annotations: [{
            text: `Time: ${t[frame].toFixed(2)} s`,
            xref: 'paper',
            yref: 'paper',
            x: 0.02,
            y: 0.98,
            showarrow: false,
            font: { size: 16 }
          }],
          scene: {
            xaxis: { title: 'x' },
            yaxis: { title: 'dx/dt' },
            camera: {
                eye: {
                x: 0.5,   // Set to 0 to align the camera with the YZ plane
                y: -2
                }
            }
        }
        }
      }, {
        transition: { duration: 0 },
        frame: { duration: 0, redraw: true }
      });

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
    animate(solution);
}