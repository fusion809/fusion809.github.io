/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param objectOfInputs Interaction parameter.
 * @param t              Time (seconds).
 * @param vars           Solution variables.
 * @return               [dx/dt, dy/dt, dz/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {sigma, rho, beta} = objectOfInputs;
    var [x, y, z] = vars;
    return [dt*sigma*(y-x), dt*(x*(rho-z) - y), dt*(x*y-beta*z)];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions from object and write to 2d array
    var {x0, y0, z0} = objectOfInputs;
    var vars0 = [[x0, y0, z0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generates a 3D phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generate3DPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y, z] = vars;

    gen3DPlot(x, y, z, "phasePlotXYZ", "Lorenz attractor phase plot")
}

/**
 * Generates a XY phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateXYPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y] = vars;

    // Generate 2D plot
    gen2DPlot(x, y, "phasePlotXY", "y against x phase plot");
}

/**
 * Generates a XZ phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateXZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var x = vars[0];
    var z = vars[2];
    
    // Generate 2D plot
    gen2DPlot(x, z, "phasePlotXZ", "z against x phase plot");
}

/**
 * Generates a YZ phase plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateYZPhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var y = vars[1];
    var z = vars[2];

    // Generate 2D plot
    gen2DPlot(y, z, "phasePlotYZ", "z against y phase plot");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing all solution values.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "y", "z"], "timePlot", "Plot of x, y and z against time")
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
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generate3DPhasePlot(solution);
    generateXYPhasePlot(solution);
    generateXZPhasePlot(solution);
    generateYZPhasePlot(solution);
    generateTimePlot(solution);
}

function animate(solution) {
  const x = solution.vars[0];
  const y = solution.vars[1];
  const z = solution.vars[2];
  const t = solution.t;

  const tracePath = {
    x: x,
    y: y,
    z: z,
    mode: 'lines',
    type: 'scatter3d',
    line: { color: 'blue', width: 4 },
    name: 'Path',
  };

  const traceMarker = {
    x: [x[0]],
    y: [y[0]],
    z: [z[0]],
    mode: 'markers',
    type: 'scatter3d',
    marker: { color: 'red', size: 5 },
    name: 'Object'
  };

  const layout = {
    margin: { l: 0, r: 0, b: 0, t: 0 },
    scene: {
      xaxis: { title: 'X' },
      yaxis: { title: 'Y' },
      zaxis: { title: 'Z' }
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
            y: [y[frame]],
            z: [z[frame]],
            mode: 'markers',
            type: 'scatter3d',
            marker: { color: 'red', size: 5 },
            name: 'Object'
          }
        ],
        layout: {
          annotations: [{
            text: `Time: ${t[frame].toFixed(2)} s`,
            xref: 'paper',
            yref: 'paper',
            x: 0.05,
            y: 0.95,
            showarrow: false,
            font: { size: 16 }
          }],
          scene: {
            xaxis: { title: 'X' },
            yaxis: { title: 'Y' },
            zaxis: { title: 'Z' },
            camera: {
                eye: {
                x: 2,   // Set to 0 to align the camera with the YZ plane
                y: 0,
                z: 0
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