/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param alpha      Interaction parameter.
 * @param beta       Interaction parameter.
 * @param gamma      Interaction parameter.
 * @param delta      Interaction parameter.
 * @param t          Time (seconds).
 * @param x          Prey population.
 * @param y          Predator population.
 * @return           [dx/dt, dy/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {alpha, beta, gamma, delta} = objectOfInputs;
    var [x, y] = vars;
    return [dt*(alpha*x - beta*x*y), dt*(delta*x*y-gamma*y)];
}

/** 
 * Solve the problem using RKF45
 *
 * @param solution       An object containing solution data.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial conditions from object and write to 2d array
    var {x0, y0} = objectOfInputs;
    var vars0 = [[x0, y0]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}

/**
 * Generate phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the phase plot.
 */
function generatePhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var [x, y] = vars;

    // Generate 2D phase plot
    gen2DPlot(x, y, "phasePlot", "Phase plot of y against x")
}

/**
 * Generate a plot of x and y against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates a plot of x and y against time.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["x", "y"], "timePlot", "Plot of x and y against time");
}

/**
 * Generate two plots:
 * - one of y and x against t; and
 * - a phase plot of y against x.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(RKF45, objectOfInputs);

    // Generate plots
    generatePhasePlot(solution);
    generateTimePlot(solution);
}

function animate(solution) {
  const x = solution.vars[0];
  const y = solution.vars[1];
  const t = solution.t;
  const padding = 1;
  const tracePath = {
    x: x,
    y: y,
    mode: 'lines',
    type: 'scatter',
    line: { color: 'blue', width: 4 },
    name: 'Path',
  };

  const traceMarker = {
    x: [x[0]],
    y: [y[0]],
    mode: 'markers',
    type: 'scatter',
    marker: { color: 'red', size: 5 },
    name: 'Object'
  };

  const layout = {
    margin: { l: 0, r: 0, b: 0, t: 0 },
    xaxis: {
        title: 'x',
        range: [Math.min(x) - padding, Math.max(x) + padding],
        showticklabels: true,   // Ensure tick labels are shown
        ticks: 'outside',       // Optional: put ticks outside the plot
        tickfont: { size: 12 }  // Optional: control label font size
    },
    yaxis: {
        title: 'y',
        range: [Math.min(y) - padding, Math.max(y) + padding],
        showticklabels: true,
        ticks: 'outside',
        tickfont: { size: 12 },
        scaleanchor: "x"
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
                x: 0.5,
                y: 0.98,
                showarrow: false,
                font: { size: 16 }
            }],
            xaxis: {
                title: 'x',
                showticklabels: true,
                ticks: 'outside',
                tickfont: { size: 12 }
            },
            yaxis: {
                title: 'y',
                showticklabels: true,
                ticks: 'outside',
                tickfont: { size: 12 },
                scaleanchor: "x"
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