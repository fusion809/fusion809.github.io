/**
 * Adjust plot height.
 * 
 * @param element       Element whose height is to be adjusted.
 * @return              None.
 */
function adjustPlotHeight(element) {
    windowInnerHeight = window.innerHeight;
    document.getElementById(element).style = "height: " + windowInnerHeight + "px;";
}

/**
 * Remove plot.
 * 
 * @param element       Plot element to be cleared.
 * @return              None.
 */
function rmPlot(element) {
    if (!!document.getElementById(element)) {
        document.getElementById(element).innerHTML = '';
        document.getElementById(element).style = '';
    }
}

/**
 * Remove XY phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXYPhasePlot() {
    rmPlot("phasePlotXY");
}

function removeR1TPlot() {
    rmPlot("R1TPlot");
}

/**
 * Remove XZ phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXZPhasePlot() {
    rmPlot("phasePlotXZ");
}

/**
 * Remove YZ phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeYZPhasePlot() {
    rmPlot("phasePlotYZ");
}

/**
 * Remove XYZ phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function remove3DPhasePlot() {
    rmPlot("phasePlotXYZ");
}

/**
 * Remove time plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeTimePlot() {
    rmPlot("timePlot");
}

/**
 * Remove error plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeErrorPlot() {
    rmPlot("errorPlot");
}

/**
 * Remove 2D phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removePhasePlot() {
    rmPlot("phasePlot");
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
 * Remove theta1 ptheta2 plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeTheta1P2PhasePlot() {
    rmPlot("phasePlotTheta1P2");
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

/**
 * Remove theta2 ptheta1 plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeTheta2P1PhasePlot() {
    rmPlot("phasePlotTheta2P1");
}

/**
 * Remove xDot vs x plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeXXDotPhasePlot() {
    rmPlot("phasePlotXXDot");
}

/**
 * Remove theta vs x plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeXThetaPhasePlot() {
    rmPlot("phasePlotXTheta");
}

/**
 * Remove thetaDot vs theta plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removeThetaThetaDotPhasePlot() {
    rmPlot("phasePlotThetaThetaDot");
}

/**
 * Remove pendulum 1 coordinates plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removePendulum1Plot() {
    rmPlot("pendulum1Plot");
}

/**
 * Remove pendulum 1 coordinates against time plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removePendulum1TimePlot() {
    rmPlot("pendulum1TimePlot");
}

/**
 * Remove pendulum 2 coordinates plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removePendulum2Plot() {
    rmPlot("pendulum2Plot");
}

/**
 * Remove pendulum 2 coordinates against time plot
 * 
 * @params         None.
 * @return         Nothing.
 */
function removePendulum2TimePlot() {
    rmPlot("pendulum2TimePlot");
}

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    // Clear HTML and CSS of the plots
    removeXYPhasePlot();
    removeXZPhasePlot();
    removeYZPhasePlot();
    remove3DPhasePlot();
    removePhasePlot();
    removeTimePlot();
    removeErrorPlot();
    removeXXDotPhasePlot();
    removeXThetaPhasePlot();
    removeThetaThetaDotPhasePlot();
    removeTheta1Theta2PhasePlot();
    removeTheta1Dtheta1PhasePlot();
    removeTheta1Dtheta2PhasePlot();
    removeTheta2Dtheta1PhasePlot();
    removeTheta2Dtheta2PhasePlot();
    removeDtheta1Dtheta2PhasePlot();
    removePendulum1Plot();
    removePendulum1TimePlot();
    removePendulum2Plot();
    removePendulum2TimePlot();
    removePendulumPlots();
}

/**
 * Remove all pendulum coordinate plots.
 * 
 * @params         None.
 * @return         Nothing.
 */
function removePendulumPlots() {
    removePendulum1Plot();
    removePendulum1TimePlot();
    removePendulum2Plot();
    removePendulum2TimePlot();
    rmPlot("pendulumPlot");
    rmPlot("pendulumTimePlot");
}

/**
 * Generate 2D plot
 * @param x        Array of x-axis values.
 * @param y        Array of y-axis values.
 * @param element  HTML element of the plot.
 * @param title    Title of the plot.
 * @return         Nothing.
 */
function gen2DPlot(x, y, element, title, xtitle, ytitle) {
    // Height and width of the plot
    adjustPlotHeight(element);

    var plotXY = {
        x: x,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    }
    var dataXY = [plotXY];

    // layout object
    var layoutXY = {
        title: title,
    };

    // Generate plot
    Plotly.newPlot(element, dataXY, layoutXY);
}

function gen2DPlotXYLabs(x, y, element, title, xtitle, ytitle) {
    // Height and width of the plot
    adjustPlotHeight(element);

    var plotXY = {
        x: x,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    }
    var dataXY = [plotXY];

    // layout object
    var layoutXY = {
        title: title,
        xaxis: {
            title: {
                text: xtitle
            }
        },
        yaxis: {
            title: {
                text: ytitle
            }
        }
    };

    // Generate plot
    Plotly.newPlot(element, dataXY, layoutXY);
}
/**
 * Generate 3D plot
 * 
 * @param x        Array storing x-values.
 * @param y        Array storing y-values.
 * @param z        Array storing z-values.
 * @param element  HTML element the plot will go in.
 * @param title    Title for the plot.
 */
function gen3DPlot(x, y, z, element, title) {
    // Height and width of plot
    adjustPlotHeight(element);

    // Plot object and data object array
    var plotXYZ = {
        x: x,
        y: y,
        z: z,
        type: 'scatter3d',
        mode: 'lines',
        opacity: 1,
        line: {
           width: 6,
           reversescale: false
        }
    };
    var dataXYZ = [plotXYZ];
   
    // layout object
    var layoutXYZ = {
       title: title
    };
   
    // Generate plot
    Plotly.newPlot(element, dataXYZ, layoutXYZ);
}

/**
 * Generate multiple time plots on the one graph
 * 
 * @param solution Solution object.
 * @param varnames Names for the variables to be displayed on the plot.
 * @param element  HTML element of the plot.
 * @param title    Title of the plot.
 * @return         Nothing.
 */
function genMultPlot(solution, varnames, element, title) {
    // Height and width of plot
    adjustPlotHeight(element);

    // Extract solution
    var {t, vars} = solution;

    // Plot object and data object array
    if (vars.length == 3) {
        var [x, y, z] = vars;

        var plotTX = {
            x: t,
            y: x,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[0]
        };
        var plotTY = {
            x: t,
            y: y,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[1]
        };
        var plotTZ = {
            x: t,
            y: z,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[2]
        };
        var dataTimePlot = [plotTX, plotTY, plotTZ];
    } else if (vars.length == 2) {
        var [x, y] = vars;

        var plotTX = {
            x: t,
            y: x,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[0]
        };
        var plotTY = {
            x: t,
            y: y,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[1]
        };
        var dataTimePlot = [plotTX, plotTY];
    } else if (vars.length == 4) {
        var [S, E, I, R] = vars;

        var plotTS = {
            x: t,
            y: S,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[0]
        };
        var plotTE = {
            x: t,
            y: E,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[1]
        };
        var plotTI = {
            x: t,
            y: I,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[2]
        };
        var plotTR = {
            x: t,
            y: R,
            type: 'scatter',
            mode: 'lines',
            opacity: 1,
            name: varnames[3]
        };
        var dataTimePlot = [plotTS, plotTE, plotTI, plotTR];
    } else {
        var N = vars.length;
        var dataTimePlot = new Array(N);
        for (let i = 0; i < N; i++) {
            var plot = {
                x: t,
                y: vars[i],
                type: 'scatter',
                mode: 'lines',
                opacity: 1,
                name: varnames[i]
            };
            dataTimePlot[i] = plot;
        }
    }
    
    // layout object
    var layoutTimePlot = {
        title: title
    };
    
    // Generate plot
    Plotly.newPlot(element, dataTimePlot, layoutTimePlot);
}

/**
 * Create pendulum animation.
 * @param objectOfInputs Object of page inputs.
 * @param solution       Solution object containing t and vars. 
 * @param label          Label of system, can be "Double pendulum", 
 *                       "Double elastic pendulum", "Elastic pendulum", 
 *                       "Simple pendulum" and "Single pendulum".
 * @param IdSuffix       Suffix for HTML IDs used by animation.
 */
function animatePendulum(objectOfInputs, solution, label, IdSuffix="") {
  const t = solution.t;
  if (label=="Double pendulum") {
    var [t1, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
  } else if (label == "Double elastic pendulum") {
    var [t1, x1, y1, x2, y2] = generatePendulumCoords(solution);
  } else if (label == "Simple pendulum" || label == "Single pendulum" || label == "Elastic pendulum") {
    var [x, y] = generatePendulumCoords(objectOfInputs, solution);
  }
  const padding = 1;
  var data;
  if (label == "Double pendulum" || label == "Double elastic pendulum") {
    var xmin = func2Vecs(Math.min, x1, x2);
    var ymin = func2Vecs(Math.min, y1, y2);
    var xmax = func2Vecs(Math.max, x1, x2);
    var ymax = func2Vecs(Math.max, y1, y2);
    const trace1 = {
        x: [], y: [],
        mode: "lines+markers",
        marker: { size: 8 },
        line: { color: "blue" },
        name: "Rod 1"
    };
    const trace2 = {
        x: [], y: [],
        mode: "lines+markers",
        marker: { size: 8 },
        line: { color: "red" },
        name: "Rod 2"
    };
    data = [trace1, trace2];
    if (isFinite(xmin)) {
        xmin = xmin - padding;
    } else {
        xmin = -10 - padding;
    }
    if (isFinite(xmax)) {
        xmax = xmax + padding;
    } else {
        xmax = 10 + padding;
    }

    if (isFinite(ymin)) {
        ymin = ymin - padding;
    } else {
        ymin = -40 - padding;
    }

    if (isFinite(ymax)) {
        ymax = ymax + padding;
    } else {
        ymax = 10 + padding;
    }
  } else {
    var xmin = Math.min(...x) - padding;
    var ymin = Math.min(...y) - padding;
    var xmax = Math.max(...x) + padding;
    var ymax = Math.max(...y) + padding;
    const trace1 = {
        x: [], y: [],
        mode: "lines+markers",
        marker: { size: 8 },
        line: { color: "blue" },
        name: "Rod"
    };
    data = [trace1];
  }

  
  const layout = {
    title: label,
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
  Plotly.newPlot("animation" + IdSuffix, data, layout).then(() => {

    let startTime = null;
    let frame = 0;
    let paused = false;
    let pauseStart = 0;
    let totalPausedTime = 0;

    const pauseButton = document.getElementById("toggleButton" + IdSuffix);
    pauseButton.addEventListener("click", () => {
      paused = !paused;
      pauseButton.textContent = paused ? "Play" : "Pause";
      if (paused) {
        pauseStart = Date.now();
      } else {
        totalPausedTime += Date.now() - pauseStart;
        if (!paused) requestAnimationFrame(animateFrame);
      }
    });

    const restartButton = document.getElementById("restartButton" + IdSuffix);
    restartButton.addEventListener("click", () => {
      startTime = null;
    });
    const addTimeButton = document.getElementById("addTimeButton" + IdSuffix);
    addTimeButton.addEventListener("click", () => {
      //startTime += objectOfInputs.Time - t[frame];'
      var Deltat = readInputs().Deltat;
      if (t[frame] + Deltat > t[t.length-1]) {
        startTime = null;
      } else {
        startTime += t[frame] - Deltat*1000
        frame = t.findIndex(x => x >= t[frame] + Deltat);
      }
    });

    let cycle = 0;
    function animateFrame(timestamp) {
      if (!startTime || frame == t.length - 1) {
        frame = 0;
        startTime = timestamp; 
      }
      if (paused) return;

      const elapsedSec = (timestamp - startTime - totalPausedTime) / 1000;


    // Advance to the frame corresponding to elapsed time
    while (frame < t.length -1 && t[frame] < elapsedSec) {
      frame++;
    }
    if (frame >= t.length) {
        frame = t.length - 1;
        cycle++;
    }
    if (label == "Double elastic pendulum" || label=="Double pendulum") {
        Plotly.animate("animation" + IdSuffix, {
        data: [
            { x: [0, x1[frame]], y: [0, y1[frame]] },
            { x: [x1[frame], x2[frame]], y: [y1[frame], y2[frame]] }
        ]
        }, {
        transition: { duration: 0 },
        frame: { duration: 0, redraw: true }
        });
    } else {
        Plotly.animate("animation" + IdSuffix, {
            data: [
                { x: [0, x[frame]], y: [0, y[frame]] }
            ]
            }, {
            transition: { duration: 0 },
            frame: { duration: 0, redraw: true }
        });
    }
    layout.annotations[0].text = `Time: ${t[frame].toFixed(2)} s`;
    Plotly.relayout('animation' + IdSuffix, layout);
    requestAnimationFrame(animateFrame);
  }

  requestAnimationFrame(animateFrame);
});
}

function range(x, padding) {
  return [Math.min(...x)-padding, Math.max(...x)+padding]
}
/**
 * Create 2D plot animation.
 * @param solution Solution object containing t and vars.
 * @param varnames Names of the variables to be animated.
 * @param timer    Position of timer.
 * @param IdSuffix The suffix of HTML element IDs for buttons and animation.
 * @param nos      The index of variables to be plotted within vars.
 * @return         Nothing.
 */
function animate2D(solution, varnames=["x", "y"], timer=[0, 0.98], IdSuffix="", nos=[0, 1]) {
  const x = solution.vars[nos[0]];
  const y = solution.vars[nos[1]];
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
    marker: { color: 'red', size: 10 },
    name: 'Object'
  };

  const layout = {
  margin: { l: 0, r: 0, b: 0, t: 0 },  // make space for labels
    xaxis: {
      title: varnames[0],
      range: range(x, padding),
      showticklabels: true,
      automargin: true,
      ticks: 'outside',
      tickfont: { size: 12 }
    },
    yaxis: {
      title: varnames[1],
      range: range(y, padding),
      showticklabels: true,
      automargin: true,
      ticks: 'outside',
      tickfont: { size: 12 }
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

  Plotly.newPlot("animation" + IdSuffix, [tracePath, traceMarker], layout).then(() => {
    let frame = 0;
    let startTime = null;
    let paused = false;
    let pauseStart = 0;
    let totalPausedTime = 0;

    const button = document.getElementById("toggleButton" + IdSuffix);
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

    const restartButton = document.getElementById("restartButton" + IdSuffix);
    restartButton.addEventListener("click", () => {
      startTime = null;
    });

    const addTimeButton = document.getElementById("addTimeButton" + IdSuffix);
    addTimeButton.addEventListener("click", () => {
      //startTime += objectOfInputs.Time - t[frame];'
      var Deltat = readInputs().Deltat;
      if (t[frame] + Deltat > t[t.length-1]) {
        startTime = null;
      } else {
        startTime += t[frame] - Deltat*1000
        frame = t.findIndex(x => x >= t[frame] + Deltat);
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
      Plotly.animate("animation" + IdSuffix, {
        data: [
          tracePath,
          {
            x: [x[frame]],
            y: [y[frame]],
            mode: 'markers',
            type: 'scatter',
            marker: { color: 'red', size: 10 },
            name: 'Object'
          }
        ],
        layout: {
            annotations: [{
                text: `Time: ${t[frame].toFixed(2)} s`,
                xref: 'paper',
                yref: 'paper',
                x: timer[0],
                y: timer[1],
                showarrow: false,
                font: { size: 16 }
            }],
            xaxis: {
                title: varnames[0],
                range: range(x, padding),
                showticklabels: true,
                ticks: 'outside',
                tickfont: { size: 12 }
            },
            yaxis: {
                title: varnames[1],
                range: range(y, padding),
                showticklabels: true,
                ticks: 'outside',
                tickfont: { size: 12 }
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

/**
 * Create a 3D plot animation.
 * @param solution Solution object containing t and vars (within which solution variables are).
 * @param view     Initial view of animation.
 * @param varnames Labels for variables within vars that you wish to plot.
 * @param nos      The index, within vars, of the variables to be plotted.
 * @param IdSuffix Suffix for animation and button HTML IDs. 
 * @return         None
 */
function animate3D(solution, view=[0.5, -2, 0.5], varnames=["x", "y", "z"], nos=[0, 1, 2], IdSuffix="") {
  const x = solution.vars[nos[0]];
  const y = solution.vars[nos[1]];
  const z = solution.vars[nos[2]];
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
    marker: { color: 'red', size: 7 },
    name: 'Object'
  };

  const layout = {
    margin: { l: 0, r: 0, b: 0, t: 0 },
    scene: {
      xaxis: { title: varnames[0] },
      yaxis: { title: varnames[1] },
      zaxis: { title: varnames[2] }
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

  Plotly.newPlot("animation" + IdSuffix, [tracePath, traceMarker], layout).then(() => {
    let frame = 0;
    let startTime = null;
    let paused = false;
    let pauseStart = 0;
    let totalPausedTime = 0;

    const toggleButton = document.getElementById("toggleButton" + IdSuffix);
    toggleButton.addEventListener("click", () => {
      paused = !paused;
      toggleButton.textContent = paused ? "Play" : "Pause";
      if (paused) {
        pauseStart = Date.now();
      } else {
        totalPausedTime += Date.now() - pauseStart;
        if (!paused) requestAnimationFrame(animateFrame);
      }
    });

    const restartButton = document.getElementById("restartButton" + IdSuffix);
    restartButton.addEventListener("click", () => {
      startTime = null;
    });

    const addTimeButton = document.getElementById("addTimeButton" + IdSuffix);
    addTimeButton.addEventListener("click", () => {
      var Deltat = readInputs().Deltat;
      if (t[frame] + Deltat > t[t.length-1]) {
        startTime = null;
      } else {
        startTime += t[frame] - Deltat*1000
        frame = t.findIndex(x => x >= t[frame] + Deltat);
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
      Plotly.animate("animation" + IdSuffix, {
        data: [
          tracePath,
          {
            x: [x[frame]],
            y: [y[frame]],
            z: [z[frame]],
            mode: 'markers',
            type: 'scatter3d',
            marker: { color: 'red', size: 7 },
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
                x: view[0],   // Set to 0 to align the camera with the YZ plane
                y: view[1],
                z: view[2]
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