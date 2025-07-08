/**
 * Right-hand side of our ODE written as a system of first-order ODEs.
 *
 * @param objectOfInputs An object containing problem parameters.
 * @param t              Time (seconds).
 * @param vars           An array of [theta1, p1, theta2, p2]
 * @return               [dtheta1/dt, dp1/dt, dtheta2/dt, dp2/dt]
 */
function f(objectOfInputs, t, vars, dt) {
    var {g, l1, l2, m1b, m1r, m2b, m2r, b1b, b1r, c1b, c1r, b2b, b2r, c2b, c2r} = objectOfInputs;
    var [theta1, dtheta1, theta2, dtheta2] = vars;
    outercoef = 1/((m2r/12+m2b)*l2 - (m2b**2*l2*Math.cos(theta1-theta2)**2)/(m1r/12 + m1b+m2b));
    mass = m1r/2 + m2r + m1b + m2b; 
    v1b = l1*dtheta1;
    v1r = v1b/2;
    v2b = Math.sqrt(l1**2*dtheta1**2+l2**2*dtheta2**2+2*l1*l2*dtheta1*dtheta2*Math.cos(theta1-theta2));
    v2r = Math.sqrt(l1**2*dtheta1**2+(l2**2*dtheta2**2)/4 + l2*dtheta1*dtheta2*Math.cos(theta1-theta2));
    drag1b = (b1b + c1b*v1b)*v1b;
    drag1r = (b1r + c1r*v1r)*v1r/2;
    drag2b = (b2b + c2b*v2b)*(l1*dtheta1 + l2*dtheta2*Math.cos(theta1-theta2));
    drag2b2 = (b2b + c2b*v2b)*(l1*dtheta1*Math.cos(theta1-theta2) + l2*dtheta2);
    drag2r = (b2r + c2r*v2r)*(l1*dtheta1 + l2*dtheta2*Math.cos(theta1-theta2)/2);
    drag2l2 = 1/4*(b2r + c2r*v2r)*(2*l1*dtheta1*Math.cos(theta1-theta2) + l2*dtheta2);
    innercoef = -m2b*Math.cos(theta1-theta2)/(m1r/12+m1b+m2b);
    inner = innercoef*(-m2b*l2*dtheta2**2*Math.sin(theta1-theta2) - g*Math.cos(theta1)*(m1r/2+m2r+m1b+m2b)-drag1b-drag2b-drag1r-drag2r);
    extra = m2b*(l1*dtheta1**2*Math.sin(theta1-theta2)-g*Math.cos(theta2))-drag2l2-drag2b2;
    d2theta2 = outercoef*(inner+extra);
    outercoef1 = 1/((m1r/12+m1b+m2b)*l1);
    innel11 = -m2b*l2*(d2theta2*Math.cos(theta1-theta2)+dtheta2**2*Math.sin(theta1-theta2));
    innel12 = -g*Math.cos(theta1)*mass - drag1b - drag2b - drag1r - drag2r;
    d2theta1 = outercoef1*(innel11 + innel12);

    
    // Return statement
    return [dt*dtheta1, dt*d2theta1, dt*dtheta2, dt*d2theta2];
}

/** 
 * Solve the problem using RKF45
 *
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               [t, vars]
 */
function RKF45(objectOfInputs) {
    // Extract initial from object and add to 2d array
    var {theta10, dtheta10, theta20, dtheta20} = objectOfInputs;
    var vars0 = [[theta10, dtheta10, theta20, dtheta20]];
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
    gen2DPlot(theta1, theta2, "phasePlotTheta1Theta2", "Phase plot of θ<sub>2</sub> against θ<sub>1</sub>");
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
    gen2DPlot(dtheta1, dtheta2, "phasePlotDtheta1Dtheta2", "Phase plot of dθ<sub>2</sub>/dt against dθ<sub>1</sub>/dt");
}

/**
 * Generates a dtheta1 against theta1 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta1Dtheta1PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta1 = vars[0];
    var dtheta1 = vars[1];
    
    // Generate 2D plot
    gen2DPlot(theta1, dtheta1, "phasePlotTheta1Dtheta1", "Phase plot of dθ<sub>1</sub>/dt against θ<sub>1</sub>");
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
    gen2DPlot(theta1, dtheta2, "phasePlotTheta1Dtheta2", "Phase plot of dθ<sub>2</sub>/dt against θ<sub>1</sub>");
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
    gen2DPlot(theta2, dtheta1, "phasePlotTheta2Dtheta1", "Phase plot of dθ<sub>1</sub>/dt against θ<sub>2</sub>");
}

/**
 * Generates a dtheta2 against theta2 phase plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTheta2Dtheta2PhasePlot(solution) {
    // Extract solution data from solution object
    var {vars} = solution;
    var theta2 = vars[2];
    var dtheta2 = vars[3];
    
    // Generate 2D plot
    gen2DPlot(theta2, dtheta2, "phasePlotTheta2Dtheta2", "Phase plot of dθ<sub>2</sub>/dt against θ<sub>2</sub>");
}

/**
 * Generates a time plot
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing.
 */
function generateTimePlot(solution) {
    // Generate time plot
    genMultPlot(solution, ["θ<sub>1</sub>", "dθ<sub>1</sub>/dt", "θ<sub>2</sub>", "dθ<sub>2</sub>/dt"], "timePlot", "Plot of θ<sub>1</sub>, dθ<sub>1</sub>/dt, θ<sub>2</sub> and dθ<sub>2</sub>/dt against time");
}

/**
 * Generate cartesian coordinates 
 * @param func           Function being used to integrate problem
 * @param objectOfInputs Problem parameters.
 */
function generatePendulumCoords(objectOfInputs, solution) {
    // Extract solution values and pendulum lengths
    var {t, vars} = solution;
    var [theta1, dtheta1, theta2, dtheta2] = vars;
    var {l1, l2} = objectOfInputs;
    var N = theta1.length;

    // Initialize arrays that will store x and y coords
    var x1 = new Array(N);
    var x2 = new Array(N);
    var y1 = new Array(N);
    var y2 = new Array(N);
    for (let i = 0; i < N; i++) {
        x1[i] = l1*Math.cos(theta1[i]);
        y1[i] = l1*Math.sin(theta1[i]);
        x2[i] = x1[i] + l2*Math.cos(theta2[i]);
        y2[i] = y1[i] + l2*Math.sin(theta2[i]);
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
function generatePendulumPlots(objectOfInputs, solution) {
    var [t, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
    adjustPlotHeight("pendulumPlot");
    adjustPlotHeight("pendulumTimePlot");
    
    // Show two pendulum bob locations on the same plot
    var plotPen1 = {
        x: x1,
        y: y1,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: "Pendulum 1 bob"
    }
    var plotPen2 = {
        x: x2,
        y: y2,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: "Pendulum 2 bob"
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
        name: "Pendulum 1 bob"
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
        name: 'Pendulum 2 bob'
    };
    var dataPen = [plotPen1Time, plotPen2Time];
    var layoutPenTime = {
        title: "Pendulum bob position against time plot"
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
    generateTheta1Dtheta1PhasePlot(solution);
    generateTheta1Dtheta2PhasePlot(solution);
    generateTheta2Dtheta1PhasePlot(solution);
    generateTheta2Dtheta2PhasePlot(solution);
    generateDtheta1Dtheta2PhasePlot(solution);
    generatePendulumPlots(objectOfInputs, solution)
    generateTimePlot(solution);
}

function func2Vecs(f, vec1, vec2) {
    var val;
    try {
        val = f(...vec1.map((v,i) => f(v, vec2[i])));
    } catch (e) {
        if (e instanceof RangeError) {
            val = vec1.concat(vec2).reduce((a, b) => f(a, b), Infinity);
        } else {
            throw e; // rethrow if it's not the expected error
        }
    }
    return val;
}
function animatePendulum(objectOfInputs, solution) {
  const t = solution.t;
  const vars = solution.vars;
  var [t1, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
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
  const data = [trace1, trace2];
  const padding = 1;
  const totalLen = objectOfInputs.l1 + objectOfInputs.l2;
    if (isFinite(xmin)) {
    xmin = xmin - padding;
  } else {
    xmin = -totalLen - padding;
  }
  if (isFinite(xmax)) {
    xmax = xmax + padding;
  } else {
    xmax = totalLen + padding;
  }

  if (isFinite(ymin)) {
    ymin = ymin - padding;
  } else {
    ymin = -totalLen - padding;
  }

  if (isFinite(ymax)) {
    ymax = ymax + padding;
  } else {
    ymax = totalLen + padding;
  }
  const layout = {
    title: "Double Pendulum",
    xaxis: { range: [xmin-padding, xmax+padding], title: "x" },
    yaxis: { range: [ymin-padding, ymax+padding], title: "y", scaleanchor: "x" },
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
  Plotly.newPlot("animation", data, layout);

  let startTime = null;
  let frame = 0;
  let cycle;
  function animateFrame(timestamp) {
    if (frame == t.length - 1 || !startTime) {
        frame = 0;
        startTime = timestamp; 
    }

    const elapsedSec = (timestamp - startTime) / 1000;

    // Advance to the frame corresponding to elapsed time
    while (frame < t.length -1 && t[frame] < elapsedSec) {
      frame++;
    }
    if (frame >= t.length) {
        frame = t.length - 1;
        cycle++;
    }

    const r1 = objectOfInputs.l1;
    const r2 = objectOfInputs.l2;
    const th1 = vars[0][frame];
    const th2 = vars[2][frame];

    const x1 = r1 * Math.cos(th1);
    const y1 = r1 * Math.sin(th1);
    const x2 = x1 + r2 * Math.cos(th2);
    const y2 = y1 + r2 * Math.sin(th2);

    Plotly.animate("animation", {
      data: [
        { x: [0, x1], y: [0, y1] },
        { x: [x1, x2], y: [y1, y2] }
      ]
    }, {
      transition: { duration: 0 },
      frame: { duration: 0, redraw: true }
    });
    layout.annotations[0].text = `Time: ${t[frame].toFixed(2)} s`;
    Plotly.relayout('animation', layout);

    // if (frame == t.length - 1) {
    //     layout.annotations[0].text = `Animation finished. 3s delay before animation is repeated.`;
    //     setTimeout(() => {
    //         requestAnimationFrame(animateFrame);  // resume animation
    //     }, 3000);
    // } else {
        requestAnimationFrame(animateFrame);
    // }
  }

  requestAnimationFrame(animateFrame);
}

function removeAnimation() {
    rmPlot("animation");
}

function removeTheta1Dtheta1PhasePlot() {
  rmPlot("phasePlotTheta1Dtheta1");
}

function removeTheta1Dtheta2PhasePlot() {
  rmPlot("phasePlotTheta1Dtheta2");
}

function removeTheta2Dtheta1PhasePlot() {
  rmPlot("phasePlotTheta2Dtheta1");
}

function removeTheta2Dtheta2PhasePlot() {
  rmPlot("phasePlotTheta2Dtheta2");
}

function removeDtheta1Dtheta2PhasePlot() {
  rmPlot("phasePlotDtheta1Dtheta2");
}