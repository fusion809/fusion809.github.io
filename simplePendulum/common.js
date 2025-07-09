/**
 * Right-hand side of our second-order ODE written as a system of 
 * first-order ODEs.
 *
 * @param g          Acceleration due to gravity in metres per second 
 * squared.
 * @param l          Length of the pendulum in metres.
 * @param t          Time.
 * @param theta      Angle of the pendulum from the positive x-axis.
 * @param thetaDot   Rate of change against time of the angle from the positive x-axis.
 * @return           [dtheta/dt, d2theta/dt2]
 */
function f(g, l, t, theta, thetaDot) {
    var thetaDDot = -g/l * Math.cos(theta);
    return [thetaDot, thetaDDot];
};

/**
 * Returns the value of theta dot squared.
 * 
 * @param g        Acceleration due to gravity.
 * @param l        Length of the pendulum.
 * @param theta0   Initial angle from the positive x-axis.
 * @param thetaDot0 Initial rate of change in the angle from the positive x-axis.
 * @param theta    Angle from the positive x-axis.
 * @return         Theta dot squared value.
 */
function thetaDotSq(g, l, theta0, thetaDot0, theta) {
    return thetaDot0**2 + 2*g/l * (Math.sin(theta0)-Math.sin(theta));
}

/**
 * The derivative with respect to theta of the theta dot squared expression
 * 
 * @param g        Acceleration due to gravity.
 * @param l        Length of the pendulum.
 * @param theta    Angle from the x-axis.
 * @return         Theta derivative of the theta dot squared expression.
 */
function thetaDotSqPrime(g, l, theta) {
    return -2*g/l * Math.cos(theta);
}

/**
 * Correction to theta according to Newton's method
 * 
 * @param g        Acceleration due to gravity.
 * @param l        Length of the pendulum rod.
 * @param theta0   Initial angle from the positive x-axis.
 * @param thetaDot0 Initial rate of change against time of the angle from the positive x-axis.
 * @param theta    Angle from the x-axis.
 * @return         Correction to theta min/max value according to Newton's  method.
 */
function newtonsCorrection(g, l, theta0, thetaDot0, theta) {
    return -thetaDotSq(g, l, theta0, thetaDot0, theta)/thetaDotSqPrime(g, l, theta);
}

/**
 * Calculate the bounds of theta in our integration by using Newton's method
 * on our expression for theta dot squared. 
 * 
 * @param objectOfInputs      An object of problem parameters.
 * @return                    Same array with the bounds for integration on theta added to the end.
 */
function thetaBounds(objectOfInputs) {
    // Parameters of the problem that are necessary to calculate the period
    var {g, l, theta0, thetaDot0, t0} = objectOfInputs;

    // Take our initial guess for thetaMax
    // Check if the problem satisfies the conditions for periodic behaviour
    if ( thetaDot0**2 + 2*g/l * (Math.sin(theta0)-1) > 0) {
        document.getElementById("integralDisplay").innerHTML = "Problem is not periodic.";
        // Too high of tf in these cases can freeze the solver up
        document.getElementById("tf").value = t0 + 10;
        // The default, 1e-11, can freeze the webpage up for these problems
        // e.g. with thetaDot0 > 5, this will likely be the case.
        document.getElementById("epsilon").value = 1e-9;
        // Add extra N to array, so that we can kill periodCalc() if we get 
        // this far
        objectOfInputs.periodic = false;
        // Clearing vars made irrelevant by non-periodicity
        document.getElementById("thetaMinDisplay").innerHTML = '';
        document.getElementById("thetaMaxDisplay").innerHTML = '';
        // Return objectOfInputs
        return objectOfInputs;
    }
    objectOfInputs.periodic = true;
    var thetaMin = Math.asin((thetaDot0**2+2*(g/l)*Math.sin(theta0))/(2*g/l));
    var thetaMax = - thetaMin - Math.PI;

    // Mention theta bounds and Newton's info in data table
    document.getElementById("thetaMinDisplay").innerHTML = thetaMin;
    document.getElementById("thetaMaxDisplay").innerHTML = thetaMax;

    // Add thetaMin and thetaMax to objectOfInputs
    objectOfInputs.thetaMin = thetaMin;
    objectOfInputs.thetaMax = thetaMax;

    return objectOfInputs;
}

/**
 * Calculates the value of theta when theta dot = 0 using Newton's method
 * then uses Chebyshev-Gauss quadrature to compute the time taken to reach this period.
 * 
 * @param solution      An object containing solution data.
 * @return              Nothing. Changes the value of integralDisplay and tf to T and 4T, respectively.
 */
function periodCalc(objectOfInputs) {
    // Parameters of the problem that are necessary to calculate the period
    var {g, l, N, t0, theta0, thetaDot0} = objectOfInputs;
    var nodes = 0;
    var integrand = 0;
    var transformedGrid = 0;
    var integral = 0;
    var period;

    // Kill the function if problem isn't periodic
    if (objectOfInputs.periodic) {
        var thetaMin = objectOfInputs.thetaMin;
        var thetaMax = objectOfInputs.thetaMax;
    } else {
        return;
    }

    // Integrate problem from thetaMin to thetaMax to calculate the period
    // in seconds
    for ( let i = 1; i < N+1; i++) {
        nodes = Math.cos((2*i-1)*Math.PI/(2*N));
        transformedGrid = (thetaMax-thetaMin)*nodes/2+(thetaMax+thetaMin)/2;
        integrand = Math.sqrt(1-nodes**2)*Math.pow(thetaDotSq(g, l, theta0, thetaDot0, transformedGrid),-1/2);
        integral += ((thetaMax-thetaMin)/2) * (Math.PI/N)*integrand;
    }
    period = 2*Math.abs(integral);

    // Change what's displayed on the page accordingly
    document.getElementById("integralDisplay").innerHTML = period;
    if ( period == Infinity ) {
        document.getElementById("tf").value = t0 + 100;
    } else {
        document.getElementById("tf").value = t0 + 4*period;
    }
    return period;
}

/**
 * Generate phase plot of theta dot against theta
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generatePhasePlot(solution) {
    // Extract solution data
    var {theta, thetaDot} = solution;

    // Generate 2D plot
    gen2DPlot(theta, thetaDot, "phasePlot", "Phase plot of theta dot against theta");
}

/**
 * Generate plot of theta and theta dot against time
 * 
 * @param solution       An object containing solution data.
 * @return               Nothing. Just generates the relevant plot.
 */
function generateTimePlot(solution) {
    // Extract solution values
    var {t, theta, thetaDot} = solution;
    var sol = {
        t: t,
        vars: [theta, thetaDot]
    }

    // Generate time plot
    genMultPlot(sol, ["theta", "theta dot"], "timePlot", "Plot of theta and theta dot against time");
}

/**
 * Tabulates solution data.
 *
 * @param objectOfInputs  An object that contains all the problem parameters.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTableSP(objectOfInputs) {
    // Define all global variables
    var solution = solveProblemSP(objectOfInputs);
    var epsilon = objectOfInputs.epsilon;
    var {t, theta, thetaDot} = solution;
    var j = 0;

    // Check the number of elements in solution
    for (m in solution) {
        j++;
    }

    // Add relevant vars if solution contains them
    if (j > 3) {
        var errorThetaDot = solution.errorThetaDot;
        var rmsErrorthetaDot = Math.sqrt(math.dot(errorThetaDot, errorThetaDot)/(errorThetaDot.length-1));
    }
    
    // Generate table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>&theta; (radians) </th>';
    tableContents += '<th>&theta; dot<br/>(radians &middot; s<sup>-1</sup>)</th>';
    if (j > 3) {
        tableContents += '<th>&theta; dot error</th>';
    }
    tableContents += '</tr>';
    for (let i = 0; i < theta.length; i++) {
        tableContents += '<tr>';
        tableContents += '<td>' + i + '</td>';
        tableContents += '<td>' + t[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + theta[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + thetaDot[i].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        if (j > 3) {
            tableContents += '<td>' + errorThetaDot[i].toExponential(5) + '</td>';
        }
        tableContents += '</tr>';
    }
    if (j > 3) {
        tableContents += '<tr>';
        tableContents += '<td colspan = "6">RMS error in &theta; dot is: ';
        tableContents += rmsErrorthetaDot.toExponential(10);
        tableContents += '</td>';
        tableContents += '</tr>';
    }
    document.getElementById('tableOutputs').innerHTML = tableContents;
};

/**
 * Generate two plots:
 * - one of thetaDot and theta against t;
 * - a phase plot of thetaDot against theta.
 * 
 * @param objectOfInputs  An object that contains all the problem parameters.
 * @return               Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    // Solve problem
    var solution = solveProblemSP(objectOfInputs);

    // Generate solution plots
    generateTimePlot(solution);
    generatePhasePlot(solution);
    // The following if statement is to ensure that errorPlot is only generated if it's present in HTML
    if (!!document.getElementById("errorPlot")) {
        generateErrorPlot(solution);
    }
};

function generatePendulumCoords(l, th) {
    var N = th.length;
    var x = new Array(N);
    var y = new Array(N);
    for (let i = 0; i < N; i++) {
        x[i] = l*Math.cos(th[i]);
        y[i] = l*Math.sin(th[i]);
    }
    return [x, y];
}

function animatePendulum(objectOfInputs, solution) {
  const t = solution.t;
  const th = solution.theta;
  var [x, y] = generatePendulumCoords(objectOfInputs.l, th);
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
  const totalLen = objectOfInputs.l;
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
    title: "Simple pendulum",
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
    var solution = solveProblemSP(readInputs());
    animate(solution);
}