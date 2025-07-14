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

function min(x) {
  var min;
  try {
    min = Math.min(...x)
  } catch(e) {
    if (e instanceof RangeError) {
      min = x.reduce((a, b) => Math.min(a, b), Infinity);
    } else {
      throw e;
    }
  }
  return min;
}
function max(x) {
  var max;
  try {
    max = Math.max(...x)
  } catch(e) {
    if (e instanceof RangeError) {
      max = x.reduce((a, b) => Math.max(a, b), -Infinity);
    } else {
      throw e;
    }
  }
  return max;
}
/**
 * Range of values for a plot. 
 * @param x Either an array or array of arrays.
 * @returns Vector of values corresponding to minimum of x - 1% padding, 
 * maximum of x + 1% padding.
 */
function range(x) {
  if (x[0].length == undefined) {
    var xmin = min(x);
    var xmax = max(x);
  } else {
    var xmin = max(x[0]);
    var xmax = min(x[0]);
    for (let i = 0; i < x.length; i++) {
      let xmini = min(x[i]);
      let xmaxi = max(x[i]);
      if (xmin > xmini) {
        xmin = xmini;
      } 
      if (xmax < xmaxi) {
        xmax = xmaxi;
      }
    }
  }
  let range = xmax - xmin;
  let padding = range * 0.02;
  return [xmin - padding, xmax+padding]
}


function addTitleAndFrmt(objectOfInputs, element, {title: title, label:label} = {}) {
  if (/[pP]lot/.test(element)) {
    var plotID = element;
  } else {
    var plotID = "animation" + element;
  }
  var contID = "container" + element;
  var contEl = document.getElementById(contID)
  var plotEl = document.getElementById(plotID)
  if (label == undefined) {
    var label = title;
  }
  if (title == undefined) {
    var title = "Animation of the " + label.toLowerCase() + ".";
  }

  var infoEl = document.getElementById("info" + element)
  if (infoEl) {
    infoEl.innerHTML = title;
  }


  infoEl.style = "padding-top: 0px; border-bottom: 1px solid black; margin-bottom: 5px; width:100%;"
  contEl.style = "border: 1px solid black; padding-bottom: 5px; padding-top: 0px; padding-left: 5px; width:100%";
  plotEl.width = objectOfInputs.Width + "px";
  plotEl.height = objectOfInputs.Height + "px";
  return plotID;
}

/**
 * Calculate scaled font size based on plot dimensions
 * @param {number} width - Plot width
 * @param {number} height - Plot height
 * @param {number} minSize - Minimum font size (default: 16)
 * @param {number} maxSize - Maximum font size (default: 24)
 * @param {number} baseWidth - Base width for minimum size (default: 800)
 * @param {number} baseHeight - Base height for minimum size (default: 600)
 * @param {number} maxWidth - Width for maximum size (default: 1920)
 * @param {number} maxHeight - Height for maximum size (default: 1080)
 * @returns {number} Scaled font size
 */
function getScaledFontSize(width, height, minSize = 16, maxSize = 24, 
                          baseWidth = 800, baseHeight = 600, 
                          maxWidth = 1920, maxHeight = 1080) {
    // Calculate scaling factors for width and height
    const widthFactor = Math.min(Math.max((width - baseWidth) / (maxWidth - baseWidth), 0), 1);
    const heightFactor = Math.min(Math.max((height - baseHeight) / (maxHeight - baseHeight), 0), 1);
    
    // Use the average of width and height factors (or you could use max/min)
    const scaleFactor = (widthFactor + heightFactor) / 2;
    
    // Linear interpolation between minSize and maxSize
    const fontSize = minSize + (maxSize - minSize) * scaleFactor;
    
    return Math.round(fontSize);
}

/**
 * Get scaled font sizes for different plot elements
 * @param {number} width - Plot width
 * @param {number} height - Plot height
 * @returns {object} Object with scaled font sizes for different elements
 */
function getPlotFontSizes(objectOfInputs) {
    return {
        title: getScaledFontSize(objectOfInputs.Width, objectOfInputs.Height, 20, 32),      // Title: 20-32px
        axisTitle: getScaledFontSize(objectOfInputs.Width, objectOfInputs.Height, 16, 24),  // Axis titles: 16-24px
        axisLabels: getScaledFontSize(objectOfInputs.Width, objectOfInputs.Height, 12, 18), // Axis labels: 12-18px
        legend: getScaledFontSize(objectOfInputs.Width, objectOfInputs.Height, 12, 16),     // Legend: 12-16px
        annotation: getScaledFontSize(objectOfInputs.Width, objectOfInputs.Height, 16, 24)  // Annotations: 16-24px
    };
}

function margin(objectOfInputs) {
  return { l: Math.max(0.04*objectOfInputs.Width, 60), 
    r: 0, 
    b: 0.20*objectOfInputs.Height, 
    t: 0.20*objectOfInputs.Height};
}
/**
 * Generate 2D plot
 * @param x        Array of x-axis values.
 * @param y        Array of y-axis values.
 * @param element  HTML element of the plot.
 * @param title    Title of the plot.
 * @return         Nothing.
 */
function gen2DPlot(x, y, element, title, xtitle="x", ytitle="y") {
  // Height and width of the plot
  adjustPlotHeight(element);
  var objectOfInputs = readInputs();

  var plotXY = {
    x: x,
    y: y,
    type: 'scatter',
    mode: 'lines',
    opacity: 1
  }
  var dataXY = [plotXY];

  // layout object
  var fontSizes = getPlotFontSizes(objectOfInputs)
  var layoutXY = {
    margin: margin(objectOfInputs),  // make space for labels
    width: objectOfInputs.Width,
    height: objectOfInputs.Height,
    title: {text: title, font: {size: fontSizes.title }},
    xaxis: {
      range: range(x),
      title: {text: xtitle, font: {size: fontSizes.axisTitle}},
      tickfont: {size: fontSizes.axisLabels}
    },
    yaxis: {
      range: range(y),
      title: {text: ytitle, font: {size: fontSizes.axisTitle}},
      tickfont: {size: fontSizes.axisLabels}
    },
    font: {size: fontSizes.legend}
  };
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
function gen3DPlot(x, y, z, element, title, {view = undefined, xtitle="x", ytitle="y", ztitle="z"} = {}) {
    // Height and width of plot
    adjustPlotHeight(element);
    
    var objectOfInputs = readInputs();

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
    var fontSizes = getPlotFontSizes(objectOfInputs)
    var layoutXYZ = {
      margin: margin(objectOfInputs),  // make space for labels
      width: objectOfInputs.Width,
      height: objectOfInputs.Height,
      title: {text: title, font: {size: fontSizes.title }},
      xaxis: {
        range: range(x),
        title: {text: xtitle, font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      yaxis: {
        range: range(y),
        title: {text: ytitle, font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      font: {size: fontSizes.legend},
      zaxis: {
        range: range(z),
        title: { text: ztitle, font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      }
    }
  if (typeof view !== "undefined") {
    layoutXYZ.scene = {
      camera: {
        eye: {
          x: view[0],   // Set to 0 to align the camera with the YZ plane
          y: view[1],
          z: view[2]
          }
        }
       }
  }
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
    var objectOfInputs = readInputs();

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
    
    var fontSizes = getPlotFontSizes(objectOfInputs);
    // layout object
    var layoutTimePlot = {
      margin: margin(objectOfInputs),  // make space for labels
      width: objectOfInputs.Width,
      height: objectOfInputs.Height,
      title: {text: title, font: {size: fontSizes.title }},
      xaxis: {
        range: range(t),
        title: {text: "t", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      yaxis: {
        range: range(vars),
        title: {text: "Variables", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      font: {size: fontSizes.legend}
    };
  Plotly.newPlot(element, dataTimePlot, layoutTimePlot);
}

function generatePendulumCoords(objectOfInputs, solution) {
  if (objectOfInputs.hasOwnProperty("l1") && solution.vars.length == 8) {
    // Extract solution values and pendulum lengths
    var {t, vars} = solution;
    var [r1, dr1, r2, dr2, theta1, dtheta1, theta2, dtheta2] = vars;
    var N = theta1.length;

    // Initialize arrays that will store x and y coords
    var x1 = new Array(N);
    var x2 = new Array(N);
    var y1 = new Array(N);
    var y2 = new Array(N);
    for (let i = 0; i < N; i++) {
        x1[i] = r1[i]*Math.cos(theta1[i]);
        y1[i] = r1[i]*Math.sin(theta1[i]);
        x2[i] = x1[i] + r2[i]*Math.cos(theta2[i]);
        y2[i] = y1[i] + r2[i]*Math.sin(theta2[i]);
    }

    // Return t and Cartesian coordinates of the pendulum bobs
    return [t, x1, y1, x2, y2];
  } else if (objectOfInputs.hasOwnProperty("l1")) {
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
  } else {
    var vars = solution.vars;
    var r = vars[0].map(item => item + objectOfInputs.l0);
    var th = vars[2];
    var N = th.length;
    var x = new Array(N);
    var y = new Array(N);
    for (let i = 0; i < N; i++) {
        x[i] = r[i]*Math.cos(th[i]);
        y[i] = r[i]*Math.sin(th[i]);
    }
    return [x, y];
  }
}

function generatePendulumPlot(objectOfInputs, solution) {
  if (objectOfInputs.hasOwnProperty("l1") && solution.vars.length == 4) {
    // This is for double pendulum
    var [t1, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
    var x=[x1, x2];
    var y=[y1, y2];
  } else if (objectOfInputs.hasOwnProperty("l1")) {
    // Double elastic pendulum
    var [t1, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
    var x=[x1, x2];
    var y=[y1, y2];
  } else {
    // single pendulum systems
    var [x, y] = generatePendulumCoords(objectOfInputs, solution);
    i
  }
  var element = "pendulumPlot"
  var fontSizes = getPlotFontSizes(objectOfInputs);
  if (x1 == undefined) {
    adjustPlotHeight("pendulumPlot");
        
    // Show two pendulum bob locations on the same plot
    var plotPen = {
        x: x,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: "Pendulum bob"
    }
    var dataPen = [plotPen];
    var layoutPen = {
      margin: margin(objectOfInputs),  // make space for labels
      width: objectOfInputs.Width,
      height: objectOfInputs.Height,
      title: {text: "Pendulum position plot.", font: {size: fontSizes.title }},
      xaxis: {
        range: range(x),
        title: {text: "x", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      yaxis: {
        range: range(y),
        title: {text: "y", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      font: {size: fontSizes.legend}
    };
  } else {
    adjustPlotHeight(element);
    
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
      margin: margin(objectOfInputs),  // make space for labels
      width: objectOfInputs.Width,
      height: objectOfInputs.Height,
      title: {text: "Pendulum position plot.", font: {size: fontSizes.title }},
      xaxis: {
        range: range(x),
        title: {text: "x", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      yaxis: {
        range: range(y),
        title: {text: "y", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      font: {size: fontSizes.legend}
    };
  }
  Plotly.newPlot(element, dataPen, layoutPen);
}

function generatePendulumTimePlot(objectOfInputs, solution) {
  if (objectOfInputs.hasOwnProperty("l1") && solution.vars.length == 4) {
    // This is for double pendulum
    var [t1, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
    var x=[x1, x2];
    var y=[y1, y2];
  } else if (objectOfInputs.hasOwnProperty("l1")) {
    // Double elastic pendulum
    var [t1, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
    var x=[x1, x2];
    var y=[y1, y2];
  } else {
    // single pendulum systems
    var [x, y] = generatePendulumCoords(objectOfInputs, solution);
    i
  }
  var element = "pendulumTimePlot";
  adjustPlotHeight(element);
  var t = solution.t
  if (x1 == undefined) {
    // Plot pendulum bob location against time plot
    var plotPenTime = {
      x: t,
      y: x,
      z: y,
      type: 'scatter3d',
      mode: 'lines',
      opacity: 1,
      line: {
        width: 6,
        reversescale: false
      },
      name: "Pendulum"
    }
    var dataPen = [plotPenTime];
  } else {
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
    var x = [x1, x2];
    var y = [y1, y2];
    var dataPen = [plotPen1Time, plotPen2Time];
  }
  var fontSizes = getPlotFontSizes(objectOfInputs);
  var layoutPenTime = {
    margin: margin(objectOfInputs),  // make space for labels
    width: objectOfInputs.Width,
    height: objectOfInputs.Height,
    title: {text: "Generate pendulum time plot.", font: {size: fontSizes.title }},
    xaxis: {
      range: range(t),
      title: {text: "t", font: {size: fontSizes.axisTitle}},
      tickfont: {size: fontSizes.axisLabels}
      },
      yaxis: {
        range: range(x),
        title: {text: "x", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      },
      font: {size: fontSizes.legend},
      zaxis: {
        range: range(y),
        title: { text: "y", font: {size: fontSizes.axisTitle}},
        tickfont: {size: fontSizes.axisLabels}
      }
  };
  Plotly.newPlot(element, dataPen, layoutPenTime);
}

function generatePendulumPlots(objectOfInputs, solution) {
  generatePendulumPlot(objectOfInputs, solution);
  generatePendulumTimePlot(objectOfInputs, solution);
}
/**
 * Update time label in animation.
 * @param layout   Layout object for animation.
 * @param t        Time vector for solution.
 * @param frame    Frame number we are up to in the animation.
 * @param IdSuffix Suffix of element IDs used.
 * @return         Nothing.
 */
function updateTimeLabel(layout, t, state, IdSuffix) {
  layout.annotations[0].text = `Time: ${t[state.frame].toFixed(2)} s`;
  Plotly.relayout('animation' + IdSuffix, layout);
}

/**
 * Pause animation via button.
 * @param state        State object for animation.
 * @param IdSuffix     Suffix of element IDs used.
 * @param animateFrame animateFrame function.
 * @return             Nothing.
 */
function pause(state, IdSuffix, animateFrame) {
  const pauseButton = document.getElementById("toggleButton" + IdSuffix);
  pauseButton.addEventListener("click", () => {
    state.paused = !state.paused;
    pauseButton.textContent = state.paused ? "Play" : "Pause";
    if (state.paused) {
      state.pauseStart = Date.now();
    } else {
      state.totalPausedTime += Date.now() - state.pauseStart;
      requestAnimationFrame(animateFrame);
    }
  });
}

/**
 * Restart animation via pressing a button.
 * @param state    State object for animation.
 * @param layout   Layout object for animation.
 * @param t        Time vector for solution being animated.
 * @param IdSuffix Suffix for element IDs used in animation.
 * @return         Nothing.
 */
function restart(state, layout, t, IdSuffix) {
  const restartButton = document.getElementById("restartButton" + IdSuffix);
  restartButton.addEventListener("click", () => {
    state.startTime =  null;
    state.frame =  0;  // Reset state.frame
    updateTimeLabel(layout, t, state, IdSuffix);
  });
}

/**
 * Add time to animation.
 * @param state    State object for animation.
 * @param layout   Layout object for animation.
 * @param t        Time vector for solution being animated.
 * @param IdSuffix Suffix for element IDs used in animation.
 * @return         Nothing.
 */
function addTime(state, layout, t, IdSuffix) {
  const addTimeButton = document.getElementById("addTimeButton" + IdSuffix);
  addTimeButton.addEventListener("click", () => {
    //startTime += objectOfInputs.Time - t[state.frame];'
    var Deltat = readInputs().Deltat;
    if (t[state.frame] + Deltat > t[t.length-1]) {
      state.startTime =  null;
      state.frame = 0;
    } else {
      state.startTime += t[state.frame] - Deltat*1000;
      let targetTime = t[state.frame] + Deltat;
      state.frame = t.reduce((prevIndex, currValue, currIndex, array) => {
        return Math.abs(currValue - targetTime) < Math.abs(array[prevIndex] - targetTime)
          ? currIndex
          : prevIndex;
      });
    }
    updateTimeLabel(layout, t, state, IdSuffix);
  });
}

/**
 * Skip to t1 in animation.
 * @param state    State object of system.
 * @param layout   Layout object of animation.
 * @param t        Time vector for solution.
 * @param IdSuffix Suffix for elements used in animation.
 * @return         Nothing.
 */
function skipTo(state, layout, t, IdSuffix) {
  const skipToButton = document.getElementById("skipToButton" + IdSuffix);
  skipToButton.addEventListener("click", () => {
    var t1 = readInputs().t1;
    if (t1 > t[t.length-1]) {
      state.startTime =  null;
      state.frame =  0;
      alert("t1 > tf, so restarting!");
    } else {
      state.startTime += (t[state.frame]-t1) * 1000;
      state.frame = t.reduce((prevIndex, currValue, currIndex, array) => {
        return Math.abs(currValue - t1) < Math.abs(array[prevIndex] - t1)
          ? currIndex
          : prevIndex;
      });
    }
    updateTimeLabel(layout, t, state, IdSuffix);
  });
}

function buttons(state, layout, t, IdSuffix, animateFrame) {
  pause(state, IdSuffix, animateFrame);
  restart(state, layout, t, IdSuffix);
  addTime(state, layout, t, IdSuffix);
  skipTo(state, layout, t, IdSuffix);
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
    var x=[x1, x2];
    var y=[y1, y2];
  } else if (label == "Double elastic pendulum") {
    var [t1, x1, y1, x2, y2] = generatePendulumCoords(objectOfInputs, solution);
    var x=[x1, x2];
    var y=[y1, y2];
  } else if (label == "Simple pendulum" || label == "Single pendulum" || label == "Elastic pendulum") {
    var [x, y] = generatePendulumCoords(objectOfInputs, solution);
    i
  }
  var data;
  var animID = addTitleAndFrmt(objectOfInputs, IdSuffix, {label: label});
  var animElem = document.getElementById(animID);
  if (animElem.frameRequest !== undefined) {
    cancelAnimationFrame(animElem.frameRequest);
  }
  if (label == "Double pendulum" || label == "Double elastic pendulum") {
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
  } else {
    const trace1 = {
        x: [], y: [],
        mode: "lines+markers",
        marker: { size: 8 },
        line: { color: "blue" },
        name: "Rod"
    };
    data = [trace1];
  }

  var fontSizes = getPlotFontSizes(objectOfInputs);
  const layout = {
    margin: margin(objectOfInputs),  // make space for labels
    width: objectOfInputs.Width,
    height: objectOfInputs.Height,
    xaxis: {
      range: range(x),
      title: {text: "x", font: {size: fontSizes.axisTitle}},
      tickfont: {size: fontSizes.axisLabels}
    },
    yaxis: {
      range: range(y),
      title: {text: "y", font: {size: fontSizes.axisTitle}},
      tickfont: {size: fontSizes.axisLabels}
    },
    font: {size: fontSizes.legend},
    showlegend: false,
    annotations: [{
      x: 0,
      y: 1.1,  // slightly above plot
      xref: 'paper',
      yref: 'paper',
      text: 'Time: 0.00 s',
      automargin: true,
      showarrow: false,
      font: { size: 16 }
    }]
  };
  
  var tScale = objectOfInputs.tScale;
  Plotly.newPlot(animID, data, layout).then(() => {

    let state = {
      startTime: null,
      frame: 0,
      paused: false,
      pauseStart: 0,
      totalPausedTime: 0,
      animating: false
    };
    let cycle = 0;
    function animateFrame(timestamp) {
      if (!state.startTime || state.frame == t.length - 1) {
        state.frame = 0;
        state.startTime = timestamp; 
      }

      if (state.paused) {
        state.animating = false; // Mark animation as stopped
        return;
      }

      const elapsedSec = tScale*(timestamp - state.startTime - state.totalPausedTime) / 1000;


      // Advance to the state.frame corresponding to elapsed time
      while (state.frame < t.length -1 && t[state.frame] < elapsedSec) {
        state.frame++;
      }
      if (state.frame >= t.length) {
          state.frame =  t.length - 1;
          cycle++;
      }
      if (label == "Double elastic pendulum" || label=="Double pendulum") {
          Plotly.animate(animID, {
          data: [
              { x: [0, x1[state.frame]], y: [0, y1[state.frame]] },
              { x: [x1[state.frame], x2[state.frame]], y: [y1[state.frame], y2[state.frame]] }
          ]
          }, {
          transition: { duration: 0 },
          frame: { duration: 0, redraw: true }
          });
      } else {
          Plotly.animate(animID, {
              data: [
                  { x: [0, x[state.frame]], y: [0, y[state.frame]] }
              ]
              }, {
              transition: { duration: 0 },
              frame: { duration: 0, redraw: true }
          });
      }
      layout.annotations[0].text = `Time: ${t[state.frame].toFixed(2)} s`;
      Plotly.relayout(animID, layout);
      animElem.frameRequest = requestAnimationFrame(animateFrame);
    }

    buttons(state, layout, t, IdSuffix, animateFrame);

    animElem.frameRequest = requestAnimationFrame(animateFrame);
  });
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
function animate2D(solution, {varnames=["x", "y"], timer=[0.05, 0.98], IdSuffix="", nos=[0, 1], title="Phase plot animation of y against x."} = {}) {
  const x = solution.vars[nos[0]];
  const y = solution.vars[nos[1]];
  var objectOfInputs = readInputs();
  var animID = addTitleAndFrmt(objectOfInputs, IdSuffix, {title: title});
  var animElem = document.getElementById(animID);
  if (animElem.frameRequest !== undefined) {
    cancelAnimationFrame(animElem.frameRequest);
  }
  const t = solution.t;
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

  var fontSizes = getPlotFontSizes(objectOfInputs);
  const layout = {
    margin: margin(objectOfInputs),  // make space for labels
    width: objectOfInputs.Width,
    height: objectOfInputs.Height,
    xaxis: {
      title: {text: varnames[0], font: {size: fontSizes.axisTitle}},
      range: range(x),
      showticklabels: true,
      automargin: true,
      ticks: 'outside',
      tickfont: { size: fontSizes.axisLabels }
    },
    yaxis: {
      title: {text: varnames[1], font: {size: fontSizes.axisTitle}},
      range: range(y),
      showticklabels: true,
      automargin: true,
      ticks: 'outside',
      tickfont: { size: fontSizes.axisLabels }
    },
    font: {size: fontSizes.legend},
    annotations: [{
      text: `Time: 0.00 s`,
      xref: 'paper',
      yref: 'paper',
      x: 0.05,
      y: 0.95,
      showarrow: false,
      textangle: 0,
      font: { size: 16 }
    }],
    showlegend: false
  };

  var tScale = objectOfInputs.tScale;
  Plotly.newPlot(animID, [tracePath, traceMarker], layout).then(() => {
    let state = {
      frame: 0,
      startTime: null,
      paused: false,
      pauseStart: 0,
      totalPausedTime: 0
    };

    let cycle = 0;
    function animateFrame(timestamp) {
      if (state.frame == t.length - 1 || !state.startTime) {
        state.frame =  0;
        state.startTime =  timestamp; 
      }
      if (state.paused) return;

      const elapsedSec = tScale*(timestamp - state.startTime - state.totalPausedTime) / 1000;

      while (state.frame < t.length - 1 && t[state.frame] < elapsedSec) {
        state.frame++;
      }
      if (state.frame >= t.length) {
        state.frame =  t.length - 1;
        cycle++;
      }
      Plotly.animate(animID, {
        data: [
          tracePath,
          {
            x: [x[state.frame]],
            y: [y[state.frame]],
            mode: 'markers',
            type: 'scatter',
            marker: { color: 'red', size: 10 },
            name: 'Object'
          }
        ],
        layout: {
            annotations: [
              {
                text: `Time: ${t[state.frame].toFixed(2)} s`,
                xref: 'paper',
                yref: 'paper',
                x: timer[0],
                y: timer[1],
                automargin: true,
                showarrow: false,
                font: { size: 16 }
              }
            ],
            margin: margin(objectOfInputs),
            xaxis: {
                title: {text: varnames[0], font: {size: fontSizes.axisTitle}},
                range: range(x),
                showticklabels: true,
                ticks: 'outside',
                tickfont: { size: fontSizes.axisLabels }
            },
            yaxis: {
                title: {text: varnames[1], font: {size: fontSizes.axisTitle}},
                range: range(y),
                showticklabels: true,
                ticks: 'outside',
                tickfont: { size: fontSizes.axisLabels }
            },
            font: {size: fontSizes.legend}
        }
      }, {
        transition: { duration: 0 },
        frame: { duration: 0, redraw: true }
      });
      
      animElem.frameRequest = requestAnimationFrame(animateFrame);
    }
    buttons(state, layout, t, IdSuffix, animateFrame);
    animElem.frameRequest = requestAnimationFrame(animateFrame);
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
function animate3D(solution, {view = [0.5, -2, 0.5], 
  varnames = ["x", "y", "z"],
  nos = [0, 1, 2], 
  IdSuffix = "",
  title = "X, Y and Z phase plot."} = {}) {
  const x = solution.vars[nos[0]];
  const y = solution.vars[nos[1]];
  const z = solution.vars[nos[2]];
  const t = solution.t;
  var objectOfInputs = readInputs();
  var animID = addTitleAndFrmt(objectOfInputs, IdSuffix, {title: title});
  var animElem = document.getElementById(animID);
  if (animElem.frameRequest !== undefined) {
    cancelAnimationFrame(animElem.frameRequest);
  }

  const tracePath = {
    x: x,
    y: y,
    z: z,
    mode: 'lines',
    type: 'scatter3d',
    opacity: objectOfInputs.Opacity,
    line: { color: 'blue', width: 4 },
    name: 'Path',
  };

  const traceMarker = {
    x: [x[0]],
    y: [y[0]],
    z: [z[0]],
    mode: 'markers',
    type: 'scatter3d',
    opacity: 1,
    marker: { color: 'red', size: 7, },
    name: 'Object'
  };

  var fontSizes = getPlotFontSizes(objectOfInputs);
  const layout = {
    margin: margin(objectOfInputs),
    width: objectOfInputs.Width,
    height: objectOfInputs.Height,
    scene: {
      xaxis: { title: {text: varnames[0], font: {size: fontSizes.axisTitle}}, tickfont: fontSizes.axisLabels },
      yaxis: { title: {text: varnames[1], font: {size: fontSizes.axisTitle}}, tickfont: fontSizes.axisLabels },
      zaxis: { title: {text: varnames[2], font: {size: fontSizes.axisTitle}}, tickfont: fontSizes.axisLabels },
      font: {size: fontSizes.legend},
      camera: {
                eye: {
                x: view[0],   // Set to 0 to align the camera with the YZ plane
                y: view[1],
                z: view[2]
                }
            }
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
    showlegend: false
  };

  var tScale = objectOfInputs.tScale;
  Plotly.newPlot(animID, [tracePath, traceMarker], layout).then(() => {
    let state = {
      frame: 0,
      startTime:  null,
      paused: false,
      pauseStart: 0,
      totalPausedTime: 0
    };

    let cycle = 0;
    function animateFrame(timestamp) {
      if (state.frame == t.length - 1 || !state.startTime) {
        state.frame =  0;
        state.startTime =  timestamp; 
      }
      if (state.paused) return;

      const elapsedSec = tScale*(timestamp - state.startTime - state.totalPausedTime) / 1000;

      while (state.frame < t.length - 1 && t[state.frame] < elapsedSec) {
        state.frame++;
      }
      if (state.frame >= t.length) {
        state.frame =  t.length - 1;
        cycle++;
      }
      Plotly.animate(animID, {
        data: [
          tracePath,
          {
            x: [x[state.frame]],
            y: [y[state.frame]],
            z: [z[state.frame]],
            mode: 'markers',
            type: 'scatter3d',
            opacity: 1,
            marker: { color: 'red', size: 7, opacity: 1, symbol: 'circle', layer: 'above' },
            name: 'Object'
          }
        ],
        layout: {
          annotations: [{
            text: `Time: ${t[state.frame].toFixed(2)} s`,
            xref: 'paper',
            yref: 'paper',
            x: 0.05,
            y: 0.95,
            //automargin: true,
            showarrow: false,
            font: { size: 16 }
          }],
          scene: {
            xaxis: { title: {text: varnames[0], font: {size: 16}} },
            yaxis: { title: {text: varnames[1], font: {size: 16}} },
            zaxis: { title: {text: varnames[2], font: {size: 16}} }
          }
          }
        }, {
          transition: { duration: 0 },
          frame: { duration: 0, redraw: true }
        });
        animElem.frameRequest  = requestAnimationFrame(animateFrame);
      }
    buttons(state, layout, t, IdSuffix, animateFrame);
    animElem.frameRequest = requestAnimationFrame(animateFrame);
  });
}