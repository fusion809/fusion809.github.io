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
    removeTheta1P1PhasePlot();
    removeTheta1P2PhasePlot();
    removeTheta2P1PhasePlot();
    removeTheta2P2PhasePlot();
    removeP1P2PhasePlot();
    removePendulum1Plot();
    removePendulum1TimePlot();
    removePendulum2Plot();
    removePendulum2TimePlot();
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
function gen2DPlot(x, y, element, title) {
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
        title: title
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
    }
    
    // layout object
    var layoutTimePlot = {
        title: title
    };
    
    // Generate plot
    Plotly.newPlot(element, dataTimePlot, layoutTimePlot);
}