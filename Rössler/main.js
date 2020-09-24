/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param a          Interaction parameter.
 * @param b          Interaction parameter.
 * @param c          Interaction parameter.
 * @param t          Time (seconds).
 * @param x          x value.
 * @param y          y value.
 * @param z          z value.
 * @return           [dx/dt, dy/dt, dz/dt]
 */
function f(a, b, c, t, x, y, z) {
    return [-y-z, x + a*y, b + z*(x-c)];
}

function readInputs() {
    var a = parseFloat(document.getElementById("a").value);
    var b = parseFloat(document.getElementById("b").value);
    var c = parseFloat(document.getElementById("c").value);
    var t0 = parseFloat(document.getElementById("t0").value);
    var tf = parseFloat(document.getElementById("tf").value);
    var x0 = parseFloat(document.getElementById("x0").value);
    var y0 = parseFloat(document.getElementById("y0").value);
    var z0 = parseFloat(document.getElementById("z0").value);
    var epsilon = parseFloat(document.getElementById("epsilon").value);
    var dtInitial = parseFloat(document.getElementById("dtInitial").value);

    var objectOfInputs = {
        a: a,
        b: b,
        c: c,
        t0: t0,
        tf: tf,
        x0: x0,
        y0: y0,
        z0: z0,
        epsilon: epsilon,
        dtInitial: dtInitial
    };
    return objectOfInputs;
}

/** 
 * Solve the problem using RK45.
 *
 * @params           None. Uses parameter values in the forum.
 * @return           Nothing. But it enters the solution values into the solution
 * object.
 */
function solveProblem(objectOfInputs) {
    // Obtain the parameters of the problem
    var {a, b, c, t0, tf, x0, y0, z0, epsilon, dtInitial} = objectOfInputs;

    // Initialize the arrays used and loop variables
    var t = [t0];
    var x = [x0];
    var y = [y0];
    var z = [z0];
    var dt = dtInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in x, y and z
        // over the step
        // 1st approx
        var k1 = dt*f(a, b, c, t[i], x[i], y[i], z[i])[0];
        var l1 = dt*f(a, b, c, t[i], x[i], y[i], z[i])[1];
        var m1 = dt*f(a, b, c, t[i], x[i], y[i], z[i])[2];
        // 2nd approx
        var k2 = dt*f(a, b, c, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4)[0];
        var l2 = dt*f(a, b, c, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4)[1];
        var m2 = dt*f(a, b, c, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4)[2];
        // 3rd approx
        var k3 = dt*f(a, b, c, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32)[0];
        var l3 = dt*f(a, b, c, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32)[1];
        var m3 = dt*f(a, b, c, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32)[2];
        // 4th approx
        var k4 = dt*f(a, b, c, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[0];
        var l4 = dt*f(a, b, c, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[1];
        var m4 = dt*f(a, b, c, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[2];
        // 5th approx
        var k5 = dt*f(a, b, c, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[0];
        var l5 = dt*f(a, b, c, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[1];
        var m5 = dt*f(a, b, c, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[2];
        // 6th approx
        var k6 = dt*f(a, b, c, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[0];
        var l6 = dt*f(a, b, c, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[1];
        var m6 = dt*f(a, b, c, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[2];

        // x1, y1 and z1 are our fourth order approximations
        var x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        var y1 = y[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        var z1 = z[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;

        // x2, y2 and z2 are our fifth order approximations
        var x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        var y2 = y[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        var z2 = z[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;

        // The following are used to correct the step size
        var Rx = Math.abs(x1-x2)/dt;
        var Ry = Math.abs(y1-y2)/dt;
        var Rz = Math.abs(z1-z2)/dt;
        var sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
        var sy = 0.84*Math.pow(epsilon/Ry, 1/4);
        var sz = 0.84*Math.pow(epsilon/Rz, 1/4);
        var R = Math.max(Rx, Ry, Rz);
        var s = Math.min(sx, sy, sz);

        // If R is less than or equal to epsilon move onto the next step
        if ( R <= epsilon ) {
            t.push(t[i]+dt);
            x.push(x1);
            y.push(y1);
            z.push(z1);
            i++;
            dt *= s;
        } else {
            dt *= s;
        }
    }

    // Write t, x, y and z to our solution object
    var solution = {
        t: t,
        x: x,
        y: y,
        z: z
    };
    return solution;
}

/**
 * Tabulates solution data.
 *
 * @params           None. Uses the entries of the solution object, however. 
 * @return           Nothing. Just populates the table with the solution values. 
 */
function fillTable(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(objectOfInputs);

    // Extract coordinate arrays from the solution object
    var {t, x, y, z} = solution;

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    tableContents += '<th>x</th>';
    tableContents += '<th>y</th>';
    tableContents += '<th>z</th>';
    tableContents += "</tr>";
    for (let j = 0; j < x.length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + x[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + y[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '<td>' + z[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        tableContents += '</tr>';
    }
    document.getElementById('tableOutputs').innerHTML = tableContents;
}

/**
 * Removes the solution table
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution table.
 */
function removeTable() {
    // Clear table content
    document.getElementById('tableOutputs').innerHTML = '';
}

/**
 * Generates a 3D phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generate3DPhasePlot(objectOfInputs) {
    // Solve problem
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {x, y, z} = solution;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXYZ").style = "height: " + windowInnerHeight + "px;";

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
        title: 'Phase plot of the solution to the Rössler equations'
    };

    // Generate plot
    Plotly.newPlot('phasePlotXYZ', dataXYZ, layoutXYZ);
}

/**
 * Remove 3D phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function remove3DPhasePlot() {
    document.getElementById("phasePlotXYZ").innerHTML = '';
    document.getElementById("phasePlotXYZ").style = '';
}

/**
 * Generates a XY phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateXYPhasePlot(objectOfInputs) {
    // Solve problem
    var solution = solveProblem();

    // Extract solution data from solution object
    var {x, y} = solution;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXY").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotXY = {
        x: x,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXY = [plotXY];

    // layout object
    var layoutXY = {
        title: "xy phase plot"
    };

    // Generate plot
    Plotly.newPlot('phasePlotXY', dataXY, layoutXY);
}

/**
 * Remove XY phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXYPhasePlot() {
    document.getElementById("phasePlotXY").innerHTML = '';
    document.getElementById("phasePlotXY").style = '';
}

/**
 * Generates a XZ phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateXZPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(objectOfInputs);
    
    // Extract solution data from solution object
    var {x, z} = solution;
    
    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotXZ").style = "height: " + windowInnerHeight + "px;";
    
    // Plot object and data object array
    var plotXZ = {
        x: x,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataXZ = [plotXZ];
    
    // layout object
    var layoutXZ = {
        title: "xz phase plot"
    };
    
    // Generate plot
    Plotly.newPlot('phasePlotXZ', dataXZ, layoutXZ);
}

/**
 * Remove XZ phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeXZPhasePlot() {
    document.getElementById("phasePlotXZ").innerHTML = '';
    document.getElementById("phasePlotXZ").style = '';
}

/**
 * Generates a YZ phase plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateYZPhasePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {y, z} = solution;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("phasePlotYZ").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotYZ = {
        x: y,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };
    var dataYZ = [plotYZ];

    // layout object
    var layoutYZ = {
        title: "yz phase plot"
    };

    // Generate plot
    Plotly.newPlot('phasePlotYZ', dataYZ, layoutYZ);
}

/**
 * Remove YZ phase plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeYZPhasePlot() {
    document.getElementById("phasePlotYZ").innerHTML = '';
    document.getElementById("phasePlotYZ").style = '';
}

/**
 * Generates a time plot
 * 
 * @params           None.
 * @return           Nothing.
 */
function generateTimePlot(objectOfInputs) {
    // Solve the problem
    var solution = solveProblem(objectOfInputs);

    // Extract solution data from solution object
    var {t, x, y, z} = solution;

    // Height and width of plot
    var windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";

    // Plot object and data object array
    var plotTX = {
        x: t,
        y: x,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'x'
    };
    var plotTY = {
        x: t,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'y'
    };
    var plotTZ = {
        x: t,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'z'
    };
    var dataTimePlot = [plotTX, plotTY, plotTZ];

    // layout object
    var layoutTimePlot = {
        title: "Time plots of the solution to the problem"
    };

    // Generate plot
    Plotly.newPlot('timePlot', dataTimePlot, layoutTimePlot);
}

/**
 * Remove time plot
 * 
 * @params           None.
 * @return           Nothing. Just removes the plot.
 */
function removeTimePlot() {
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("timePlot").style = '';
}

/**
 * Generate five plots:
 * - The first is a 3D phase plot of x, y and z.
 * - The second is a 2D phase plot of y against x.
 * - The third is a 2D phase plot of z against x.
 * - The fourth is a 2D phase plot of z against y.
 * - The fifth is a plot of x, y and z against time.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots(objectOfInputs) {
    generate3DPhasePlot(objectOfInputs);
    generateXYPhasePlot(objectOfInputs);
    generateXZPhasePlot(objectOfInputs);
    generateYZPhasePlot(objectOfInputs);
    generateTimePlot(objectOfInputs);
};

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
    // Clear HTML and CSS of the plots
    document.getElementById("timePlot").innerHTML = '';
    document.getElementById("phasePlotXYZ").innerHTML = '';
    document.getElementById("phasePlotXY").innerHTML = '';
    document.getElementById("phasePlotXZ").innerHTML = '';
    document.getElementById("phasePlotYZ").innerHTML = '';
    document.getElementById("timePlot").style = '';
    document.getElementById("phasePlotXYZ").style = '';
    document.getElementById("phasePlotXY").style = '';
    document.getElementById("phasePlotXZ").style = '';
    document.getElementById("phasePlotYZ").style = '';
};