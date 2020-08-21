/**
 * Right-hand side of our second-order ODE written as a simple of first-order
 * ODEs.
 *
 * @param sigma      Interaction parameter.
 * @param rho        Interaction parameter.
 * @param beta       Interaction parameter.
 * @param t          Time (seconds).
 * @param x          x value.
 * @param y          y value.
 * @param z          z value.
 * @return           [dx/dt, dy/dt, dz/dt]
 */
function f(sigma, rho, beta, t, x, y, z) {
    return [sigma*(y-x), x*(rho-z) - y, x*y-beta*z];
}

// Initialize our global variables
var solution = {
    t: [],
    x: [],
    y: [],
    z: []
};
var epsilon;

/** 
 * Solve the problem using RK45.
 *
 * @params           None. Uses parameter values in the forum.
 * @return           Nothing. But it enters the solution values into the solution
 * object.
 */
function solveProblem() {
    // Obtain the parameters of the problem
    sigma = parseFloat(document.getElementById("sigma").value);
    rho = parseFloat(document.getElementById("rho").value);
    beta = parseFloat(document.getElementById("beta").value);
    t0 = parseFloat(document.getElementById("t0").value);
    tf = parseFloat(document.getElementById("tf").value);
    x0 = parseFloat(document.getElementById("x0").value);
    y0 = parseFloat(document.getElementById("y0").value);
    z0 = parseFloat(document.getElementById("z0").value);
    epsilon = parseFloat(document.getElementById("epsilon").value);
    dtInitial = parseFloat(document.getElementById("dtInitial").value);

    // Initialize the arrays used and loop variables
    t = [t0];
    x = [x0];
    y = [y0];
    z = [z0];
    dt = dtInitial;
    i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        // Step size, as dictated by the method
        dt = Math.min(dt, tf-t[i]);

        // Runge-Kutta-Fehlberg approximations of the change in x and y
        // over the step
        k1 = dt*f(sigma, rho, beta, t[i], x[i], y[i], z[i])[0];
        l1 = dt*f(sigma, rho, beta, t[i], x[i], y[i], z[i])[1];
        m1 = dt*f(sigma, rho, beta, t[i], x[i], y[i], z[i])[2];

        k2 = dt*f(sigma, rho, beta, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4)[0];
        l2 = dt*f(sigma, rho, beta, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4)[1];
        m2 = dt*f(sigma, rho, beta, t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4)[2];

        k3 = dt*f(sigma, rho, beta, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32)[0];
        l3 = dt*f(sigma, rho, beta, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32)[1];
        m3 = dt*f(sigma, rho, beta, t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32)[2];

        k4 = dt*f(sigma, rho, beta, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[0];
        l4 = dt*f(sigma, rho, beta, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[1];
        m4 = dt*f(sigma, rho, beta, t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197)[2];

        k5 = dt*f(sigma, rho, beta, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[0];
        l5 = dt*f(sigma, rho, beta, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[1];
        m5 = dt*f(sigma, rho, beta, t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104)[2];

        k6 = dt*f(sigma, rho, beta, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[0];
        l6 = dt*f(sigma, rho, beta, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[1];
        m6 = dt*f(sigma, rho, beta, t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40)[2];

        // x1, y1 and z1 are our fourth order approximations
        x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
        y1 = y[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
        z1 = z[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5;

        // x2, y2 and z2 are our fifth order approximations
        x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
        y2 = y[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;
        z2 = z[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55;

        // The following are used to correct the step size
        Rx = Math.abs(x1-x2)/dt;
        Ry = Math.abs(y1-y2)/dt;
        Rz = Math.abs(z1-z2)/dt;
        sx = 0.84*Math.pow(epsilon/Rx, 1/4);                
        sy = 0.84*Math.pow(epsilon/Ry, 1/4);
        sz = 0.84*Math.pow(epsilon/Rz, 1/4);
        R = Math.max(Rx, Ry, Rz);
        s = Math.min(sx, sy, sz);
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

    // Write t, x and y to our solution object
    solution = {
        t: t,
        x: x,
        y: y,
        z: z
    };
}

/**
 * Tabulates solution data.
 *
 * @params           None. Uses the entries of the solution object, however. 
 * @return           Nothing. Just populates the table with the solution values. 
 */
function fillTable() {
    if ( solution["t"].length == 0) {
        alert('You haven\'t solved the problem yet! Press the "Solve the problem" button before pressing the "Tabulate the solution" button again.');
        return
    }
    t = solution["t"];
    x = solution["x"];
    y = solution["y"];
    z = solution["z"];
    document.getElementById('tableOutputs').innerHTML = '';
    tableContents = '<tr>';
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
    document.getElementById('tableOutputs').innerHTML = '';
}

/**
 * Generate two plots:
 * - one of y and x against t; and
 * - a phase plot of y against x.
 * 
 * @params           None.
 * @return           Nothing. Just generates the plots.
 */
function generatePlots() {
    if ( solution["t"].length == 0) {
        alert('You haven\'t solved the problem yet! Before pressing the "Plot the solution" button again, press the "Solve the problem" button.');
        return
    }
    t = solution["t"];
    x = solution["x"];
    y = solution["y"];
    z = solution["z"];

    // Height and width of plots
    windowInnerWidth  = window.innerWidth;
    windowInnerHeight = window.innerHeight;
    document.getElementById("timePlot").style = "height: " + windowInnerHeight + "px;";
    document.getElementById("phasePlotXYZ").style = "height: " + windowInnerHeight + "px;";
    document.getElementById("phasePlotXY").style = "height: " + windowInnerHeight + "px;";
    document.getElementById("phasePlotXZ").style = "height: " + windowInnerHeight + "px;";
    document.getElementById("phasePlotYZ").style = "height: " + windowInnerHeight + "px;";

    // Characteristics of the x and y against time plot
    var plot1 = {
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

    var plot2 = {
        x: t,
        y: x,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'x'
    };

    var plot3 = {
        x: t,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'y'
    };

    var plot4 = {
        x: t,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1,
        name: 'z'
    };

    var plotxy = {
        x: x,
        y: y,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };

    var plotxz = {
        x: x,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };

    var plotyz = {
        x: y,
        y: z,
        type: 'scatter',
        mode: 'lines',
        opacity: 1
    };

    var layout1 = {
        title: 'Phase plot of the solution to the Lorenz equations'
    };

    var layout2 = {
        title: "Time plots of the solution to the problem"
    };

    var layoutxy = {
        title: "xy phase plot"
    };

    var layoutxz = {
        title: "xz phase plot"
    };

    var layoutyz = {
        title: "yz phase plot"
    };

    data1 = [plot1];
    data2 = [plot2, plot3, plot4];
    dataxy = [plotxy];
    dataxz = [plotxz];
    datayz = [plotyz];

    // Generate plots
    Plotly.newPlot('timePlot', data2, layout2);
    Plotly.newPlot('phasePlotXYZ', data1, layout1);
    Plotly.newPlot('phasePlotXY', dataxy, layoutxy);
    Plotly.newPlot('phasePlotXZ', dataxz, layoutxz);
    Plotly.newPlot('phasePlotYZ', datayz, layoutyz);
};

/**
 * Removes solution plots
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution plots.
 */
function removePlots() {
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