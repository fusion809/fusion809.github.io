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
 * Return an object of parameter values
 * 
 * @param elements      Elements as extracted using document.querySelectorAll('form')[0]
 * @return              ObjectOfInputs
 */
var formDataExtractor = elements => [].reduce.call(elements, (objectOfInputs, element) => {
    // If it is a number, but in string form, convert it to a number
    if (!isNaN(element.value)) {
        objectOfInputs[element.name] = parseFloat(element.value);
    } else {
        objectOfInputs[element.name] = element.value;
    }
    return objectOfInputs;
}, {});

/**
 * Read inputs from the form and enter them into an object.
 * 
 * @params              None.
 * @return              An object containing all the form inputs.
 */
function readInputs() {
    // Form elements
    var form = document.querySelectorAll('form')[0];
    var objectOfInputs = formDataExtractor(form);
    return objectOfInputs;
}

/**
 * Add 1d arrays contained with a 2d array
 * 
 * @param arr      A 2d array with each column being a different array to be added.
 * @return         1d array whose elements are the sum of corresponding elements of the 1d subarrays of arr.
 */
function arrAdd(arr) {
    var finalArr = new Array(n1);
    var n1 = arr.length;
    var n2 = arr[0].length;

    for (let i = 0; i < n2 ; i++) {
        finalArr[i] = arr[0][i];
        for (let j = 1; j < n1; j++) {
            finalArr[i] += arr[j][i];
        }
    }

    return finalArr;
}

/**
 * Array to be multiplied by a scalar
 * 
 * @param arr      Array to be multiplied
 * @param scalar   Scalar by which array is to be multiplied.
 */
function arrMult(arr, scalar) {
    var n = arr.length;
    var finalArr = new Array(n);

    for (let j = 0; j < n; j++) {
        finalArr[j] = arr[j]*scalar;
    }

    return finalArr;
}

/**
 * Divide arr by scalar
 * 
 * @param arr      The array whose elements are to be divided.
 * @param scalar   The scalar by which arr's elements are to be divided.
 * @return         Divided array.
 */
function arrDiv(arr, scalar) {
    return arrMult(arr, 1/scalar);
}

function approxRKF45(dt, objectOfInputs, t, vars, i) {
    // Initialize variables
    var K = [[]];
    var X1 = [];
    var X2 = [];
    var Rarr = [];
    var sarr = [];
    var epsilon = objectOfInputs.epsilon;

    // K[i][j], i goes from 0 to 5 and represents k1, k2, k3, k4, k5, k6
    // j goes from 0 to the dimensionality of the problem - 1
    K[0] = arrMult(f(objectOfInputs, t[i], vars[i]), dt);
    K[1] = arrMult(f(objectOfInputs, t[i] + dt/4, arrAdd([vars[i], arrDiv(K[0], 4)])), dt);
    K[2] = arrMult(f(objectOfInputs, t[i] + 3*dt/8, arrAdd([vars[i], arrMult(K[0], 3/32), arrMult(K[1], 9/32)])), dt);
    K[3] = arrMult(f(objectOfInputs, t[i] + 12*dt/13, arrAdd([vars[i], arrMult(K[0], 1932/2197), arrMult(K[1], -7200/2197), arrMult(K[2], 7296/2197)])), dt);
    K[4] = arrMult(f(objectOfInputs, t[i]+dt, arrAdd([vars[i], arrMult(K[0], 439/216), arrMult(K[1], -8), arrMult(K[2], 3680/513), arrMult(K[3], -845/4104)])), dt);
    K[5] = arrMult(f(objectOfInputs, t[i] + dt/2, arrAdd([vars[i], arrMult(K[0], -8/27), arrMult(K[1], 2), arrMult(K[2], -3544/2565), arrMult(K[3], 1859/4104), arrMult(K[4], -11/40)])), dt);

    // X1[j] are our fourth order approxs & X2[j] are our 5th order approxs
    for (let j = 0; j < K[0].length; j++) {
        X1[j] = vars[i][j] + 25*K[0][j]/216+1408*K[2][j]/2565+2197*K[3][j]/4104-K[4][j]/5;
        X2[j] = vars[i][j] + 16*K[0][j]/135+6656*K[2][j]/12825+28561*K[3][j]/56430-9*K[4][j]/50+2*K[5][j]/55;
        Rarr[j] = Math.abs(X1[j]-X2[j])/dt;
        sarr[j] = 0.84*Math.pow(epsilon/Rarr[j], 1/4);  
    }

    // R is our error estimate and s is how much our step size needs to
    // be adjusted by.
    var R = Math.max(...Rarr);
    var s = Math.min(...sarr);
    // If R is less than or equal to epsilon move onto the next step
    if ( R <= epsilon ) {
        t.push(t[i]+dt);
        vars.push(X1);
        i++;
        dt *= s;
    } else {
        dt *= s;
    }

    return [dt, t, vars, i];
}

/**
 * Solves the problem using RKF45.
 * 
 * @param objectOfInputs An object containing all the problem parameters.
 * @param vars0          Initial conditions of the problem in a 2D array (e.g. [[x0, y0]]).
 * @return               [t, vars] where vars is an array of solution values.
 */
function RKF45Body(objectOfInputs, vars0) {
    // Initialize vars
    var {t0, tf, dtInitial} = objectOfInputs;
    var vars = vars0;
    var t = [t0];
    var dt = dtInitial;
    var i = 0;

    // Loop over each step until we reach the endpoint
    while ( t[i] < tf ) {
        dt = Math.min(dt, tf-t[i]);
        [dt, t, vars, i] = approxRKF45(dt, objectOfInputs, t, vars, i);
    }

    // Transpose vars
    vars = vars[0].map((_, colIndex) => vars.map(row => row[colIndex]));

    // Return t and vars
    return [t, vars];
}

/** 
 * Solve the problem using RKF45.
 *
 * @param func           RKF45 function to be used.
 * @param objectOfInputs An object containing all the problem parameters.
 * @return               A solution object.
 */
function solveProblem(func, objectOfInputs) {
    // Solve the problem
    var [t, vars] = func(objectOfInputs);
    // Write t, x and xDot to our solution object
    var solution = {
        t: t,
        vars: vars
    };
    return solution;
}

/**
 * Tabulates solution data.
 *
 * @param objectOfInputs An object containing all the form parameters. 
 * @param headings       Headings for each dependent variable column of the table.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTable(objectOfInputs, headings) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    var {t, vars} = solution;
    if (vars.length != headings.length) {
        alert("libs/common.js#fillTable: vars length and headings length do not match!")
    }
    var epsilon = objectOfInputs.epsilon;

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    for (let j = 0; j < headings.length ; j++) {
        tableContents += '<th>' + headings[j] + '</th>';
    }
    tableContents += "</tr>";
    for (let j = 0; j < vars[0].length; j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        for (let k = 0; k < vars.length ; k++) {
            tableContents += '<td>' + vars[k][j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        }
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
};