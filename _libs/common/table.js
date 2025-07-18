/**
 * Tabulates solution data.
 *
 * @param objectOfInputs An object containing all the form parameters. 
 * @param headings       Headings for each dependent variable column of the table.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTable(objectOfInputs, headings, solution=undefined) {
    // Solve the problem
    if (solution == undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    if (solution.varsLen() != headings.length) {
        alert("libs/common.js#fillTable: vars length and headings length do not match!")
    }
    var epsilon = objectOfInputs.epsilon;

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<table><tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    for (let j = 0; j < headings.length ; j++) {
        tableContents += '<th>' + headings[j] + '</th>';
    }
    tableContents += "</tr>";
    for (let j = 0; j < solution.tLen(); j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + solution.t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        for (let k = 0; k < solution.varsLen() ; k++) {
            tableContents += '<td>' + solution.extract(k)[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        }
        tableContents += '</tr>';
    }
    tableContents += "</table>"
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

function generateAllOutputs(objectOfInputs=undefined, solution=undefined) {
    if (objectOfInputs==undefined) {
        var objectOfInputs = readInputs();
    }
    if (solution==undefined) {
        var solution = solveProblem(RKF45, objectOfInputs);
    }
    generateTable(objectOfInputs, solution);
    generatePlots(objectOfInputs, solution);
    try {
        generateAnimations(objectOfInputs, solution);
    } catch(e) {
        generateAnimation(objectOfInputs, solution);
    }
}

function removeAllOutputs() {
    removeTable();
    removePlots();
    try {
        removeAnimations()
    } catch(e) {
        removeAnimation();
    }
}