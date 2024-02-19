var objectOfInputs = readInputs();
var {dtInitial, t0, x0, y0, z0} = objectOfInputs;
var dt = dtInitial;
var [t, x, y, z] = [[t0], [x0], [y0], [z0]];
var vars = [[x0, y0, z0]];
var i = 0;
var [dt, t, vars, i] = approxRKF45(dt, objectOfInputs, [objectOfInputs.t0], vars, i);
