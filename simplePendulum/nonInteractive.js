// RHS of our ODE system
function f(t, theta, dtheta) {
    return [dtheta, -9.8*Math.cos(theta)];
}

// Our initial conditions and time bounds
var t0 = 0.0;
// The period of our problem, ellipk(1/2)/sqrt(2.45)
var tf = 2.3690497221753453;
var theta0 = 0;
var dtheta0 = 0;

// Details of our numerical approximation
// epsilon is our error tolerance
var epsilon = 1e-12;
var dtInitial = 0.1;

// Initialize arrays, dt and i
var t = [t0];
var theta = [theta0];
var dtheta = [dtheta0];
var dt = dtInitial;
var i = 0;

while ( t[i] < tf ) {
    dt = Math.min(dt, tf-t[i]);

    // Runge-Kutta-Fehlberg approximators
    k1 = dt*f(t[i], theta[i], dtheta[i])[0];
    l1 = dt*f(t[i], theta[i], dtheta[i])[1];
    k2 = dt*f(t[i]+dt/4, theta[i]+k1/4, dtheta[i]+l1/4)[0];
    l2 = dt*f(t[i]+dt/4, theta[i]+k1/4, dtheta[i]+l1/4)[1];
    k3 = dt*f(t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, dtheta[i]+3*l1/32+9*l2/32)[0];
    l3 = dt*f(t[i]+3*dt/8, theta[i]+3*k1/32+9*k2/32, dtheta[i]+3*l1/32+9*l2/32)[1];
    k4 = dt*f(t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, dtheta[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[0];
    l4 = dt*f(t[i]+12*dt/13, theta[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, dtheta[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197)[1];
    k5 = dt*f(t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, dtheta[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[0];
    l5 = dt*f(t[i]+dt, theta[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, dtheta[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104)[1];
    k6 = dt*f(t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, dtheta[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[0];
    l6 = dt*f(t[i]+dt/2, theta[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, dtheta[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40)[1];

    // theta1 and dtheta1 are our fourth order approximations
    theta1 = theta[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5;
    dtheta1 = dtheta[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5;
    // theta2 and dtheta2 are our fifth order approximations
    theta2 = theta[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55;
    dtheta2 = dtheta[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55;

    // The following are used to correct the step size
    Rtheta = Math.abs(theta1-theta2)/dt;
    Rdtheta = Math.abs(dtheta1-dtheta2)/dt;
    stheta = 0.84*Math.pow(epsilon/Rtheta, 1/4);                
    sdtheta = 0.84*Math.pow(epsilon/Rdtheta, 1/4);
    R = Math.max(Rtheta, Rdtheta);
    s = Math.min(stheta, sdtheta);
    if ( R <= epsilon ) {
        t.push(t[i]+dt);
        theta.push(theta1);
        dtheta.push(dtheta1);
        i++;
        dt *= s;
    } else {
        dt *= s;
    }
}

// Put our data into a table
document.write("<table border='1px'>")
document.write("<tr>")
document.write("<th>Index</th>")
document.write("<th>t</th>")
document.write("<th>theta</th>")
document.write("<th>dtheta</th>")
document.write("</tr>")
for (let j = 0; j < theta.length; j++) {
    document.write("<tr>")
    document.write("<td>" + j + "</td>")
    document.write("<td>" + t[j] + "</td>");
    document.write("<td>" + theta[j] + "</td>");
    document.write("<td>" + dtheta[j] + "</td>");
    document.write("</tr>");
}
document.write("</table>");

// Create our plot using Plotly
var plot1 = {
    x: t,
    y: theta,
    type: 'scatter',
    name: "Theta"
};

var plot2 = {
    x: t,
    y: dtheta,
    type: 'scatter',
    name: "dtheta"
};

var layout1 = {
    title: 'dtheta and theta against time plots',
    xaxis: {
        title: 'Time'
    }
}

data1 = [plot1, plot2];

var plot3 = {
    x: theta,
    y: dtheta,
    type: 'scatter',
    name: "Phase plot"
};

var layout2 = {
    title: "Phase plot of dtheta against theta",
    xaxis: {
        title: "theta"
    },
    yaxis: {
        title: "dtheta"
    }
};

data2 = [plot3];

document.write("<div id='myDiv' style='width:1000px; height:700px;'></div>");
Plotly.newPlot('myDiv', data1, layout1);
document.write("<div id='phasePlot' style='width:1000px; height:700px;'></div>");
Plotly.newPlot('phasePlot', data2, layout2);