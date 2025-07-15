@def title="Double elastic pendulum problem solver"
@def hassim=true;
@def params = (g=(val=9.81, desc="Acceleration due to gravity in \\(\\mathrm{m}\\cdot \\mathrm{s}^{-2}\\)."),l1=(val=1, desc="Rest length of pendulum 1 in metres."),l2=(val=1, desc="Rest length of pendulum 2 in metres."),m1=(val=1, desc="Mass of pendulum bob 1 in kilograms."),m2=(val=1, desc="Mass of pendulum bob 2 in kilograms."),k1=(val=10, desc="Coefficient for first pendulum spring."),k2=(val=10, desc="Coefficient for second pendulum spring."),b1=(val=0, desc="Linear dissipation coefficient for pendulum 1."),b2=(val=0, desc="Linear dissipation coefficient for pendulum 2."),c1=(val=0, desc="Quadratic dissipation coefficient for pendulum 1."),c2=(val=0, desc="Quadratic dissipation coefficient for pendulum 2."),t0=(val=0, desc="Starting time for the simulation in seconds (s)."),tf=(val=120, desc="End time for the simulation in seconds."),r10=(val=1, desc="Initial value of \\(r_1\\) in metres."),dr10=(val=1, desc="Initial value of \\(\\dot{r}_1\\) in \\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\)."),r20=(val=1, desc="Initial value of \\(r_2\\) in metres."),dr20=(val=1, desc="Initial value of \\(\\dot{r}_2\\) in \\(\\mathrm{m}\\cdot \\mathrm{s}^{-1}\\)."),theta10=(val=1, desc="Initial value of \\(\\theta_1\\) in radians."),dtheta10=(val=1, desc="Initial value of \\(\\dot{\\theta}_1\\) in radians per second."),theta20=(val=1, desc="Initial value of \\(\\theta_2\\) in radians."),dtheta20=(val=1, desc="Initial value of \\(\\dot{\\theta}_2\\) in radians per second."),epsilon=(val=1e-6,))

@def labels=["Tabulate the solution","Remove the solution table","Generate pendulum position plots","Remove pendulum position plots","Generate a time plot of all variables (+time derivatives)","Remove time plot of all variables (+time derivatives)","Generate a \\(r_1\\) against time plot","Remove \\(r_1\\) against time plot","Generate a \\(\\dot{r}_1\\) against time plot","Remove \\(\\dot{r}_1\\) against time plot","Generate a phase plot of \\(\\dot{r}_1\\) vs \\(r_1\\)","Remove phase plot of \\(\\dot{r}_1\\) vs \\(r_1\\)","Generate a \\(r_2\\) against time plot","Remove \\(r_2\\) against time plot","Generate a \\(\\dot{r}_2\\) against time plot","Remove \\(\\dot{r}_2\\) against time plot","Generate a phase plot of \\(\\dot{r}_2\\) vs \\(r_2\\)","Remove phase plot of \\(\\dot{r}_2\\) vs \\(r_2\\)","Generate a \\(r_2\\) against \\(r_1\\) phase plot","Remove \\(r_2\\) against \\(r_1\\) phase plot","Generate a \\(\\theta_1\\) against time plot","Remove \\(\\theta_1\\) against time plot","Generate a \\(\\dot{\\theta}_1\\) against time plot","Remove \\(\\dot{\\theta}_1\\) against time plot","Generate a phase plot of \\(\\dot{\\theta}_1\\) vs \\(\\theta_1\\)","Remove phase plot of \\(\\dot{\\theta}_1\\) vs \\(\\theta_1\\)","Generate a \\(\\theta_2\\) against time plot","Remove \\(\\theta_2\\) against time plot","Generate a \\(\\dot{\\theta}_2\\) against time plot","Remove \\(\\dot{\\theta}_2\\) against time plot","Generate a phase plot of \\(\\dot{\\theta}_2\\) vs \\(\\theta_2\\)","Remove phase plot of \\(\\dot{\\theta}_2\\) vs \\(\\theta_2\\)","Generate a \\(\\theta_2\\) against \\(\\theta_1\\) phase plot","Remove \\(\\theta_2\\) against \\(\\theta_1\\) phase plot","Generate all solution plots","Remove all plots","Generate an animation of the system","Remove system animation","Generate a \\(r_1\\) phase plot animation","Remove \\(r_1\\) phase plot animation","Generate a \\(r_2\\) phase plot animation","Remove \\(r_2\\) phase plot animation","Generate a \\(\\theta_1\\) phase plot animation","Remove \\(\\theta_1\\) phase plot animation","Generate a \\(\\theta_2\\) phase plot animation","Remove \\(\\theta_2\\) phase plot animation","Generate all animations","Remove all animations"];

@def ids=["tableOutputs", "pendulumPlot", "pendulumTimePlot", "timePlot", "plotR1T", "plotDr1T", "plotDr1R1", "plotR2T", "plotDr2T", "plotDr2R2", "plotR2R1", "plotTheta1T", "plotDtheta1T", "plotDtheta1Theta1", "plotTheta2T", "plotDtheta2T", "plotDtheta2Theta2", "plotTheta1Theta2", "animation", "animationR1Phase", "animationR2Phase", "animationTheta1Phase", "animationTheta2Phase", "animation"]
@def funcs=["generateTable()","removeTable()","generatePendulumPlots(readInputs(), solveProblem(RKF45, readInputs()))","removePendulumPlots()","generateTimePlot(solveProblem(RKF45, readInputs()))","removeTimePlot()","generateR1TPlot(solveProblem(RKF45, readInputs()))","removeR1TPlot()","generateDr1TPlot(solveProblem(RKF45, readInputs()))","removeDr1TPlot()","generateDr1R1Plot(solveProblem(RKF45, readInputs()))","removeDr1R1Plot()","generateR2TPlot(solveProblem(RKF45, readInputs()))","removeR2TPlot()","generateDr2TPlot(solveProblem(RKF45, readInputs()))","removeDr2TPlot()","generateDr2R2Plot(solveProblem(RKF45, readInputs()))","removeDr2R2Plot()","generateR1R2PhasePlot(solveProblem(RKF45, readInputs()))","removeR1R2PhasePlot()","generateTheta1TPlot(solveProblem(RKF45, readInputs()))","removeTheta1TPlot()","generateDtheta1TPlot(solveProblem(RKF45, readInputs()))","removeDtheta1TPlot()","generateDtheta1Theta1Plot(solveProblem(RKF45, readInputs()))","removeDtheta1Theta1Plot()","generateTheta2TPlot(solveProblem(RKF45, readInputs()))","removeTheta2TPlot()","generateDtheta2TPlot(solveProblem(RKF45, readInputs()))","removeDtheta2TPlot()","generateDtheta2Theta2Plot(solveProblem(RKF45, readInputs()))","removeDtheta2Theta2Plot()","generateTheta1Theta2PhasePlot(solveProblem(RKF45, readInputs()))","removeTheta1Theta2PhasePlot()","generatePlots(readInputs())","rmPlots()","generateAnimation()","removeAnimation()","generateR1PhaseAnimation()","removeR1PhaseAnimation()","generateR2PhaseAnimation()","removeR2PhaseAnimation()","generateTheta1PhaseAnimation()","removeTheta1PhaseAnimation()","generateTheta2PhaseAnimation()","removeTheta2PhaseAnimation()","generateAnimations()","removeAnimations()"]

This webpage uses the [Runge-Kutta-Fehlberg fourth-order method with fifth-order error checking (RKF45)](/RKF45/) to approximate the solution to the problem of the double elastic pendulum. 

The equations that will be integrated here are derived in [this article](/doubleElasticPendulum/derivation/).
~~~
<figure>
<img src="/doubleElasticPendulum/Double elastic pendulum.svg" width="500px"></img>
<figcaption style="font-weight: bold; font-size: 18px;">Figure 1: Diagram of the double elastic pendulum.</b></figcaption>
</figure>
~~~

            The ordinary differential equation system being solved is

\begin{aligned}
        \begin{bmatrix}
\ddot{r}_1 \\
\ddot{r}_2 \\
\ddot{\theta}_1 \\
\ddot{\theta}_2
        \end{bmatrix} &= \begin{bmatrix}
        1 & \dfrac{m_2\cos{\Delta}}{m_1+m_2} & 0 & -\dfrac{m_2r_2\sin{\Delta}}{m_1+m_2} \\
        \cos{\Delta} & 1 & r_1\sin{\Delta} & 0 \\
        0 & \dfrac{m_2\sin{\Delta}}{(m_1+m_2)r_1} & 1 & \dfrac{m_2r_2\cos{\Delta}}{(m_1+m_2)r_1} \\
        -\dfrac{\sin{\Delta}}{r_2} & 0 & \dfrac{r_1\cos{\Delta}}{r_2} & 1
        \end{bmatrix}^{-1} \begin{bmatrix}
        r_1\dot{\theta}_1^2-g\sin{\theta_1} + \dfrac{m_2}{m_1+m_2}\left(r_2\dot{\theta}_2^2\cos{\Delta} + 2\dot{r}_2\dot{\theta}_2\sin{\Delta}\right)  + \dfrac{Q_{r_1}-k_1(r_1-l_1)}{m_1+m_2}\\
        r_2\dot{\theta}_2^2-g\sin{\theta_2} + r_1\dot{\theta}_1^2\cos{\Delta} - 2\dot{r}_1\dot{\theta}_1\sin{\Delta} + \dfrac{Q_{r_2}-k_2(r_2-l_2)}{m_2} \\
        -\dfrac{2\dot{r}_1\dot{\theta}_1}{r_1} - \dfrac{g\cos{\theta_1}}{r_1} - \dfrac{m_2}{(m_1+m_2)r_1}\left[2\dot{r}_2\dot{\theta}_2\cos{\Delta} -r_2\dot{\theta}_2^2\sin{\Delta}\right] + \dfrac{Q_{\theta_1}}{(m_1+m_2)r_1^2} \\
        -\dfrac{2\dot{r}_2\dot{\theta}_2}{r_2}- \dfrac{g\cos{\theta_2}}{r_2} - \dfrac{2\dot{r}_1\dot{\theta}_1\cos{\Delta}}{r_2} - \dfrac{r_1\dot{\theta}_1^2\sin{\Delta}}{r_2} + \dfrac{Q_{\theta_2}}{m_2r_2^2}
        \end{bmatrix}.
\end{aligned}
~~~
{{ insert template.html}}
~~~