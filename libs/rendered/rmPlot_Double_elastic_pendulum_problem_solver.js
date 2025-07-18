function removePendulumPlots() {
  removeTimePlot()
  removeR1TPlot()
  removeDr1TPlot()
  removeDr1R1Plot()
  removeR2TPlot()
  removeDr2TPlot()
  removeDr2R2Plot()
  removeR1R2PhasePlot()
  removeTheta1TPlot()
  removeDtheta1TPlot()
  removeDtheta1Theta1Plot()
  removeTheta2TPlot()
  removeDtheta2TPlot()
  removeDtheta2Theta2Plot()
  removeTheta1Theta2PhasePlot()
}
function removeTimePlot() {
  rmPlot("pendulumTimePlot")
}
function removeR1TPlot() {
  rmPlot("timePlot")
}
function removeDr1TPlot() {
  rmPlot("plotR1T")
}
function removeDr1R1Plot() {
  rmPlot("plotDr1T")
}
function removeR2TPlot() {
  rmPlot("plotDr1R1")
}
function removeDr2TPlot() {
  rmPlot("plotR2T")
}
function removeDr2R2Plot() {
  rmPlot("plotDr2T")
}
function removeR1R2PhasePlot() {
  rmPlot("plotDr2R2")
}
function removeTheta1TPlot() {
  rmPlot("plotR2R1")
}
function removeDtheta1TPlot() {
  rmPlot("plotTheta1T")
}
function removeDtheta1Theta1Plot() {
  rmPlot("plotDtheta1T")
}
function removeTheta2TPlot() {
  rmPlot("plotDtheta1Theta1")
}
function removeDtheta2TPlot() {
  rmPlot("plotTheta2T")
}
function removeDtheta2Theta2Plot() {
  rmPlot("plotDtheta2T")
}
function removeTheta1Theta2PhasePlot() {
  rmPlot("plotDtheta2Theta2")
}
function removeAnimation() {
  rmPlot("animation")
}
function removeR1PhaseAnimation() {
  rmPlot("animationR1Phase")
}
function removeR2PhaseAnimation() {
  rmPlot("animationR2Phase")
}
function removeTheta1PhaseAnimation() {
  rmPlot("animationTheta1Phase")
}
function removeTheta2PhaseAnimation() {
  rmPlot("animationTheta2Phase")
}
function removeAnimations() {
  removeAnimation()
  removeR1PhaseAnimation()
  removeR2PhaseAnimation()
  removeTheta1PhaseAnimation()
  removeTheta2PhaseAnimation()
}
