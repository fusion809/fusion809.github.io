function removeDtheta1Theta1PhasePlot() {
  rmPlot("phasePlotDtheta1Theta1")
}
function removeDtheta2Theta2PhasePlot() {
  rmPlot("phasePlotDtheta2Theta2")
}
function removeDtheta3Theta3PhasePlot() {
  rmPlot("phasePlotDtheta3Theta3")
}
function removeTheta1Theta2PhasePlot() {
  rmPlot("phasePlotTheta1Theta2")
}
function removeTheta1Dtheta2PhasePlot() {
  rmPlot("phasePlotDtheta2Theta1")
}
function removeTheta2Dtheta1PhasePlot() {
  rmPlot("phasePlotDtheta1Theta2")
}
function removeDtheta1Dtheta2PhasePlot() {
  rmPlot("phasePlotDtheta2Dtheta1")
}
function removePendulumPlots() {
}
function removeTimePlot() {
  rmPlot("pendulumTimePlot")
}
function removePlots() {
  removeDtheta1Theta1PhasePlot()
  removeDtheta2Theta2PhasePlot()
  removeDtheta3Theta3PhasePlot()
  removeTheta1Theta2PhasePlot()
  removeTheta1Dtheta2PhasePlot()
  removeTheta2Dtheta1PhasePlot()
  removeDtheta1Dtheta2PhasePlot()
  removePendulumPlots()
  removeTimePlot()
}
function removeAnimation() {
  rmPlot("animation")
}
function removeTheta1PhaseAnimation() {
  rmPlot("animationTheta1Phase")
}
function removeTheta2PhaseAnimation() {
  rmPlot("animationTheta2Phase")
}
function removeTheta3PhaseAnimation() {
  rmPlot("animationTheta3Phase")
}
function removeThetasPhaseAnimation() {
  rmPlot("animationThetasPhase")
}
function removeAnimations() {
  removeAnimation()
  removeTheta1PhaseAnimation()
  removeTheta2PhaseAnimation()
  removeTheta3PhaseAnimation()
  removeThetasPhaseAnimation()
}
