function removeRThetaPhasePlot() {
  rmPlot("phasePlotRTheta")
}
function removeDrRPhasePlot() {
  rmPlot("phasePlotDrR")
}
function removeDthetaThetaPhasePlot() {
  rmPlot("phasePlotDthetaTheta")
}
function removeTimePlot() {
  rmPlot("timePlot")
}
function removePendulumPlot() {
  rmPlot("pendulumPlot")
}
function removePendulumTimePlot() {
  rmPlot("pendulumTimePlot")
}
function removePendulumPlots() {
  removePendulumPlot()
  removePendulumTimePlot()
}
function removePlots() {
  removeRThetaPhasePlot()
  removeDrRPhasePlot()
  removeDthetaThetaPhasePlot()
  removeTimePlot()
  removePendulumPlots()
}
function removeAnimation() {
  rmPlot("animation")
}
function removeThetaPhaseAnimation() {
  rmPlot("animationThetaPhase")
}
function removeRPhaseAnimation() {
  rmPlot("animationRPhase")
}
function removeAnimations() {
  removeAnimation()
  removeThetaPhaseAnimation()
  removeRPhaseAnimation()
}
