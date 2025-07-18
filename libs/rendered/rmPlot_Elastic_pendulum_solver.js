function removeZThetaPhasePlot() {
  rmPlot("phasePlotZTheta")
}
function removeDzZPhasePlot() {
  rmPlot("phasePlotDzZ")
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
  removeZThetaPhasePlot()
  removeDzZPhasePlot()
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
function removeZPhaseAnimation() {
  rmPlot("animationZPhase")
}
function removeAnimations() {
  removeAnimation()
  removeThetaPhaseAnimation()
  removeZPhaseAnimation()
}
