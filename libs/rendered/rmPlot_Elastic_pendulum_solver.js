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
  removeZThetaPhasePlot()
  removeDzZPhasePlot()
  removeDthetaThetaPhasePlot()
  removeTimePlot()
  removePendulumPlot()
  removePendulumTimePlot()
  removePlots()
}
function removePlots() {
  removeZThetaPhasePlot()
  removeDzZPhasePlot()
  removeDthetaThetaPhasePlot()
  removeTimePlot()
  removePendulumPlot()
  removePendulumTimePlot()
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
