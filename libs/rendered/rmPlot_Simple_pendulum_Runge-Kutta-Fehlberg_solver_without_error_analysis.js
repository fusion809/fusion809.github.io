function removePhasePlot() {
  rmPlot("phasePlot")
}
function removeTimePlot() {
  rmPlot("timePlot")
}
function removePlots() {
  removePhasePlot()
  removeTimePlot()
}
function removeAnimation() {
  rmPlot("animation")
}
function removePhaseAnimation() {
  rmPlot("animationPhase")
}
function removeAnimations() {
  removeAnimation()
  removePhaseAnimation()
}
