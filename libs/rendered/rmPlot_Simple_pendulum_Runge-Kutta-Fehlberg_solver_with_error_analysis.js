function removePhasePlot() {
  rmPlot("phasePlot")
}
function removeTimePlot() {
  rmPlot("timePlot")
}
function removeErrorPlot() {
  rmPlot("errorPlot")
}
function removePlots() {
  removePhasePlot()
  removeTimePlot()
  removeErrorPlot()
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
