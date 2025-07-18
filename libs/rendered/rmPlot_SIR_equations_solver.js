function removeSIRPhasePlot() {
  rmPlot("phasePlotSIR")
}
function removeSIPhasePlot() {
  rmPlot("phasePlotSI")
}
function removeSRPhasePlot() {
  rmPlot("phasePlotSR")
}
function removeIRPhasePlot() {
  rmPlot("phasePlotIR")
}
function removeTimePlot() {
  rmPlot("timePlot")
}
function removePlots() {
  removeSIRPhasePlot()
  removeSIPhasePlot()
  removeSRPhasePlot()
  removeIRPhasePlot()
  removeTimePlot()
}
function removeAnimation() {
  rmPlot("animation")
}
