function removeXPhasePlot() {
  rmPlot("xPhasePlot")
}
function removeYPhasePlot() {
  rmPlot("yPhasePlot")
}
function removeXYPlot() {
  rmPlot("xYPlot")
}
function removeTimePlot() {
  rmPlot("timePlot")
}
function removePlots() {
  removeXPhasePlot()
  removeYPhasePlot()
  removeXYPlot()
  removeTimePlot()
}
function removeXAnimation() {
  rmPlot("animationX")
}
function removeYAnimation() {
  rmPlot("animationY")
}
function removeXYAnimation() {
  rmPlot("animationXY")
}
function removeAnimations() {
  removeXAnimation()
  removeYAnimation()
  removeXYAnimation()
}
