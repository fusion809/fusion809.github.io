/**
 * Tabulates solution data.
 *
 * @param objectOfInputs An object containing all the form parameters. 
 * @param headings       Headings for each dependent variable column of the table.
 * @return               Nothing. Just populates the table with the solution values. 
 */
function fillTable(objectOfInputs, headings) {
    // Solve the problem
    var solution = solveProblem(RKF45, objectOfInputs);
    if (solution.varsLen() != headings.length) {
        alert("libs/common.js#fillTable: vars length and headings length do not match!")
    }
    var epsilon = objectOfInputs.epsilon;

    // Write to table
    document.getElementById('tableOutputs').innerHTML = '';
    var tableContents = '<table><tr>';
    tableContents += '<th>Index</th>';
    tableContents += '<th>t (seconds)</th>';
    for (let j = 0; j < headings.length ; j++) {
        tableContents += '<th>' + headings[j] + '</th>';
    }
    tableContents += "</tr>";
    for (let j = 0; j < solution.tLen(); j++) {
        tableContents += '<tr>';
        tableContents += '<td>' + j + '</td>';
        tableContents += '<td>' + solution.t[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        for (let k = 0; k < solution.varsLen() ; k++) {
            tableContents += '<td>' + solution.extract(k)[j].toFixed(Math.ceil(Math.log10(1/epsilon))) + '</td>';
        }
        tableContents += '</tr>';
    }
    tableContents += "</table>"
    document.getElementById('tableOutputs').innerHTML = tableContents;
}

/**
 * Removes the solution table
 * 
 * @params           None.
 * @return           Nothing. Just removes the solution table.
 */
function removeTable() {
    // Clear table content
    document.getElementById('tableOutputs').innerHTML = '';
}

function toggleTOC() {
  const toc = document.getElementById("toc-panel");
  if (!toc) {
    // Panel not yet loaded; retry shortly
    setTimeout(toggleTOC, 10);
    return;
  }
  toc.style.display = (toc.style.display === "none" || toc.style.display === "") ? "block" : "none";
}

function addH1ToTOC(tocPanel, headers) {
    headers.forEach(h => {
        if (!h.id) h.id = h.textContent.trim().replace(/\s+/g, "-").toLowerCase();
        let label = h.textContent.trim();
        if (!label) return;
        if (h.tagName === "H1") {
            var text = '<span id="h1_toc">' + label + '</span>'
        }
        tocPanel.insertAdjacentHTML("afterbegin", `<a href="#${h.id}">${text}</a>`);
    });
}
function addHeadersToTOC(tocList, headers) {
    headers.forEach(h => {
        if (!h.id) h.id = h.textContent.trim().replace(/\s+/g, "-").toLowerCase();
        const li = document.createElement("li");
        let label = h.textContent.trim();
        if (!label) return;
        if (h.tagName === "H1") {
            return;
        } else {
            var id = (h.tagName).toLowerCase() + "_toc";
            var text = `<span id="${id}">` + label + '</span>';
        }
        li.innerHTML = `<a href="#${h.id}">${text}</a>`;
        tocList.appendChild(li);
    });
}

function addTablesToTOC(tocList, tables) {
    tables.forEach(table => {
        const id = table.id;
        const suffix = id;
        const infoEl = document.getElementById("info" + suffix);
        var label;
        if (id == "parameterForm") {
            label = "Simulation parameter form"
        } else if (id == "buttontable") {
            label = "Table of plot, tabulation and animation trigger buttons."
        } else if (infoEl) {
            label = infoEl.textContent.trim();
        } else {
            var suff = id.replace(/([A-Z])/g, (match) => ' ' + match.toLowerCase()).trim();
            suff = suff.replace(/ theta/g, ", θ");
            suff = suff.replace(/minimum and maximum/, "minimum and maximum, ");
            suff = suff + ".";
            label = suff.charAt(0).toUpperCase() + suff.slice(1)
        }
        if (!table.id) return;
        const anchorId = "anchor-" + suffix;
        if (!document.getElementById(anchorId)) {
            const anchor = document.createElement("a");
            anchor.id = anchorId;
            table.parentNode.insertBefore(anchor, table);
            
        }
        const li = document.createElement("li");
        li.innerHTML = `<a href="#${anchorId}">${label.trim()}</a>`;
        tocList.appendChild(li);
    });
}

function generateContLabel(id, content) {
    // If infoEl is missing or empty, fallback to a cleaned-up version of the ID
    const rawSuffix = id.replace(/^container/, "");
    if (rawSuffix == "tableOutputs") {
        label = "Table of outputs"
    } else if (id == "container") {
        label = "System animation"
    } else if (/animationSIR/.test(content)) {
        label = "Animation of S, I and R phase plot."
    } else if (/animationSER/.test(content)) {
        label = "Animation of S, E and R phase plot."
    } else if (/animationEIR/.test(content)) {
        label = "Animation of E, I and R phase plot."
    } else if (/animationSEI/.test(content)) {
        label = "Animation of S, E and I phase plot."
    } else if (/animation/.test(content)) {
        let suff = rawSuffix.replace(/([A-Z])/g, (match) => ' ' + match.toLowerCase()).trim();
        suff = suff.replace(/([a-zA-Z])(\d)$/g, '$1<sub>$2</sub>');
        suff = suff.replace(/(\d)/, "<sub>$1</sub>");
        suff = suff.replace(/theta/g, "θ");
        label = `Animation of ${suff} plot`
    } else if (/errorPlot/.test(content)) {
        label = "Error plot"
    } else if (/phasePlotXYZ/.test(content)) {
        label = "3D phase plot of x, y and z."
    } else if (/phasePlotXY/.test(content)) {
        label = "Phase plot of y vs x."
    } else if (/phasePlotXZ/.test(content)) {
        label = "Phase plot of z vs x."
    } else if (/phasePlotYZ/.test(content)) {
        label = "Phase plot of z vs y."
    } else if (/phasePlotSIR/.test(content)) {
        label = "3D phase plot of S, I and R."
    } else if (/phasePlotEIR/.test(content)) {
        label = "3D phase plot of E, I and R."
    } else if (/phasePlotSER/.test(content)) {
        label = "3D phase plot of S, E and R."
    } else if (/phasePlotSEI/.test(content)) {
        label = "3D phase plot of S, E and I."
    } else if (/phasePlotSI/.test(content)) {
        label = "Phase plot of I vs S."
    } else if (/phasePlotSR/.test(content)) {
        label = "Phase plot of R vs S."
    } else if (/phasePlotIR/.test(content)) {
        label = "Phase plot of R vs I."
    } else if (/phasePlot/.test(content)) {
        let suff = rawSuffix.replace(/phasePlot/, "").replace(/([A-Z])/g, (match) => ' ' + match.toLowerCase()).trim();
        if (suff) {    
            suff = suff.replace(/\s([a-zA-Z]*\d)$/, ' vs $1').replace(/\s([a-zA-Z]*)$/, ' vs $1');
            suff = suff.replace(/([a-zA-Z])(\d)$/g, '$1<sub>$2</sub>');
            suff = suff.replace(/d([a-zA-Z0-9]*)/, "d$1/dt")
            suff = suff.replace(/(\d)/, "<sub>$1</sub>");
            suff = suff.replace(/theta/g, "θ");
            label = `Phase plot of ${suff}`
        } else {
            label = "Phase plot"
        }
    } else if (/pendulum.*Plot/.test(content) || /timePlot/.test(content)) {
        suff = rawSuffix.replace(/([A-Z])/g, (match) => ' ' + match.toLowerCase()).trim();
        label = suff.charAt(0).toUpperCase() + suff.slice(1)
    } else if (/[pP]lot/.test(content)) {
        let suff = rawSuffix.replace(/([A-Z])/g, (match) => ' ' + match.toLowerCase()).trim();
        suff = suff.charAt(0).toUpperCase() + suff.slice(1)
        suff = suff.replace(/[pP]lot /, "Plot of ");
        // Replace digits after a letter with <sub>digit>
        suff = suff.replace(/\s([a-zA-Z]*\d)$/, ' vs $1').replace(/\s([a-zA-Z]*)$/, ' vs $1');
        suff = suff.replace(/d([a-zA-Z0-9]*)/, "d$1/dt")
        suff = suff.replace(/([a-zA-Z])(\d)$/g, '$1<sub>$2</sub>');
        suff = suff.replace(/theta/g, "θ");
        // Add "vs" before the last variable with subscript (or number)
        label = suff.replace(/(\d)/, "<sub>$1</sub>")
    } else if (rawSuffix) {
        label = `Figure ${rawSuffix}`
    } else {
        label = id;
    }
    return label;
}
function addContainersToTOC(tocList, containers) {
    containers.forEach(div => {
        const id = div.id;
        var suffix = id.replace(/^container/, "");
        const infoEl = document.getElementById("info" + suffix);
        let label = infoEl?.textContent.trim();
        const content = document.getElementById(id).innerHTML;
        if (!label) label = generateContLabel(id, content);
        if (!div.id) return;
        const anchorId = "anchor-" + suffix;
        if (!document.getElementById(anchorId)) {
            const anchor = document.createElement("a");
            anchor.id = anchorId;
            div.parentNode.insertBefore(anchor, div);
        }
        const li = document.createElement("li");
        li.innerHTML = `<a href="#${anchorId}">${label.replace(/container/gi, "animation").trim()}</a>`;
        tocList.appendChild(li);
    });
}
function createTOC() {
    const tocPanel = document.getElementById("toc-panel");
    tocPanel.innerHTML = '<ul id="toc-list" style="list-style: none; padding-left: 0;"></ul>';
    const headers = [...document.querySelectorAll('h1, h2, h3:not(div.page-foot h3):not(div[id^="container"] h3), h4, h5')];
    const containers = [...document.querySelectorAll("div[id^='container']")];
    const tables = [...document.querySelectorAll("table")];
    addH1ToTOC(tocPanel, headers);
    const tocList = document.getElementById("toc-list");
    addHeadersToTOC(tocList, headers);
    addTablesToTOC(tocList, tables);
    addContainersToTOC(tocList, containers)
}

document.addEventListener("DOMContentLoaded", function () {
  createTOC();

  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      const id = entry.target.getAttribute("id");
      const tocLinks = document.querySelectorAll("#toc-list a");

      if (entry.isIntersecting && id) {
        tocLinks.forEach(link => link.classList.remove("active"));
        const active = [...tocLinks].find(link => link.getAttribute("href") === `#${id}`);
        if (active) active.classList.add("active");
      }
    });
  }, {
    threshold: 0.1
  });

  // 👇 Observe heading and table anchors directly
  document.querySelectorAll("h1, h2, h3, h4, table").forEach(el => {
    if (!el.id) {
      el.id = el.textContent.trim().replace(/\s+/g, "-").toLowerCase();
    }
    observer.observe(el);
  });

  // 👇 Observe all anchors inserted for containers
  document.querySelectorAll("a[id^='anchor-']").forEach(anchor => {
    observer.observe(anchor);
  });
  const visibleSections = entries
  .filter(entry => entry.isIntersecting)
  .sort((a, b) => a.boundingClientRect.top - b.boundingClientRect.top);

if (visibleSections.length > 0) {
  const id = visibleSections[0].target.getAttribute("id");
  // highlight matching link...
}
});