/* Number and enhance code blocks marked with \codelisting{label}{caption}. */
(function () {
  function nextPre(marker) {
    var walker = document.createTreeWalker(document.body, NodeFilter.SHOW_ELEMENT);
    walker.currentNode = marker;
    while (walker.nextNode()) {
      if (walker.currentNode.tagName === "PRE") return walker.currentNode;
    }
    return null;
  }

  function enhance() {
    var markers = Array.prototype.slice.call(document.querySelectorAll(".code-listing-marker"));
    var listings = [];
    markers.forEach(function (marker) {
      var pre = nextPre(marker);
      if (!pre) return;
      var captionNodes = Array.prototype.map.call(marker.childNodes, function (node) {
        return node.cloneNode(true);
      });
      var label = marker.getAttribute("data-label") || "";
      var listingElement = document.createElement("div");
      listingElement.className = "code-listing";
      if (label) listingElement.id = label;
      var button = document.createElement("button");
      button.type = "button";
      button.className = "code-listing-toggle";
      button.textContent = "Show code";
      button.setAttribute("aria-expanded", "false");
      pre.hidden = true;
      button.addEventListener("click", function () {
        pre.hidden = !pre.hidden;
        button.textContent = pre.hidden ? "Show code" : "Hide code";
        button.setAttribute("aria-expanded", String(!pre.hidden));
      });
      listingElement.addEventListener("code-listing:reveal", function () {
        if (pre.hidden) button.click();
      });
      var parent = pre.parentNode;
      parent.insertBefore(listingElement, pre);
      listingElement.appendChild(button);
      listingElement.appendChild(pre);
      var listingCaption = document.createElement("div");
      listingCaption.className = "code-listing-caption";
      listings.push({ element: listingElement, caption: listingCaption, nodes: captionNodes, label: label });
      listingElement.appendChild(listingCaption);
      marker.remove();
    });

    listings.forEach(function (listing, index) {
      var number = index + 1;
      listing.caption.appendChild(document.createTextNode("Code Listing " + number + ": "));
      listing.nodes.forEach(function (node) { listing.caption.appendChild(node); });
      if (typeof renderMathInElement === "function") {
        renderMathInElement(listing.caption, {
          throwOnError: false,
          delimiters: [
            { left: "$$", right: "$$", display: true },
            { left: "$", right: "$", display: false },
            { left: "\\(", right: "\\)", display: false },
            { left: "\\[", right: "\\]", display: true }
          ]
        });
      }
      if (listing.label) listing.element.id = listing.label;
      document.querySelectorAll(".code-listing-ref").forEach(function (ref) {
        if (ref.getAttribute("data-label") === listing.label) {
          var numberLink = ref.querySelector(".code-listing-number");
          if (!numberLink) return;
          numberLink.href = "#" + listing.label;
          numberLink.textContent = number;
          numberLink.addEventListener("click", function () {
            listing.element.dispatchEvent(new Event("code-listing:reveal"));
          });
        }
      });
    });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", enhance);
  else enhance();
})();
