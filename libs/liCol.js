function splitList($list, n) {
    var i, k;
    var colHeight = Math.ceil($list.children().length / n)
    var colWidth = Math.floor(100 / n) + "%"

    for (i = 0; i < n; i++) {
        $list.append("<ul class='liCol'></ul>")
        for (k = 0; k < colHeight; k++) {
            $list.children("li").eq(0).appendTo(".liCol:last")          
        }   
    }

    $(".liCol").css("width", colWidth)
    $list.show() // list originally hidden to avoid displaying before ready
}