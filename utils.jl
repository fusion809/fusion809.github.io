function hfun_bar(vname)
  val = Meta.parse(vname[1])
  return round(sqrt(val), digits=2)
end

function hfun_m1fill(vname)
  var = vname[1]
  return pagevar("index", var)
end

function lx_baz(com, _)
  # keep this first line
  brace_content = Franklin.content(com.braces[1]) # input string
  # do whatever you want here
  return uppercase(brace_content)
end

 function render_form(params::Vector{String})
     inputs = join(["<tr><td><label for=\"$p\">\\($p\\):</label></td>\n<td><input name=\"$p\" id=\"$p\" name=\"$p\" value=\"{{$p}}\"/></td>\n<td>{{ {$p}_explanation }} </td>\n</tr>" for p in params], "\n")
     return "<form>\n$inputs\n</form>"
end