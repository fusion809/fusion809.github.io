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

function add_params_bef(params_val, params_desc, params_list, var_next, var, val, explanation)
  idx = findfirst(==(var_next), params_list);
  insert!(params_val, idx, val);
  insert!(params_desc, idx, explanation)
  insert!(params_list, idx, var)
end

function add_params_end(params_val, params_desc, params_list, var, val, explanation)
  push!(params_val, val);
  push!(params_desc, explanation)
  push!(params_list, var)
end

"""
  move_to_end!(val::Vector, list::Vector{String}, key::String)

Moves the element of list that equals key and its corresponding element of val to the end. 
"""
function move_to_end!(val::Vector, list::Vector{String}, key::String)
    i = findfirst(==(key), list)
    if isnothing(i)
        @warn "Key $key not found in params_list"
        return
    end
    push!(val,  popat!(val,  i))
    push!(list,  popat!(list,  i))
end

"""
  move3_to_end!(val::Vector{Number}, list::Vector{String}, desc::Vector{String}, key::String)

Moves the element of list that equals key and its corresponding element of val to the end. 
"""
function move3_to_end!(val::Vector{Real}, list::Vector{String}, desc::Vector{Any}, key::String)
    i = findfirst(==(key), list)
    if isnothing(i)
        @warn "Key $key not found in params_list"
        return
    end
    push!(val,  popat!(val,  i))
    push!(desc,  popat!(desc,  i))
    push!(list,  popat!(list,  i))
end

function calls_animate3D()
  path = locvar(:fd_rpath)
  dir = splitdir(path)[1]
  js_files = filter(f -> endswith(f, ".js"), readdir(dir; join=true))
  for jsf in js_files
    content = read(jsf, String)
    if occursin("animate3D", content)
      return true
    end
  end
  return false
end

function hfun_params_render()
  params = locvar("params");
  params_desc = [v.desc for v in values(params) if hasproperty(v, :desc)]
  if isempty(params_desc)
    params_desc = [];
  end
  params_val = [v.val for v in values(params) if hasproperty(v, :val)]
  params_list = collect(String.(keys(params)));
  for key in (:alpha, :gamma, :sigma, :rho, :beta, :lambda, :mu, :a, :b, :c)
    if haskey(params, key) && !hasproperty(getfield(params, key), :desc)
        push!(params_desc, "Problem parameter")
    end
  end
  if haskey(params, :tf) && !hasproperty(params.tf, :desc)
    push!(params_desc, "End time for the simulation in seconds (s)")
  end
  for key in (:x0, :y0, :z0, :S0, :E0, :I0, :R0)
    if haskey(params, key) && !hasproperty(getfield(params, key), :desc)
      key = replace(String(key), r"\d+" => "");
      push!(params_desc, "Initial \\($key\\) coordinate.")
    end
  end
  if !haskey(params, :epsilon)
    push!(params_list, "epsilon")
    push!(params_desc, "Error tolerance.")
    push!(params_val, 1e-11)
  end
  if haskey(params, :epsilon) && !hasproperty(params.epsilon, :desc)
    push!(params_desc, "Error tolerance.")
  end
  if !haskey(params, :tolType)
    push!(params_list, "tolType")
    push!(params_desc, "Tolerance type, can be either absolute (0) or relative (1).")
    push!(params_val, 0)
  end
  if haskey(params, :tolType) && !hasproperty(params.tolType, :desc)
    push!(params_desc, "Tolerance type, can be either absolute (0) or relative (1).")
  end
  add_params_end(params_val, params_desc, params_list, "hInitial", 0.1, "Initial step size.");
  add_params_end(params_val, params_desc, params_list, "hMin", 1e-8, "Minimum allowed step size.");
  add_params_end(params_val, params_desc, params_list, "Deltat", 0.2*(params.tf.val), "Time increment for skipping ahead in animation.");
  if !haskey(params, :t1)
    add_params_end(params_val, params_desc, params_list, "t1", 0.9*(params.tf.val), "Time you want to skip ahead to in animation when you press the skip button.");
  elseif (!hasproperty(params.t1, :desc))
    push!(params_desc, "Time you want to skip ahead to in animation when you press the skip button.")
    move_to_end!(params_val, params_list, "t1")
  end
  if !haskey(params, :Width)
    push!(params_list, "Width")
    push!(params_desc, "Width (in px) of Plotly windows used for plotting and animation below.")
    push!(params_val, 800)
  end
  if !haskey(params, :Height)
    push!(params_list, "Height")
    push!(params_desc, "Height (in px) of Plotly windows used for plotting and animation below.")
    push!(params_val, 600)
  end
  tScale_desc = "Proportion of animation time passed per real time. \\(t_{\\mathrm{Scale}}=1.0\\) means animation and real time match. \\(t_{\\mathrm{Scale}}<1.0\\) means the animation is going more slowly than real time. \\(t_{\\mathrm{Scale}}>1.0\\) means it is going more rapidly.";
  if !haskey(params, :tScale)
    push!(params_list, "tScale")
    push!(params_desc, tScale_desc)
    push!(params_val, 1.0)
  elseif !hasproperty(params.tScale, :desc)
    idx = findfirst(==("tScale"), params_list);
    insert!(params_desc, idx, tScale_desc)
    move3_to_end!(params_val, params_list, params_desc, "tScale");
  end
  if calls_animate3D()
    opacity_desc="Opacity of the lines in the 3D phase space animation. Customizable in case you need to tweak it in order to see the red dot marker."
    if !haskey(params, :Opacity)
      push!(params_list, "Opacity")
      push!(params_desc, opacity_desc)
      push!(params_val, 0.2)
    elseif !hasproperty(params.Opacity, :desc)
      idx = findfirst(==("Opacity"), params_list);
      insert!(params_desc, idx, opacity_desc)
      move3_to_end!(params_val, params_list, params_desc, "Opacity");
    end
  end

  params_list_len = length(params_list)
  params_desc_len = length(params_desc)
  params_val_len = length(params_val)
  if (params_list_len != params_desc_len || params_list_len != params_val_len)
    println("There is a mismatch in the lengths of params_list ($params_list_len), params_desc ($params_desc_len) or params_val ($params_val_len)! Exiting params_render().")
    println("params_val = $params_val")
    println("params_desc = $params_desc")
    println("params_list = $params_list")
    return
  end
  HTML = """
  <form name="requiredData">
    <table id="parameterForm">
      <caption style="font-size: 20px; font-weight: bold; text-align: left; border: 1px solid black; padding: 5px;">Simulation parameter form.</caption>
        <tr style="border: 0px solid black;">
          <th style="border: 1px solid black;">Parameter</th>
          <th style="border: 1px solid black;">Value</th>
          <th style="border: 1px solid black;">Explanation</th>
        </tr>
  """
  for p in 1:length(params_list)
    param_name = params_list[p];
    param_val = params_val[p];
    param_name = String(param_name);
    param_desc = params_desc[p];
    matches = collect(eachmatch(r"\d", param_name))
    numbers = [parse(Int64, m.match) for m in matches]
    if (occursin.(r"dtheta", param_name))
      if (length(numbers) == 2)
        first = numbers[1];
        second = numbers[2];
        param_name_latex="\\dot{\\theta}_{$first,$second}"
      elseif (length(numbers) == 1)
        first = numbers[1];
        param_name_latex="\\dot{\\theta}_{$first}"
      else
        param_name_latex="\\dot{\\theta}"
      end
    elseif (param_name == "Width" || param_name == "Height" || param_name == "Opacity")
      param_name_latex="\\mathrm{$param_name}"
    elseif (param_name == "tScale")
      param_name_latex="t_\\mathrm{Scale}"
    elseif (occursin.(r"theta", param_name))
      param_name_latex="\\theta"
      if (length(numbers) > 1)
        param_name_latex *= "_{"
        for i in 1:length(numbers)
          if i > 1
            param_name_latex *= ","
          end
          no = numbers[i];
          param_name_latex *= "$no"
        end
        param_name_latex *= "}"
      elseif (length(numbers) == 1)
        no = numbers[1]
        param_name_latex = "\\theta_{$no}"
      end
    elseif (occursin.(r"alpha", param_name))
      param_name_latex="\\alpha"
    elseif (occursin.(r"beta", param_name))
      param_name_latex="\\beta"
    elseif (occursin.(r"gamma", param_name))
      param_name_latex="\\gamma"
    elseif (occursin.(r"delta", param_name))
      param_name_latex="\\delta"
    elseif (occursin.(r"epsilon", param_name))
      param_name_latex="\\epsilon"
    elseif (occursin.(r"sigma", param_name))
      param_name_latex="\\sigma"
    elseif (occursin.(r"rho", param_name))
      param_name_latex="\\rho"
    elseif (occursin.(r"Deltat", param_name))
      param_name_latex="\\Delta t"
    elseif (occursin.(r"lambda", param_name))
      param_name_latex="\\Lambda"
    elseif (occursin.(r"mu", param_name))
      param_name_latex="\\mu" 
    elseif (occursin.(r"tf", param_name))
      param_name_latex="t_f"
    elseif (occursin.(r"dx0", param_name))
      param_name_latex="\\dot{x}_0"
    elseif (occursin.(r"dz0", param_name))
      param_name_latex="\\dot{z}_0"
    elseif (occursin.(r"dr", param_name))
      first = numbers[1];
      second = numbers[2];
      param_name_latex="\\dot{r}_{$first,$second}"
    elseif (occursin.(r"omega", param_name))
      param_name_latex="\\omega"
    elseif (param_name == "mb1")
    elseif (occursin.(r"[a-z][1-9][rb]", param_name))
      first = numbers[1];
      param_name_base = replace(param_name, r"\d+[rb]" => "");
      param_name_sub = replace(param_name, r"^[a-z]" => "");
      param_name_sub = replace(param_name_sub, r"\d+" => "");
      param_name_latex = "$param_name_base" * "_{" * "$first" * "$param_name_sub" * "}";
    elseif (occursin.(r"[a-z][rb]", param_name))
      param_name_base = replace(param_name, r"[rb]$" => "");
      param_name_sub = replace(param_name, r"^[a-z]" => "");
      param_name_latex = "$param_name_base" * "_{" * "$param_name_sub" * "}";
    elseif (length(numbers) == 1)
      first = numbers[1];
      param_name_cleaned = replace(param_name, r"\d+" => "");
      param_name_latex="$param_name_cleaned" * "_" * "$first";
    elseif (length(numbers) == 2)
      first = numbers[1];
      second = numbers[2];
      param_name_cleaned = replace(param_name, r"\d+" => "");
      param_name_latex="$param_name_cleaned" * "_" * "{$first,$second}";
    elseif (occursin.(r"hInitial", param_name))
      param_name_latex="h_{\\mathrm{Initial}}"
    elseif (occursin.(r"hMin", param_name))
      param_name_latex="h_{\\mathrm{Min}}"
    elseif (occursin.(r"tolType", param_name))
      param_name_latex="\\mathrm{Tol}_{\\mathrm{type}}"
    else
      param_name_latex=param_name
    end

    if (locvar("SP") !== nothing)
      HTML *= """
    <tr>
        <td><label for="$param_name">\\($param_name_latex\\):</label></td>
        <td><input type="Number"
                        id="$param_name"
                        name="$param_name"
                        value="$param_val"
                        onclick="periodCalc(thetaBounds(readInputs()))"
                        onchange="periodCalc(thetaBounds(readInputs()))"></td>
                    <td>$param_desc</td>
    </tr>
     """ 
    else
    HTML *= """
    <tr>
        <td><label for="$param_name">\\($param_name_latex\\):</label></td>
        <td><input type="Number"
                        id="$param_name"
                        name="$param_name"
                        value="$param_val"></td>
                    <td>$param_desc</td>
    </tr>
     """
    end
  end

  HTML *= """
    </table>
  </form>
  """
  return HTML
end

function getButtonVars()
  if (locvar("type") == "attractor")
    ids = ["tableOutputs", "phasePlotXYZ", "phasePlotXY", "phasePlotXZ", "phasePlotYZ", "timePlot", "phasePlotXYZ", "animation", "animation"]
    funcs = ["generateTable()", "removeTable()", "generate3DPhasePlot(solveProblem(RKF45, readInputs()))", "remove3DPhasePlot()", "generateXYPhasePlot(solveProblem(RKF45, readInputs()))", "removeXYPhasePlot()", "generateXZPhasePlot(solveProblem(RKF45, readInputs()))", "removeXZPhasePlot()", "generateYZPhasePlot(solveProblem(RKF45, readInputs()))", "removeYZPhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))", "removeTimePlot()", "generatePlots(readInputs())", "removePlots()", "generateAnimation()", "removeAnimation()","generateAllOutputs()", "removeAllOutputs()"]
    labels = ["Tabulate the solution", "Remove the table", "Generate 3D phase plot", "Remove \\(x\\), \\(y\\) and \\(z\\) phase plot", "Generate \\(x\\) and \\(y\\) phase plot", "Remove \\(x\\) and \\(y\\) phase plot", "Generate \\(x\\) and \\(z\\) phase plot", "Remove \\(x\\) and \\(z\\) phase plot", "Generate \\(y\\) and \\(z\\) phase plot", "Remove \\(y\\) and \\(z\\) phase plot", "Generate time plot for \\(x\\), \\(y\\) and \\(z\\)", "Remove time plot", "Generate all solution plots", "Remove all plots", "Generate an animation", "Remove animation", "Generate all outputs", "Remove all outputs"]
  elseif (locvar("type") == "2D")
    ids = ["tableOutputs", "phasePlot", "timePlot", "phasePlot", "animation", "animation"];
    funcs = ["generateTable()", "removeTable()", "generatePhasePlot(solveProblem(RKF45, readInputs()))", "removePhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))", "removeTimePlot()", "generatePlots(readInputs())", "removePlots()", "generateAnimation()", "removeAnimation()","generateAllOutputs()", "removeAllOutputs()"];
    labels = ["Tabulate the solution", "Remove the solution table", "Generate a phase plot", "Remove phase plot", "Generate a time plot", "Remove time plot", "Generate all solution plots", "Remove all solution plots", "Generate an animation", "Remove animation", "Generate all outputs", "Remove all outputs"];
  else
    ids = locvar("ids")
    funcs = locvar("funcs")
    labels = locvar("labels")
  end
  return ids, funcs, labels
end
function hfun_button_render()
  ids, funcs, labels = getButtonVars();
  if (2*length(ids) != length(labels) || 2*length(ids) != length(funcs))
    println("Labels or funcs are not twice the length of ids. Exiting button_render().")
    return;
  end

  HTML="""
  <table id="buttontable" style="border: 0px;">
  """

  for i=1:length(ids)
    id = ids[i];
    func1 = funcs[2*i-1];
    func2 = funcs[2*i];
    label1 = labels[2*i-1];
    label2 = labels[2*i];
    if (occursin(r"animation", id))
      id = replace(id, "animation" => "container")
    end
    HTML *= """
    <tr style="border: 0px;">
      <td style="border: 0px;">
        <a href="#$id">
          <button type="button" onclick="$func1">$label1</button>
        </a>
      </td>
      <td style="border: 0px;">
        <button type="button" onclick="$func2">$label2</button>
      </td>
    </tr>
    """;
  end
  HTML *= """</table>"""
  return HTML
end

function hfun_output_render()
  ids, funcs, labels = getButtonVars();
  HTML = """
  <table>
  """;
  ids = unique(ids);
  for i in 1:length(ids)
    if (occursin(r"animation*", ids[i]))
      idSuffix = replace(ids[i], r"^animation" => "")
      id = ids[i]
      HTML *= """
      <div id="container$idSuffix" class="button-bar">
      <h2 id="info$idSuffix" style="border: 0px; padding: 0; margin: 0; overflow: none;"></h2>
      <button id="toggleButton$idSuffix">Pause</button>
      <button id="restartButton$idSuffix">Restart</button>
      <button id="addTimeButton$idSuffix">Add time \\(\\Delta t\\) to animation.</button>
      <button id="skipToButton$idSuffix">Skip to \\(t_1\\)</button>
      <div id="animation$idSuffix"></div>
      </div>
      """
    else
      id = ids[i]
      HTML *= """
      <div id="container$id">
      <div id="$id"></div>
      </div>
      """
    end
  end
  HTML *= """
  </table>
  """

end

function resolve_inserts(html::String)::String
    insert_pattern = r"\{\{\s*insert\s+([a-zA-Z0-9_./\-]+)\s*\}\}"

    return replace(html, insert_pattern => s -> begin
        filename = match(insert_pattern, s).captures[1]
        filepath = joinpath(@__DIR__, "_layout", filename)
        if isfile(filepath)
            included = read(filepath, String)
            resolve_inserts(included)  # recursive resolution
        else
            println("$filepath was not found in resolve_inserts!")
        end
    end)
end

# Only need footer to include content pages list if it is not on the homepage.
function hfun_render_footlist()
    rpath = locvar("fd_rpath");
    if (rpath == "index.md")
        return """"""
    else
      mainfile = joinpath(@__DIR__, "_layout", "pages_list.html")
      content = read(mainfile, String)
      expanded = resolve_inserts(content)

      return """
      <button type='button' class='collapsible'>Show content pages list</button>
      <div class='content'><p>
      $expanded
      </p></div>
      """
    end
end

function hfun_render_js()
  
  if (locvar("type")=="attractor")
    artTitle = locvar("title")
    artTitle = replace(artTitle, r" solver"=>"")
    if (occursin(r"Chen", artTitle))
      title = "Chen system"
      view = "[0, 2, 0]";
    elseif (occursin(r"Lorenz", artTitle))
      title = "Lorenz system"
      view = "[2, 0, 0]";
    elseif (occursin(r"Rabinovich&ndash;Fabrikant", artTitle))
      title = "Rabinovich&ndash;Fabrikant system"
      view = "[2, 0, 0]";
    elseif (occursin(r"R&ouml;ssler", artTitle))
      title = "R&ouml;ssler system"
    else
      title = artTitle;
    end
    if (@isdefined(view))
      viewStr1 = ", {view: $view}"
      viewStr2 = "view: $view, "
    else
      viewStr1 = ""
      viewStr2 = ""
    end
    titleOneWord = split(title)[1];
    targetFile = "_libs/rendered/attractor_$titleOneWord.js"
    baseFile = "_libs/common/attractor.js"
    if !isfile(targetFile) || stat(baseFile).mtime > stat(targetFile).mtime
      # Read the content of the copied file
      cp(baseFile, targetFile, force=true);
      content = read(targetFile, String)

      # Replace characters or strings
      new_content = replace(content, "\$viewStr1" => "$viewStr1")
      new_content = replace(new_content, "\$viewStr2" => "$viewStr2")
      new_content = replace(new_content, "\$title" => "$title")

      # Write the modified content back to the file
      write(targetFile, new_content)
    end
    HTML = """
    <script src="/libs/rendered/attractor_$titleOneWord.js"></script>
    """
elseif (locvar("vars") !== nothing)
  varsList = locvar("vars")
  conds = ""
  for (i, name) in enumerate(varsList)
    conds *= name * "0"
    if i != length(varsList)
      conds *= ", "
    end
  end
  content = """
/** 
* Solve the problem using RKF45
*
* @param objectOfInputs An object containing all the problem parameters.
* @return               [t, vars]
*/
RKF45 = function(objectOfInputs) {
    // Extract initial conditions from object and enter it into RKF45Body
    var {$conds} = objectOfInputs;
    var vars0 = [[$conds]];
    var [t, vars] = RKF45Body(f, objectOfInputs, vars0);
    return [t, vars];
}
  """
  title = replace(locvar("title"), " " => "_")
  targetFile = "_libs/rendered/RKF45_$title.js"
  write(targetFile, content)
  HTML = """<script src="/libs/rendered/RKF45_$title.js"></script>"""
  else
    HTML = """"""
  end
  return HTML
end