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

function hfun_params_render()
  params = locvar("params");
  params_desc = [v.desc for v in values(params) if hasproperty(v, :desc)]
  if isempty(params_desc)
    params_desc = [];
  end
  params_el = collect(values(params));
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
    push!(params_desc, "Absolute error tolerance.")
    push!(params_val, 1e-11)
  end
  if haskey(params, :epsilon) && !hasproperty(params.epsilon, :desc)
    push!(params_desc, "Absolute error tolerance.")
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
  if !haskey(params, :width)
    push!(params_list, "Width")
    push!(params_desc, "Width (in px) of Plotly windows used for plotting and animation below.")
    push!(params_val, 800)
  end
  if !haskey(params, :height)
    push!(params_list, "Height")
    push!(params_desc, "Height (in px) of Plotly windows used for plotting and animation below.")
    push!(params_val, 600)
  end
  if !haskey(params, :delay)
    push!(params_list, "Delay")
    push!(params_desc, "Proportion of animation time passed per real time. Delay=1.0 means animation and real time match. Delay<1 means the animation is going more slowly than real time. Delay>1.0 means it is going more rapidly.")
    push!(params_val, 1.0)
  end

  if (length(params_list) != length(params_desc) || length(params_list) != length(params_val))
    println("There is a mismatch in the lengths of params_list, params_desc or params_val! Exiting params_render().")
    return
  end
  HTML = """"""
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
    elseif (param_name == "Width" || param_name == "Height" || param_name=="Delay")
      param_name_latex="\\mathrm{$param_name}"
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
    else
      param_name_latex=param_name
    end

    if (locvar("SP") !== "nothing")
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

  return HTML
end

function getButtonVars()
  if (locvar("type") == "attractor")
    ids = ["tableOutputs", "phasePlotXYZ", "phasePlotXY", "phasePlotXZ", "phasePlotYZ", "timePlot", "phasePlotXYZ", "animation"]
    funcs = ["generateTable()", "removeTable()", "generate3DPhasePlot(solveProblem(RKF45, readInputs()))", "remove3DPhasePlot()", "generateXYPhasePlot(solveProblem(RKF45, readInputs()))", "removeXYPhasePlot()", "generateXZPhasePlot(solveProblem(RKF45, readInputs()))", "removeXZPhasePlot()", "generateYZPhasePlot(solveProblem(RKF45, readInputs()))", "removeYZPhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))", "removeTimePlot()", "generatePlots(readInputs())", "removePlots()", "generateAnimation()", "removeAnimation()"]
    labels = ["Tabulate the solution", "Remove the table", "Generate 3D phase plot", "Remove \\(x\\), \\(y\\) and \\(z\\) phase plot", "Generate \\(x\\) and \\(y\\) phase plot", "Remove \\(x\\) and \\(y\\) phase plot", "Generate \\(x\\) and \\(z\\) phase plot", "Remove \\(x\\) and \\(z\\) phase plot", "Generate \\(y\\) and \\(z\\) phase plot", "Remove \\(y\\) and \\(z\\) phase plot", "Generate time plot for \\(x\\), \\(y\\) and \\(z\\)", "Remove time plot", "Generate all solution plots", "Remove all plots", "Generate an animation", "Remove animation"]
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
  animation_ids = [];
  other_ids = [];
  HTML = """
  <table id="tableOutputs">
  """;
  ids = unique(ids);
  for i in 1:length(ids)
    if (occursin(r"animation*", ids[i]))
      idSuffix = replace(ids[i], r"^animation" => "")
      id = ids[i]
      HTML *= """
      <button id="toggleButton$idSuffix">Pause</button>
      <button id="restartButton$idSuffix">Restart</button>
      <button id="addTimeButton$idSuffix">Add time \\(\\Delta t\\) to animation.</button>
      <button id="skipToButton$idSuffix">Skip to \\(t_1\\)</button>
      <div id="animation$idSuffix"></div>
      """
    else
      id = ids[i]
      HTML *= """
      <div id="$id"></div>
      """
    end
  end
  HTML *= """
  </table>
  """

end