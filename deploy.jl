using Franklin;
include("utils.jl")
println("Running deploy.jl!")

function get_hassim_value(content::AbstractString)
    m = match(r"@def\s+hassim\s*=\s*(true|false)", content)
    if m === nothing
        return false;
    else
        return m.captures[1] == "true"
    end
end

function eager_generate_js()
  md_dir = joinpath(@__DIR__, ".")

  for (root, _, files) in walkdir(md_dir)
    for file in files
      if endswith(file, ".md")
        md_path = joinpath(root, file)
        content = read(md_path, String)
        println("md_path = " * md_path)
        title = match(r"title *= *\"([^\"]+)\"", content)
        type  = match(r"type *= *\"([^\"]+)\"", content)
        vars  = match(r"vars *= *\[([^\]]+)\]", content)
        hassim = get_hassim_value(content);
        if isnothing(type) || isnothing(vars)
          continue
        end

        # extract values
        page_title = isnothing(title) ? "Untitled" : title.captures[1]
        page_type  = type.captures[1]
        page_vars  = split(vars.captures[1], r"\s*,\s*")
        if (page_type == "2D")
          ids = ["tableOutputs", "phasePlot", "timePlot", "phasePlot", "animation", "animation"];
          funcs = ["generateTable()", "removeTable()", "generatePhasePlot(solveProblem(RKF45, readInputs()))", "removePhasePlot()", "generateTimePlot(solveProblem(RKF45, readInputs()))", "removeTimePlot()", "generatePlots(readInputs())", "removePlots()", "generateAnimation()", "removeAnimation()","generateAllOutputs()", "removeAllOutputs()"];
        elseif hassim
          ids   = split(match(r"ids *= *\[([^\]]+)\]", content).captures[1], r"\s*,\s*")
          funcs   = split(match(r"funcs *= *\[([^\]]+)\]", content).captures[1], r"\s*,\s*")
        end

        try
          render_js(page_title, page_type, page_vars)
          render_rmPlot(funcs, ids, page_title)
        catch e
          @warn "Error generating JS for $file" exception = e
        end
      end
    end
  end
end

eager_generate_js()
optimize()