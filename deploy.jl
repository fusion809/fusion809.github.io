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
        else
          ids = String[]
          funcs = String[]
        end

        try
          render_js(page_title, page_type, page_vars)
          isempty(funcs) || render_rmPlot(funcs, ids, page_title)
        catch e
          @warn "Error generating JS for $file" exception = e
        end
      end
    end
  end
end

function generate_fidello_saved_outputs()
  listings = [
    ("fidello", "fidello/code/findiff.jl"),
    ("fidello", "fidello/code/findiff4th.jl"),
  ]

  for (page_path, source_path) in listings
    source_file = joinpath(@__DIR__, source_path)
    isfile(source_file) || error("Required Fidello listing source is missing: $source_file")

    source_name = splitext(basename(source_file))[1]
    output_dir = joinpath(@__DIR__, "__site", "assets", page_path, "code", "output")
    output_file = joinpath(output_dir, source_name * ".out")
    mkpath(output_dir)

    open(output_file, "w") do output
      redirect_stdout(output) do
        include_string(Main, read(source_file, String), source_file)
      end
    end

    filesize(output_file) > 0 || error("Fidello listing produced empty output: $source_file")
  end
end

eager_generate_js()
generate_fidello_saved_outputs()
fidello_source = joinpath(@__DIR__, "fidello", "index.md")
fidello_output = joinpath(@__DIR__, "__site", "fidello", "index.html")
isfile(fidello_source) || error("Required Fidello source page is missing: $fidello_source")
build_started_at = time()
optimize(prerender=false)

if !isfile(fidello_output)
  error("Franklin did not generate the required Fidello page: $fidello_output")
end
if stat(fidello_output).mtime < build_started_at
  error("Franklin left a stale Fidello page instead of rebuilding it: $fidello_output")
end

fidello_html = read(fidello_output, String)
if occursin("could not find the relevant output file", fidello_html)
  error("Fidello page contains an unresolved \\output command: $fidello_output")
end

for name in ("findiff", "findiff4th")
  output_file = joinpath(@__DIR__, "__site", "assets", "fidello", "code", "output", name * ".out")
  isfile(output_file) && filesize(output_file) > 0 ||
    error("Required Fidello listing output is missing or empty: $output_file")
end
