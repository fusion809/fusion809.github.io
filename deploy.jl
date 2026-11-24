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

function fidello_saved_output_listings()
  fidello_dir = joinpath(@__DIR__, "fidello")
  output_names = Set{String}()
  source_paths = Dict{String,String}()
  output_pattern = r"\\output\{\./code/([^}]+)\}"
  listing_pattern = r"\\codelisting\{([^}]+)\}\{[^\n]*?\}\{\./code/([^}]+)\}\{([^}]+\.jl)\}"

  for (root, _, files) in walkdir(fidello_dir)
    for file in files
      endswith(file, ".md") || continue
      content = read(joinpath(root, file), String)
      union!(output_names, (m.captures[1] for m in eachmatch(output_pattern, content)))
      for m in eachmatch(listing_pattern, content)
        source_paths[m.captures[2]] = m.captures[3]
      end
    end
  end

  missing = setdiff(output_names, Set(keys(source_paths)))
  isempty(missing) || error("Fidello output commands have no matching codelisting source: $(join(sort!(collect(missing)), ", "))")

  return [("fidello", name, source_paths[name]) for name in sort!(collect(output_names))]
end

function generate_fidello_saved_outputs()
  listings = fidello_saved_output_listings()

  for (page_path, output_name, source_path) in listings
    source_file = joinpath(@__DIR__, source_path)
    isfile(source_file) || error("Required Fidello listing source is missing: $source_file")

    output_dir = joinpath(@__DIR__, "__site", "assets", page_path, "code", "output")
    output_file = joinpath(output_dir, output_name * ".out")
    mkpath(output_dir)

    open(output_file, "w") do output
      redirect_stdout(output) do
        listing_module = Module(gensym(:FidelloListing))
        Base.include_string(listing_module, read(source_file, String), source_file)
      end
    end

    filesize(output_file) > 0 || error("Fidello listing produced empty output: $source_file")

    if output_name == "chebyshev"
      source = read(source_file, String)
      n_match = match(r"(?m)^\s*n\s*=\s*(\d+)\s*$", source)
      isnothing(n_match) && error("Could not determine Chebyshev grid size from $source_file")
      expected_rows = parse(Int, n_match.captures[1]) + 1
      output_header = first(eachline(output_file))
      startswith(output_header, "$(expected_rows)×5 DataFrame") ||
        error("Chebyshev output is stale: expected $expected_rows rows for n=$(expected_rows - 1), got '$output_header'")
    end
  end

  return [name for (_, name, _) in listings]
end

function append_reference_lists(site_dir)
  citation_pattern = r"<span data-bib-key=\"([^\"]+)\" data-bib-entry=\"([^\"]+)\">([\s\S]*?)</span>"
  preview_script = r"<script>\(function\(\)\{\s*if \(window\.__bibReferenceListInstalled[\s\S]*?</script>"
  page_locator = r"(<a class=\"citation\" href=\"#[^\"]+\">)([^<]*)(</a>\s*\))\s*(?:&#123;|\{)([^{}&]+)(?:&#125;|\})"
  processed = 0
  for (root, _, files) in walkdir(site_dir)
    for file in files
      endswith(file, ".html") || continue
      path = joinpath(root, file)
      html = read(path, String)
      html = replace(html, preview_script => "")
      if occursin("<section class=\"reference-list\"", html)
        html = replace(html, "<h2 id=\"reference-list-title\">Reference list</h2>" => "<h1 id=\"reference-list-title\">Reference list</h1>")
        toc_link = "<li><a href=\"#reference-list-title\">Reference list</a></li>"
        if !occursin("href=\"#reference-list-title\"", html)
          html = replace(html, "</ol></div>" => toc_link * "</ol></div>"; count=1)
        end
        write(path, html)
        processed += 1
        continue
      end
      matches = collect(eachmatch(citation_pattern, html))
      isempty(matches) && continue
      cited = Dict{String,String}()
      reference_id(key) = "reference-" * bytes2hex(codeunits(key))
      for m in matches
        key, entry = String(base64decode(m.captures[1])), String(base64decode(m.captures[2]))
        cited[key] = entry
      end
      # Link each in-text citation to its one bibliography entry.
      html = replace(html, citation_pattern => m -> begin
        captures = match(citation_pattern, m).captures
        key = String(base64decode(captures[1]))
        "<a class=\"citation\" href=\"#$(reference_id(key))\">$(captures[3])</a>"
      end)
      # Franklin escapes literal braces in prose as entities, and the locator
      # follows the parenthesized citation in rendered HTML.
      html = replace(html, page_locator => m -> begin
        captures = match(page_locator, m).captures
        captures[1] * captures[2] * ", p. " * strip(captures[4]) * captures[3]
      end)
      sorted = sort(collect(cited); by=pair -> lowercase(strip(replace(last(pair), r"<[^>]+>" => ""))))
      items = join("<li id=\"$(reference_id(first(pair)))\">" * last(pair) * "</li>" for pair in sorted)
      section = "<section class=\"reference-list\" aria-labelledby=\"reference-list-title\"><h1 id=\"reference-list-title\">Reference list</h1><ol>$(items)</ol></section>"
      toc_link = "<li><a href=\"#reference-list-title\">Reference list</a></li>"
      html = replace(html, "</ol></div>" => toc_link * "</ol></div>"; count=1)
      if occursin("<div class=\"page-foot\">", html)
        html = replace(html, "<div class=\"page-foot\">" => section * "<div class=\"page-foot\">"; count=1)
      elseif occursin("<!-- CONTENT ENDS HERE -->", html)
        html = replace(html, "<!-- CONTENT ENDS HERE -->" => section * "<!-- CONTENT ENDS HERE -->"; count=1)
      elseif occursin("</article>", html)
        html = replace(html, "</article>" => section * "</article>"; count=1)
      else
        html = replace(html, "</body>" => section * "</body>"; count=1)
      end
      write(path, html)
      processed += 1
    end
  end
  return processed
end

eager_generate_js()
fidello_output_names = generate_fidello_saved_outputs()
fidello_source = joinpath(@__DIR__, "fidello", "index.md")
fidello_output = joinpath(@__DIR__, "__site", "fidello", "index.html")
isfile(fidello_source) || error("Required Fidello source page is missing: $fidello_source")
build_started_at = time()
optimize(prerender=false, suppress_errors=false)

if !isfile(fidello_output)
  error("Franklin did not generate the required Fidello page: $fidello_output")
end
if stat(fidello_output).mtime < build_started_at
  error("Franklin left a stale Fidello page instead of rebuilding it: $fidello_output")
end

fidello_html = read(fidello_output, String)
for name in fidello_output_names
  output_file = joinpath(@__DIR__, "__site", "assets", "fidello", "code", "output", name * ".out")
  isfile(output_file) && filesize(output_file) > 0 ||
    error("Required Fidello listing output is missing or empty: $output_file")

  output_error = string('`', raw"\output{./code/", name,
                        "}`: could not find the relevant output file.")
  occursin(output_error, fidello_html) &&
    error("Fidello page could not resolve \\output for '$name': $fidello_output")
end
reference_list_pages = append_reference_lists(joinpath(@__DIR__, "__site"))
println("Added reference lists to $reference_list_pages generated page(s).")
