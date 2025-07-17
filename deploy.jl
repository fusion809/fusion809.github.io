using Franklin

include("utils.jl")

function eager_generate_js()
  # Directory containing Markdown files
  md_dir = joinpath(@__DIR__, ".")  # assuming deploy.jl is in scripts/

  # Look for Markdown files
  for (root, _, files) in walkdir(md_dir)
    for file in files
      if endswith(file, ".md")
        md_path = joinpath(root, file)
        content = read(md_path, String)

        # Try to extract local vars like 'title', 'type', 'vars'
        title = match(r"title *= *\"([^\"]+)\"", content)
        type  = match(r"type *= *\"([^\"]+)\"", content)
        vars  = match(r"vars *= *\[([^\]]+)\]", content)

        if isnothing(type) && isnothing(vars)
          continue  # skip if there's no JS-relevant content
        end

        # Register the page variables for `locvar` to work (emulate context)
        setglobalvar!(Franklin.GLOBAL_VARS, Dict(
          "title" => isnothing(title) ? "Untitled" : title.captures[1],
          "type"  => isnothing(type)  ? ""        : type.captures[1],
          "vars"  => isnothing(vars)  ? nothing   : split(vars.captures[1], r"\s*,\s*")
        ))

        try
          hfun_render_js()  # will generate files if needed
        catch e
          @warn "Error calling hfun_render_js() on $file" exception=e
        end
      end
    end
  end
end

# Run file generation before optimizing
eager_generate_js()
optimize()