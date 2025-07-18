using Franklin;
include("utils.jl")
println("Running deploy.jl!")

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
        if isnothing(type) || isnothing(vars)
          continue
        end

        # extract values
        page_title = isnothing(title) ? "Untitled" : title.captures[1]
        page_type  = type.captures[1]
        page_vars  = split(vars.captures[1], r"\s*,\s*")
        println("page_title = " * page_title)
        println("page_type = " * page_type)
        println("page_vars[1] = " * page_vars[1])

        try
          render_js(page_title, page_type, page_vars)
        catch e
          @warn "Error generating JS for $file" exception = e
        end
      end
    end
  end
end

eager_generate_js()
optimize()