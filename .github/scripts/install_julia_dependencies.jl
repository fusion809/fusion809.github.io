using Pkg

function imported_package_name(node)
    if node isa Symbol
        return node == Symbol(".") ? nothing : String(node)
    elseif node isa Expr && node.head == :.
        isempty(node.args) && return nothing
        return imported_package_name(first(node.args))
    elseif node isa Expr && node.head == :(:)
        isempty(node.args) && return nothing
        return imported_package_name(first(node.args))
    elseif node isa Expr && node.head == :as
        isempty(node.args) && return nothing
        return imported_package_name(first(node.args))
    end
    return nothing
end

function collect_imported_packages!(packages::Set{String}, node)
    node isa Expr || return packages

    if node.head in (:using, :import)
        for imported in node.args
            name = imported_package_name(imported)
            isnothing(name) || push!(packages, name)
        end
        return packages
    end

    foreach(child -> collect_imported_packages!(packages, child), node.args)
    return packages
end

function scan_julia_packages(root::AbstractString = pwd())
    packages = Set{String}()
    skipped_dirs = Set([".git", "__site", "node_modules"])

    for (directory, subdirectories, filenames) in walkdir(root)
        filter!(name -> !(name in skipped_dirs), subdirectories)
        for filename in filenames
            (endswith(filename, ".jl") || endswith(filename, ".md")) || continue
            path = joinpath(directory, filename)
            source = read(path, String)

            if endswith(filename, ".jl")
                parsed = Meta.parseall(source; filename=path)
                collect_imported_packages!(packages, parsed)
            elseif endswith(filename, ".md")
                # Julia snippets in Markdown can import packages too, even if
                # their code is not stored in a separate .jl file.
                for block in eachmatch(r"(?ms)^```julia[^\n]*\n(.*?)^```", source)
                    parsed = Meta.parseall(block.captures[1]; filename=path)
                    collect_imported_packages!(packages, parsed)
                end
            end

            # Plots.jl loads optional backends by function call; these do not
            # appear as `using` or `import` statements for the backend package.
            occursin(r"\bplotlyjs\s*\(", source) && push!(packages, "PlotlyJS")
        end
    end

    stdlibs = Set(
        name for name in readdir(Sys.STDLIB)
        if isdir(joinpath(Sys.STDLIB, name))
    )
    filter!(name -> !(name in stdlibs) && name != "FunctionIntegrator", packages)
    return sort!(collect(packages))
end

function install_scanned_julia_dependencies(root::AbstractString = pwd())
    packages = scan_julia_packages(root)
    println("Julia packages discovered in repository scripts: ", join(packages, ", "))

    specs = [Pkg.PackageSpec(name=name) for name in packages]
    push!(specs, Pkg.PackageSpec(
        url="https://github.com/fusion809/FunctionIntegrator.jl",
        rev="v0.6.1",
    ))

    Pkg.add(specs)
    Pkg.instantiate()
    return nothing
end
