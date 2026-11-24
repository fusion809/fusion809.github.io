using SHA
import Franklin
using Base64

# Keep GR headless during Franklin's build-time evaluation of code listings.
ENV["GKSwstype"] = "100"
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

function lx_figref(com, _)
  brace_content = strip(Franklin.content(com.braces[1]))
  return "~~~<a class=\"figref\" href=\"#$brace_content\"></a>~~~"
end

function lx_Figref(com, _)
  brace_content = strip(Franklin.content(com.braces[1]))
  return "~~~<a class=\"figref-full\" href=\"#$brace_content\"></a>~~~"
end

# BibTeX-style citations.  The citation macro emits a compact build-time
# marker; deploy.jl turns these into APA author-date citations and appends the
# page's cited reference list after Franklin has rendered the page.
function bibtex_entries(source::AbstractString)
  entries = Dict{String,Dict{String,String}}()
  for m in eachmatch(r"(?is)@([a-z]+)\s*([({])\s*([^,\s]+)\s*,", source)
    kind, key = lowercase(m.captures[1]), strip(m.captures[3])
    start = m.offset + ncodeunits(m.match)
    openchar = only(m.captures[2])
    closechar = openchar == '{' ? '}' : ')'
    depth, quoted, escaped, stop = 1, false, false, 0
    for pos in eachindex(source)
      pos < start && continue
      c = source[pos]
      if escaped
        escaped = false
      elseif c == '\\'
        escaped = true
      elseif c == '"'
        quoted = !quoted
      elseif !quoted && c == openchar
        depth += 1
      elseif !quoted && c == closechar
        depth -= 1
        if depth == 0
          stop = pos
          break
        end
      end
    end
    stop == 0 && error("Unclosed BibTeX entry '$key'")
    body = source[start:prevind(source, stop)]
    fields = Dict{String,String}()
    for fm in eachmatch(r"(?is)([\w-]+)\s*=\s*(\{(?:[^{}]|\{[^{}]*\})*\}|\"(?:[^\"\\]|\\.)*\"|[^,]+)", body)
      value = strip(fm.captures[2])
      value = (startswith(value, "{") && endswith(value, "}")) ||
              (startswith(value, "\"") && endswith(value, "\"")) ? chop(value; head=1, tail=1) : value
      value = replace(value, r"\\[a-zA-Z]+\s*" => "", "~" => " ", r"\s+" => " ")
      value = replace(value, "{" => "", "}" => "")
      fields[lowercase(fm.captures[1])] = strip(value)
    end
    entries[key] = merge(Dict("_type" => kind), fields)
  end
  return entries
end

html_escape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;", "'" => "&#39;")

function bib_authors(entry)
  names = strip.(split(get(entry, "author", get(entry, "editor", "Unknown author")), r"\s+and\s+"))
  surnames, formatted = String[], String[]
  for name in names
    pieces = strip.(split(name, ","))
    if length(pieces) > 1
      surname, given = first(pieces), join(pieces[2:end], " ")
    else
      words = split(strip(replace(name, "{" => "", "}" => "")))
      surname = length(words) > 1 ? last(words) : join(words, " ")
      given = length(words) > 1 ? join(words[1:end-1], " ") : ""
    end
    push!(surnames, surname)
    initials = join((uppercase(first(strip(w))) * "." for w in split(given) if !isempty(strip(w))), " ")
    push!(formatted, isempty(initials) ? surname : surname * ", " * initials)
  end
  intext = length(surnames) == 1 ? surnames[1] : length(surnames) == 2 ? join(surnames, " & ") : surnames[1] * " et al."
  list = length(formatted) > 20 ? join(formatted[1:19], ", ") * ", … " * last(formatted) : length(formatted) > 1 ? join(formatted[1:end-1], ", ") * ", & " * last(formatted) : only(formatted)
  return intext, list
end

function bib_apa_entry(key, entry)
  _, authors = bib_authors(entry)
  year = get(entry, "year", "n.d.")
  title = get(entry, "title", "Untitled")
  kind = get(entry, "_type", "")
  italic = s -> "<em>" * html_escape(s) * "</em>"
  linked_title = function(text, url; is_italic=false)
    rendered = is_italic ? italic(text) : html_escape(text)
    isempty(url) && return rendered
    return "<a href=\"$(html_escape(url))\">$(rendered)</a>"
  end
  title_html = linked_title(title, get(entry, "url", ""); is_italic=true)
  if kind == "article"
    journal = get(entry, "journal", "")
    vol = get(entry, "volume", "")
    issue = get(entry, "number", "")
    pages = get(entry, "pages", "")
    journalpart = isempty(journal) ? "" : italic(journal) * (isempty(vol) ? "" : ", " * italic(vol)) * (isempty(issue) ? "" : "(" * html_escape(issue) * ")") * (isempty(pages) ? "" : ", " * html_escape(replace(pages, "--" => "–"))) * ". "
    doi = get(entry, "doi", "")
    doi_url = isempty(doi) ? "" : "https://doi.org/" * replace(doi, r"^https?://(dx\.)?doi.org/" => "")
    return "$(html_escape(authors)) ($(html_escape(year))). $(linked_title(title, get(entry, "url", ""))). $(journalpart)" * (isempty(doi_url) ? "" : "<a href=\"$(html_escape(doi_url))\">$(html_escape(doi_url))</a>")
  elseif kind == "book"
    publisher = get(entry, "publisher", "")
    return "$(html_escape(authors)) ($(html_escape(year))). $(title_html)." * (isempty(publisher) ? "" : " " * html_escape(publisher) * ".")
  elseif kind in ("inbook", "incollection")
    chapter = kind == "incollection" ? get(entry, "title", "Untitled") : get(entry, "chapter", "")
    booktitle = get(entry, "booktitle", get(entry, "title", ""))
    booktitle_html = italic(booktitle)
    editor = get(entry, "editor", "")
    editor_text = isempty(editor) ? "" : html_escape(bib_authors(Dict("author" => editor))[2]) * (occursin(r"\s+and\s+", editor) ? " (Eds.), " : " (Ed.), ")
    pages = get(entry, "pages", "")
    page_text = isempty(pages) ? "" : " (pp. $(html_escape(replace(pages, "--" => "–"))))."
    publisher = get(entry, "publisher", "")
    chapter_text = isempty(chapter) ? "" : linked_title(chapter, get(entry, "url", "")) * ". "
    in_book = isempty(booktitle) ? "" : "In " * editor_text * booktitle_html * (isempty(page_text) ? "." : page_text)
    return "$(html_escape(authors)) ($(html_escape(year))). $(chapter_text)$(in_book)" * (isempty(publisher) ? "" : " " * html_escape(publisher) * ".")
  elseif kind in ("inproceedings", "conference")
    booktitle = get(entry, "booktitle", "")
    publisher = get(entry, "publisher", "")
    return "$(html_escape(authors)) ($(html_escape(year))). $(linked_title(title, get(entry, "url", "")))." * (isempty(booktitle) ? "" : " In " * italic(booktitle) * ".") * (isempty(publisher) ? "" : " " * html_escape(publisher) * ".")
  end
  return "$(html_escape(authors)) ($(html_escape(year))). $(linked_title(title, get(entry, "url", "")))."
end

function inherited_bib_setting(dir, name)
  project_root = normpath(@__DIR__)
  current = normpath(dir)
  while startswith(current, project_root)
    index_file = joinpath(current, "index.md")
    if isfile(index_file)
      pattern = Regex("(?m)^\\s*@def\\s+" * name * "\\s*=\\s*([^\\n;]+)")
      m = match(pattern, read(index_file, String))
      if m !== nothing
        value = strip(m.captures[1])
        if (startswith(value, "\"") && endswith(value, "\"")) ||
           (startswith(value, "'") && endswith(value, "'"))
          value = chop(value; head=1, tail=1)
        end
        return value, current
      end
    end
    current == project_root && break
    parent = dirname(current)
    parent == current && break
    current = parent
  end
  return nothing, normpath(dir)
end

function bib_citation(com, _)
  keys = strip.(split(Franklin.content(com.braces[1]), ","))
  bibfile = Franklin.locvar("ref")
  page_path = normpath(joinpath(@__DIR__, String(Franklin.locvar(:fd_rpath))))
  page_dir = dirname(page_path)
  inherited_ref, inherited_dir = inherited_bib_setting(page_dir, "ref")
  isnothing(bibfile) && (bibfile = inherited_ref)
  isnothing(bibfile) && throw(ArgumentError("\\cite requires @def ref=\"file.bib\" in this page or its directory index.md"))
  style = Franklin.locvar("style")
  if isnothing(style)
    inherited_style, _ = inherited_bib_setting(page_dir, "style")
    style = something(inherited_style, "apa")
  end
  style = lowercase(string(style))
  style == "apa" || throw(ArgumentError("Unsupported citation style '$style'; currently supported: apa"))
  root = isnothing(Franklin.locvar("ref")) ? inherited_dir : page_dir
  path = normpath(joinpath(root, string(bibfile)))
  isfile(path) || throw(ArgumentError("BibTeX file not found: $path"))
  entries = bibtex_entries(read(path, String))
  citations = String[]
  for key in keys
    haskey(entries, key) || throw(ArgumentError("Citation key '$key' not found in $bibfile"))
    author, _ = bib_authors(entries[key])
    year = get(entries[key], "year", "n.d.")
    push!(citations, "<span data-bib-key=\"$(base64encode(key))\" data-bib-entry=\"$(base64encode(bib_apa_entry(key, entries[key])))\">$(html_escape(author)), $(html_escape(year))</span>")
  end
  preview_script = raw"""
  <script>(function(){
    if (window.__bibReferenceListInstalled) return;
    window.__bibReferenceListInstalled = true;
    function decode(value){ return new TextDecoder().decode(Uint8Array.from(atob(value), function(c){ return c.charCodeAt(0); })); }
    function referenceId(key){ return 'reference-' + Array.from(new TextEncoder().encode(key), function(b){ return b.toString(16).padStart(2, '0'); }).join(''); }
    function renderBibReferences(){
      var markers = Array.from(document.querySelectorAll('.franklin-content span[data-bib-key]'));
      var entries = new Map();
      markers.forEach(function(marker){
        var key = decode(marker.getAttribute('data-bib-key'));
        entries.set(key, decode(marker.getAttribute('data-bib-entry')));
        var label = marker.innerHTML;
        // Franklin inserts this preview script between the citation and any
        // following page locator, so inspect the text immediately after it.
        var preview = marker.nextSibling;
        while (preview && !(preview.nodeType === Node.ELEMENT_NODE && preview.tagName === 'SCRIPT' && preview.textContent.indexOf('__bibReferenceListInstalled') >= 0)) {
          preview = preview.nextSibling;
        }
        var following = preview ? preview.nextSibling : null;
        while (following && following.nodeType === Node.TEXT_NODE && !following.textContent.trim()) {
          following = following.nextSibling;
        }
        if (following && following.nodeType === Node.TEXT_NODE) {
          var locator = following.textContent.match(/^\s*\{([^{}]+)\}/);
          if (locator) {
            label += ', p. ' + locator[1].trim();
            following.textContent = following.textContent.replace(locator[0], '');
          }
        }
        marker.outerHTML = '<a class="citation" href="#' + referenceId(key) + '">' + label + '</a>';
      });
      if (!entries.size || document.querySelector('.reference-list')) return;
      var section = document.createElement('section');
      section.className = 'reference-list';
      section.setAttribute('aria-labelledby', 'reference-list-title');
      section.innerHTML = '<h1 id="reference-list-title">Reference list</h1><ul style="list-style:none;margin:0;padding-left:0;"></ul>';
      Array.from(entries.entries()).sort(function(a,b){ return a[1].replace(/<[^>]+>/g,'').localeCompare(b[1].replace(/<[^>]+>/g,'')); }).forEach(function(item){
        var li = document.createElement('li'); li.id = referenceId(item[0]); li.style.paddingLeft = '2em'; li.style.textIndent = '-2em'; li.innerHTML = item[1]; section.querySelector('ul').appendChild(li);
      });
      var pageFoot = document.querySelector('.page-foot');
      if (pageFoot && pageFoot.parentNode) pageFoot.parentNode.insertBefore(section, pageFoot);
      else (document.querySelector('.franklin-content') || document.body).appendChild(section);
      var toc = document.querySelector('.franklin-toc > ol');
      if (toc && !toc.querySelector('a[href="#reference-list-title"]')) {
        var tocItem = document.createElement('li');
        tocItem.innerHTML = '<a href="#reference-list-title">Reference list</a>';
        toc.appendChild(tocItem);
      }
      var nav = document.querySelector('#toc-list');
      if (nav && !nav.querySelector('a[href="#reference-list-title"]')) {
        var navItem = document.createElement('li');
        navItem.innerHTML = '<a href="#reference-list-title">Reference list</a>';
        nav.appendChild(navItem);
      }
    }
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', renderBibReferences);
    else renderBibReferences();
  })();</script>
  """
  return "(" * join(citations, "; ") * ")" * preview_script
end

@eval Franklin function lx_cite(com::LxCom, lxdefs)
  return getproperty(utils_module(), :bib_citation)(com, lxdefs)
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
    elseif (occursin.(r"kappa", param_name))
      param_name_latex="\\kappa"
    elseif (occursin.(r"rho", param_name))
      param_name_latex="\\rho"
    elseif (occursin.(r"Deltat", param_name))
      param_name_latex="\\Delta t"
    elseif (occursin.(r"lambda", param_name))
      param_name_latex="\\lambda"
    elseif (occursin.(r"mu", param_name))
      param_name_latex="\\mu" 
    elseif (occursin.(r"tf", param_name))
      param_name_latex="t_f"
    elseif (occursin.(r"dx0", param_name))
      param_name_latex="\\dot{x}_0"
    elseif (occursin.(r"dy0", param_name))
      param_name_latex="\\dot{y}_0"
    elseif (occursin.(r"dz0", param_name))
      param_name_latex="\\dot{z}_0"
    elseif (occursin.(r"dr", param_name))
      first = numbers[1];
      if (length(numbers) > 1)
        second = numbers[2];
        param_name_latex="\\dot{r}_{$first,$second}"
      else
        param_name_latex="\\dot{r}_{$first}"
      end
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

function latest_input_time()
    # Get base files (utils and base JS files)
    base_files = [
        "utils.jl",
        "_libs/common/RKF45_base.js",
        "_libs/common/attractor.js"
    ]

    # Get all Markdown files recursively in the pwd
    md_dir = joinpath(@__DIR__, ".")
    md_files = [];
    for (root, _, files) in walkdir(md_dir)
      for file in files
        if endswith(file, ".md")
          push!(md_files, root * "/" * file)
        end
      end
    end

    # Only keep those defining title, type, and vars
    filtered_md_files = filter(md_files) do f
        content = read(f, String)
        occursin("@def title", content) &&
        occursin("@def type", content) &&
        occursin("@def vars", content)
    end

    # Combine with base files
    all_inputs = vcat(base_files, filtered_md_files)

    # Return latest mtime
    return maximum(stat(f).mtime for f in all_inputs if isfile(f))
end

function render_js(pageTitle, type, vars)
  title = replace(pageTitle, " " => "_")
  fileTitle = replace(pageTitle, " " => "_")
  if (type=="attractor")
    artTitle = pageTitle
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
    targetFile = "_libs/rendered/attractor_$fileTitle.js"
    baseFile = "_libs/common/attractor.js"
  elseif (vars !== nothing)
    varsList = vars
    conds = ""
    for (i, name) in enumerate(varsList)
      conds *= string(name, "0")
      if i != length(varsList)
        conds *= ", "
      end
    end
    targetFile = "_libs/rendered/RKF45_$fileTitle.js"
    baseFile = "_libs/common/RKF45_base.js"
  else
    return;
  end
  if !isdir("_libs/rendered")
    mkpath("_libs/rendered/")
  end
  if !isdir("__site/libs/rendered/")
    mkpath("__site/libs/rendered/")
  end
  if !isfile(targetFile) || stat(baseFile).mtime > stat(targetFile).mtime
      println("Copying $baseFile to $targetFile...")
      cp(baseFile, targetFile, force=true)
      content = read(targetFile, String)

      # Replace characters or strings
      if (@isdefined(conds))
        content = replace(content, "\$conds" => "$conds")
      end

      if @isdefined(viewStr1)
        content = replace(content, "\$viewStr1" => "$viewStr1")
        content = replace(content, "\$viewStr2" => "$viewStr2")
      end

      content = replace(content, "\$title" => "$title")
      if (!occursin("attractor", baseFile))
        content = replace(content, r"\"" => "")
      end
      new_hash = sha256(content)
      old_hash = isfile(targetFile) ? sha256(read(targetFile, String)) : ""

      # Write the modified content back to the file
      if (new_hash != old_hash)
        write(targetFile, content)
        cp(targetFile, "__site/libs/rendered/RKF45_$fileTitle.js", force=true)
      end
    end
end


function render_rmPlot(funcs, ids, title)
  title = replace(title, " " => "_");
  filepath1 = "_libs/rendered/rmPlot_$title.js"
  filepath2 = "__site/libs/rendered/rmPlot_$title.js"
  if isfile(filepath1) && stat(filepath1).mtime > latest_input_time() &&
    isfile(filepath2) && stat(filepath2).mtime > latest_input_time()
      return
  end
  src = """""";
  for k in [r"[pP]lot", r"[aA]nimation"]
    rmFuncs = filter(f -> occursin("remove", f) && occursin(k, f), funcs)
    plotIDs = filter(f -> occursin(k, f), ids)
    if (length(rmFuncs) < length(plotIDs))
      plotIDs = plotIDs[1:length(rmFuncs)]
    end
    for i in eachindex(rmFuncs)
      rmFunc = rmFuncs[i]
      plotID = plotIDs[i]
      if !occursin(k * "s", rmFunc)
        src *= """
        function $rmFunc {
          rmPlot("$plotID")
        }
        """
      elseif occursin("removePendulumPlots", rmFunc)
        src *= """
        function $rmFunc {
        """
        for j in eachindex(rmFuncs)
          rmFuncj = rmFuncs[j]
          if (i != j) && occursin(r"[pP]endulum", rmFuncj)
            println("rmFuncj = " * rmFuncj)
            src *= """
              $rmFuncj
            """
          end
        end
        src *= """
        }
        """
      elseif occursin("removePlots", rmFunc)
        src *= """
        function $rmFunc {
        """
        for j in eachindex(rmFuncs)
          rmFuncj = rmFuncs[j]
          if (i != j) && !occursin(r"[pP]endulum.*Plot\(\)", rmFuncj)

            src *= """
              $rmFuncj
            """
          end
        end
        src *= """
        }
        """
      elseif occursin("removeAnimations", rmFunc)
        src *= """
        function $rmFunc {
        """
        for j in eachindex(rmFuncs)
          rmFuncj = rmFuncs[j]
          if (i != j)
            src *= """
              $rmFuncj
            """
          end
        end
        src *= """
        }
        """
      end
    end
  end
  write(filepath1, src)
  write(filepath2, src)
end

"""
    lx_includemd(com, _)

Include a Markdown fragment file into the current page during parsing,
ensuring all headings are properly registered in the Table of Contents (\\toc).
The path is relative to the site root folder.

Usage in a `.md` file:
    \\includemd{LFS/_motivations.md}
"""
function lx_includemd(com, _)
    rpath = strip(Franklin.content(com.braces[1]))
    fpath = isabspath(rpath) ? rpath : joinpath(Franklin.path(:folder), rpath)
    cur_p = get(Franklin.FD_ENV, :CUR_PATH, "")
    if !isempty(cur_p)
        Franklin.register_md_dependency(cur_p, rpath)
    end
    content = read(fpath, String)
    return replace(content, r"\\juliacode\{([^}]+)\}" => m -> begin
        source_path = strip(match(r"\\juliacode\{([^}]+)\}", m).captures[1])
        source_file = isabspath(source_path) ? source_path : joinpath(Franklin.path(:folder), source_path)
        read(source_file, String)
    end)
end

"""
    lx_juliacode(com, _)

Insert a Julia source file into the surrounding fenced code block. The path is
relative to the site root, like `\\includemd`.
"""
function lx_juliacode(com, _)
    rpath = strip(Franklin.content(com.braces[1]))
    fpath = isabspath(rpath) ? rpath : joinpath(Franklin.path(:folder), rpath)
    cur_p = get(Franklin.FD_ENV, :CUR_PATH, "")
    if !isempty(cur_p)
        Franklin.register_md_dependency(cur_p, rpath)
    end
    return read(fpath, String)
end

function lx_codelisting(com, lxdefs)
    label = strip(Franklin.content(com.braces[1]))
    caption = Franklin.content(com.braces[2])
    eqref_pattern = r"\\eqref\{([^{}]+)\}"
    caption = replace(caption, eqref_pattern => matched -> begin
        captures = match(eqref_pattern, matched).captures
        ref = strip(captures[1])
        number = get(Franklin.PAGE_EQREFS, Franklin.refstring(ref), nothing)
        isnothing(number) && return String(matched)
        "<a class=\"eqref\" href=\"#$(html_escape(ref))\">($(number))</a>"
    end)
    rpath = strip(Franklin.content(com.braces[3]))
    source_path = strip(Franklin.content(com.braces[4]))
    source_file = isabspath(source_path) ? source_path : joinpath(Franklin.path(:folder), source_path)

    cur_p = get(Franklin.FD_ENV, :CUR_PATH, "")
    if !isempty(cur_p)
        Franklin.register_md_dependency(cur_p, source_path)
    end

    source = read(source_file, String)
    fenced = "```julia:$(rpath)\n$(source)\n```"
    code_html = Franklin.resolve_code_block(SubString(fenced))
    marker = "<span class=\"code-listing-marker\" data-label=\"$(label)\">$(caption)</span>"
    return marker * "\n" * code_html
end

filter!(p -> p.first != raw"\codelisting", Franklin.LX_INTERNAL_COMMANDS)
push!(Franklin.LX_INTERNAL_COMMANDS,
      raw"\codelisting" => Franklin.LxDef(raw"\codelisting", 4, Franklin.subs("")))
Franklin.GLOBAL_LXDEFS[raw"\codelisting"] = Franklin.LxDef(raw"\codelisting", 4, Franklin.subs(""))
@eval Franklin function lx_codelisting(com, lxdefs)
    helper = getproperty(Franklin.utils_module(), :lx_codelisting)
    return Base.invokelatest(helper, com, lxdefs)
end

"""
    hfun_include_md(params)

Include a Markdown fragment file into the current page.
The path is relative to the site root folder.
Files whose names start with `_` are not rendered as standalone pages by
Franklin, making them safe to use as pure include fragments.

Usage in a `.md` file:
    {{include_md LFS/_motivations.md}}
"""
function hfun_include_md(params)
    rpath   = strip(params[1])
    fpath   = joinpath(Franklin.path(:folder), rpath)
    cur_p   = get(Franklin.FD_ENV, :CUR_PATH, "")
    if !isempty(cur_p)
        Franklin.register_md_dependency(cur_p, rpath)
    end
    content = read(fpath, String)
    return fd2html(content, internal=true)
end

function hfun_render_js()
  title = locvar("title")
  type = locvar("type")
  vars = locvar("vars")
  render_js(title, type, vars)
  if (type != "attractor")
    ids, funcs, labels = getButtonVars();
    render_rmPlot(funcs, ids, title)
  end
  title = replace(title, " " => "_")
  if (type == "attractor")
    HTML = """<script src="/libs/rendered/attractor_$title.js"></script><script src="/libs/rendered/rmPlot_$title.js"></script>"""
  else
    HTML = """<script src="/libs/rendered/RKF45_$title.js"></script><script src="/libs/rendered/rmPlot_$title.js"></script>"""
  end
  return HTML
end

"""
Override built-in hfun_toc to resolve equation references and bibliography references in TOC titles.
"""
function hfun_toc(params::Vector{String})::String
    if length(params) != 2
        throw(Franklin.HTMLFunctionError("I found a {{toc ...}} block and expected 2 parameters but got $(length(params)). Verify."))
    end
    isempty(Franklin.PAGE_HEADERS) && return ""

    min = 0
    max = 100
    try
        min = parse(Int, params[1])
        max = parse(Int, params[2])
    catch
        throw(Franklin.HTMLFunctionError("I found a {{toc min max}} but couldn't parse min/max. Verify."))
    end

    headers = filter(p -> min ≤ p.second[3] ≤ max, Franklin.PAGE_HEADERS)
    isempty(headers) && return ""

    levels = [h[3] for (rs, h) in headers]
    unique_levels = sort(unique(levels))
    level_map = Dict(lvl => i for (i, lvl) in enumerate(unique_levels))
    new_levels = [level_map[lvl] for lvl in levels]
    open_lists = 0
    curskip = 0
    inner = ""

    for (li, (rs, h)) in enumerate(headers)
        lvl = h[3]
        skip = new_levels[li]
        if skip ≤ curskip
            inner *= "</li>"
            for _ in 1:(curskip - skip)
                inner *= "</ol></li>"
                open_lists -= 1
            end
            inner *= "<li>"
        else
            for _ in 1:(skip - curskip)
                inner *= "<ol><li>"
                open_lists += 1
            end
        end

        title = String(h[1])
        # Replace {{href EQR ...}} with the equation number from PAGE_EQREFS
        title = replace(title, r"\{\{\s*href\s+EQR\s+([^\}]+?)\s*\}\}" => (m -> begin
            mat = match(r"\{\{\s*href\s+EQR\s+([^\}]+?)\s*\}\}", m)
            if mat !== nothing
                ref = strip(mat.captures[1])
                string(get(Franklin.PAGE_EQREFS, ref, m))
            else
                m
            end
        end))
        # Replace {{href BIBR ...}} with the citation from PAGE_BIBREFS
        title = replace(title, r"\{\{\s*href\s+BIBR\s+([^\}]+?)\s*\}\}" => (m -> begin
            mat = match(r"\{\{\s*href\s+BIBR\s+([^\}]+?)\s*\}\}", m)
            if mat !== nothing
                ref = strip(mat.captures[1])
                string(get(Franklin.PAGE_BIBREFS, ref, m))
            else
                m
            end
        end))
        # Strip <span class="eqref"> or <span class="bibref"> tags in the TOC link text
        title = replace(title, r"<span class=\"(eqref|bibref)\">([\s\S]*?)</span>" => s"\2")

        inner *= Franklin.html_ahref_key(rs, title)
        curlvl = lvl
        curskip = skip
    end

    for _ in 1:open_lists
        inner *= "</li></ol>"
    end
    return "<div class=\"franklin-toc\">" * inner * "</div>"
end

# ==============================================================================
# Multi-line numbered \begin{align} support
# ==============================================================================

# Use a zero-argument environment definition; Franklin's utility-environment
# parser otherwise incorrectly requires an argument after `\begin{align}`.
# `\def` is a zero-argument notation macro here, not TeX's definition primitive.
def_macro = Franklin.LxDef(raw"\def", 0, SubString(raw"\stackrel{\text{def}}{=}"))
filter!(p -> p.first != raw"\def", Franklin.LX_INTERNAL_COMMANDS)
push!(Franklin.LX_INTERNAL_COMMANDS, raw"\def" => def_macro)
Franklin.GLOBAL_LXDEFS[raw"\def"] = def_macro

align_def = Franklin.LxDef("align", 0, "" => "")
filter!(p -> p.first != "align", Franklin.LX_INTERNAL_ENVIRONMENTS)
push!(Franklin.LX_INTERNAL_ENVIRONMENTS, "align" => align_def)
Franklin.GLOBAL_LXDEFS["align"] = align_def
align_star_def = Franklin.LxDef("align*", 0, "" => "")
filter!(p -> p.first != "align*", Franklin.LX_INTERNAL_ENVIRONMENTS)
push!(Franklin.LX_INTERNAL_ENVIRONMENTS, "align*" => align_star_def)
Franklin.GLOBAL_LXDEFS["align*"] = align_star_def

function normalize_align_pair_markers(s::AbstractString)
    line = String(s)
    runs = Tuple{Int,Int,Int}[] # marker offset, marker length, comma split point
    for m in eachmatch(r"&{2,}(?==)", line)
        prefix_end = prevind(line, m.offset)
        prefix = prefix_end < firstindex(line) ? "" : SubString(line, firstindex(line), prefix_end)
        comma = findlast(',', prefix)
        isnothing(comma) && continue
        push!(runs, (m.offset, length(m.match), nextind(line, comma)))
    end
    isempty(runs) && return line

    starts = Dict(run[1] => run for run in runs)
    insertions = Set(run[3] for run in runs)
    out = IOBuffer()
    i = firstindex(line)
    while i <= lastindex(line)
        i in insertions && write(out, '&')
        if haskey(starts, i)
            _, marker_length, _ = starts[i]
            write(out, '&')
            i = nextind(line, i, marker_length)
        else
            write(out, line[i])
            i = nextind(line, i)
        end
    end
    return String(take!(out))
end

function render_align_environment(com, numbered::Bool)
    rows = Franklin.split_align_lines(Franklin.content(com))
    rows = normalize_align_pair_markers.(rows)
    html = IOBuffer()
    anchors = IOBuffer()
    tagged_rows = String[]
    alignment_pairs = maximum((cld(count(==('&'), row) + 1, 2) for row in rows); init=1)
    has_secondary_alignment = alignment_pairs > 1

    for row in rows
        line = strip(row)
        isempty(line) && continue

        unnumbered = !numbered || occursin(r"\\(nonumber|notag)", line)
        line = replace(line, r"\\(nonumber|notag)(?:\{\})?" => "")
        line = replace(line, r"\\def(?![A-Za-z])" => raw"\stackrel{\text{def}}{=}")

        label_match = match(r"\\label\{([^}]+)\}", line)
        number = nothing
        if !unnumbered
            Franklin.PAGE_EQREFS[Franklin.PAGE_EQREFS_COUNTER] += 1
            number = Franklin.PAGE_EQREFS[Franklin.PAGE_EQREFS_COUNTER]
        end
        if !isnothing(label_match)
            label = Franklin.refstring(label_match.captures[1])
            if !isnothing(number)
                Franklin.PAGE_EQREFS[label] = number
            end
            write(anchors, "<a id=\"$label\" class=\"anchor\"></a>")
            line = replace(line, label_match.match => "")
        end

        tag = isnothing(number) ? (numbered ? "\\notag" : "") : "\\tag{$number}"
        push!(tagged_rows, strip(line) * " " * tag)
    end

    write(html, String(take!(anchors)))
    env = has_secondary_alignment ? (numbered ? "alignat" : "alignedat") :
          (numbered ? "align" : "aligned")
    class = numbered ? "nonumber numbered-align" : "nonumber"
    display = if has_secondary_alignment
        "\\[\\begin{$env}{$alignment_pairs}\n" * join(tagged_rows, " " * "\\\\" * "\n") * "\n\\end{$env}\\]"
    else
        "\\[\\begin{$env}\n" * join(tagged_rows, " " * "\\\\" * "\n") * "\n\\end{$env}\\]"
    end
    write(html, "<div class=\"$class\">" * display * "</div>")
    return String(take!(html))
end

env_align(com, _) = render_align_environment(com, true)
env_align_star(com, _) = render_align_environment(com, false)

Core.eval(Franklin, quote
    function env_align(com, lxdefs)
        return getfield(utils_module(), :env_align)(com, lxdefs)
    end
    function var"env_align*"(com, lxdefs)
        return getfield(utils_module(), :env_align_star)(com, lxdefs)
    end
end)

"""
Split LaTeX content on `\\\\` outside of nested environments (like pmatrix, bmatrix, etc.)
and curly braces.
"""
# Override Franklin's convert_math_block directly inside module Franklin
Core.eval(Franklin, quote
    function split_align_lines(s::AbstractString)
        lines = String[]
        buf = IOBuffer()
        i = firstindex(s)
        n = lastindex(s)
        brace_depth = 0
        env_depth = 0

        while i <= n
            if s[i] == '\\'
                rem = SubString(s, i)
                if startswith(rem, "\\\\") && brace_depth == 0 && env_depth == 0
                    push!(lines, String(take!(buf)))
                    i = nextind(s, i, 2)
                    continue
                elseif startswith(rem, r"^\\begin\{")
                    m = match(r"^\\begin\{[^}]+\}", rem)
                    env_depth += 1
                    write(buf, m.match)
                    i = nextind(s, i, ncodeunits(m.match))
                    continue
                elseif startswith(rem, r"^\\end\{")
                    m = match(r"^\\end\{[^}]+\}", rem)
                    env_depth = max(0, env_depth - 1)
                    write(buf, m.match)
                    i = nextind(s, i, ncodeunits(m.match))
                    continue
                else
                    if i < n && (s[nextind(s, i)] == '{' || s[nextind(s, i)] == '}')
                        write(buf, s[i])
                        write(buf, s[nextind(s, i)])
                        i = nextind(s, i, 2)
                        continue
                    end
                end
            elseif s[i] == '{'
                brace_depth += 1
            elseif s[i] == '}'
                brace_depth = max(0, brace_depth - 1)
            end

            write(buf, s[i])
            i = nextind(s, i)
        end

        rem = String(take!(buf))
        if !isempty(strip(rem))
            push!(lines, rem)
        end
        return lines
    end

    function convert_math_block_with_alignment(β::OCBlock, lxdefs::Vector{LxDef})::String
        pm = get(MATH_BLOCKS_PARENS, β.name) do
            throw(MathBlockError("Unrecognised math block name."))
        end

        inner = chop(β.ss, head=pm[1], tail=pm[2])
        htmls = IOBuffer()

        # Franklin may wrap the display in an internal MATH_I block. Extract
        # the unique environment from the original block instead of relying on
        # a regex over Franklin's chopped or wrapped math string.
        source = String(β.ss)
        align_begin = raw"\begin{_franklin_align_}"
        align_end = raw"\end{_franklin_align_}"
        begin_range = findfirst(align_begin, source)
        end_range = isnothing(begin_range) ? nothing :
            findnext(align_end, source, nextind(source, last(begin_range)))
        align_body = if isnothing(begin_range) || isnothing(end_range)
            nothing
        else
            body_start = nextind(source, last(begin_range))
            body_end = prevind(source, first(end_range))
            body_start > body_end ? "" : source[body_start:body_end]
        end

        if β.name != :MATH_A && !isnothing(align_body)
            lines = getfield(utils_module(), :normalize_align_pair_markers).(split_align_lines(align_body))
            new_lines = String[]
            last_num = PAGE_EQREFS[PAGE_EQREFS_COUNTER]
            alignment_pairs = maximum((cld(count(==('&'), line) + 1, 2) for line in lines); init=1)
            has_secondary_alignment = alignment_pairs > 1

            for line in lines
                trimmed = strip(line)
                if isempty(trimmed)
                    continue
                end

                is_unnumbered = occursin(r"\\(nonumber|notag)", line)
                line_clean = replace(line, r"\\(nonumber|notag)(?:\{\})?" => "")

                if is_unnumbered
                    push!(new_lines, strip(line_clean) * " \\notag")
                else
                    PAGE_EQREFS[PAGE_EQREFS_COUNTER] += 1
                    curnum = PAGE_EQREFS[PAGE_EQREFS_COUNTER]
                    last_num = curnum

                    m_lbl = match(MATH_LABEL_PAT, line_clean)
                    if !isnothing(m_lbl)
                        name = refstring(m_lbl.captures[1])
                        write(htmls, "<a id=\"$name\" class=\"anchor\"></a>")
                        line_clean = replace(line_clean, MATH_LABEL_PAT => "")
                        PAGE_EQREFS[name] = curnum
                    end

                    push!(new_lines, strip(line_clean) * " \\tag{$curnum}")
                end
            end

            env = has_secondary_alignment ? "alignat" : "align"
            begin_env = has_secondary_alignment ? "\\begin{$env}{$alignment_pairs}" : "\\begin{$env}"
            write(htmls, "<div class=\"nonumber\" style=\"counter-reset: eqnum $last_num;\">\\[$begin_env\n" * join(new_lines, " \\\\\n") * "\n\\end{$env}\\]</div>")
            return String(take!(htmls))
        end

        if β.name ∉ [:MATH_A, :MATH_I]
            PAGE_EQREFS[PAGE_EQREFS_COUNTER] += 1
            matched = match(MATH_LABEL_PAT, inner)
            if !isnothing(matched)
                name = refstring(matched.captures[1])
                write(htmls, "<a id=\"$name\" class=\"anchor\"></a>")
                inner = replace(inner, MATH_LABEL_PAT => "")
                PAGE_EQREFS[name] = PAGE_EQREFS[PAGE_EQREFS_COUNTER]
            end
        end

        inner = replace(inner, r"\\begin\{_franklin_align_\}" => raw"\begin{align}")
        inner = replace(inner, r"\\end\{_franklin_align_\}" => raw"\end{align}")

        write(htmls, pm[3], convert_md_math(inner, lxdefs, from(β)), pm[4])
        return String(take!(htmls))
    end

    # Route math blocks through the site converter explicitly. This specialized
    # method prevents Franklin's precompiled generic block converter from
    # bypassing the per-line align handling above.
    function convert_juliacode_block(β::OCBlock, lxdefs::Vector{LxDef})
        lang, rpath, code = parse_fenced_block(β.ss)
        include_match = match(r"\\juliacode\{([^}]+)\}", code)
        if lang != "julia" || isnothing(rpath) || isnothing(include_match)
            return convert_block_default(β, lxdefs)
        end

        source_path = strip(include_match.captures[1])
        full_path = isabspath(source_path) ? source_path : joinpath(Franklin.path(:folder), source_path)
        code = read(full_path, String)

        fenced = "```$(lang):$(rpath)\n$(code)\n```"
        return resolve_code_block(SubString(fenced))
    end

    function convert_block_default(β::AbstractBlock, lxdefs::Vector{LxDef})::String
        β isa HTML_SPCH && return ifelse(isempty(β.r), β.ss, β.r)
        βn = β.name
        βn ∈ MD_HEADER && return convert_header(β, lxdefs)
        βn == :CODE_INLINE && return html_code_inline(stent(β) |> htmlesc)
        βn == :CODE_BLOCK_LANG && return resolve_code_block(β.ss)
        βn == :CODE_BLOCK! && return resolve_code_block(β.ss; shortcut=true)
        βn == :CODE_REPL && return resolve_code_block(β.ss; repl=true)
        βn == :CODE_PKG && return resolve_code_block(β.ss; pkg=true)
        βn == :CODE_HELP && return resolve_code_block(β.ss; help=true)
        βn == :CODE_SHELL && return resolve_code_block(β.ss; shell=true)
        βn == :CODE_BLOCK && return html_code(stent(β), "{{fill lang}}")
        βn == :CODE_BLOCK_IND && return convert_indented_code_block(β.ss)
        βn == :ESCAPE && return chop(β.ss, head=3, tail=3)
        βn == :FOOTNOTE_REF && return convert_footnote_ref(β)
        βn == :FOOTNOTE_DEF && return convert_footnote_def(β, lxdefs)
        βn == :LINK_DEF && return ""
        βn == :DOUBLE_BRACE && return β.ss
        βn == :HORIZONTAL_RULE && return "<hr />"
        βn ∈ MATH_BLOCKS_NAMES && return convert_math_block(β, lxdefs)
        if βn == :DIV
            raw_cont = stent(β)
            cont = convert_md(raw_cont, lxdefs; isrecursive=true, has_mddefs=false, nostripp=true) |> simplify_ps
            divname = chop(otok(β).ss, head=2, tail=0)
            return html_div(replace(divname, "," => " "), cont)
        end
        return ""
    end

    function convert_block(β::AbstractBlock, lxdefs::Vector{LxDef})::AS
        if β isa OCBlock && β.name ∈ MATH_BLOCKS_NAMES
            return convert_math_block_with_alignment(β, lxdefs)
        elseif β isa OCBlock && startswith(β.ss, "```") && occursin(r"\\juliacode\{[^}]+\}", β.ss)
            lang, rpath, code = parse_fenced_block(β.ss)
            if lang == "julia" && !isnothing(rpath) && occursin(r"\\juliacode\{[^}]+\}", code)
                return convert_juliacode_block(β, lxdefs)
            end
        end
        return convert_block_default(β, lxdefs)
    end

    # `convert_inter_html` is precompiled by Franklin and can retain its
    # original `convert_block` call after the site's converter is installed.
    # Override the insertion loop and use latest-world dispatch so included
    # Julia source blocks reach the site converter before Franklin evaluates
    # their literal `\\juliacode{...}` contents.
    function convert_inter_html(ihtml::AS,
                                blocks::Vector{<:AbstractBlock},
                                lxdefs::Vector{LxDef})::String
        allmatches = collect(eachmatch(INSERT_PAT, ihtml))
        isempty(allmatches) && return ihtml

        strlen = lastindex(ihtml)
        htmls = IOBuffer()
        head = 1
        for (i, m) in enumerate(allmatches)
            leftp = !isnothing(m.captures[1])
            lefts = length(m.captures[2]) > 1 ? " " : ""
            rightp = !isnothing(m.captures[3])
            prev = prevind(ihtml, m.offset)
            write(htmls, subs(ihtml, head:prev))
            leftp && !rightp && write(htmls, "<p>")
            write(htmls, lefts)
            write(htmls, Base.invokelatest(convert_block, blocks[i], lxdefs))
            rightp && !leftp && write(htmls, "</p>")
            head = nextind(ihtml, m.offset, length(m.match))
        end
        head <= strlen && write(htmls, subs(ihtml, head:strlen))
        return String(take!(htmls))
    end

    # ==============================================================================
    # Included Markdown dependency tracking and automatic rebuilding
    # ==============================================================================

    if !isdefined(Franklin, :MD_INCLUDED_IN)
        const MD_INCLUDED_IN = Dict{String, Set{String}}()
    end

    function norm_site_rpath(p::AbstractString)
        s = replace(strip(p), '\\' => '/')
        s = lstrip(s, '/')
        if startswith(s, "./")
            s = s[3:end]
        end
        return normpath(s)
    end

    function register_md_dependency(parent_rpath::AbstractString, child_rpath::AbstractString)
        parent_norm = norm_site_rpath(parent_rpath)
        child_norm  = norm_site_rpath(child_rpath)
        if isempty(parent_norm) || isempty(child_norm) || parent_norm == child_norm
            return
        end
        s = get!(MD_INCLUDED_IN, child_norm) do
            Set{String}()
        end
        push!(s, parent_norm)
    end

    function scan_md_includes!()
        root = try path(:folder) catch; pwd() end
        site_dir = try path(:site) catch; joinpath(root, "__site") end
        for (dir, _, files) in walkdir(root)
            if (!isempty(site_dir) && startswith(dir, site_dir)) || occursin(".git", dir) || occursin("node_modules", dir)
                continue
            end
            for file in files
                if endswith(file, ".md")
                    fpath = joinpath(dir, file)
                    rel_parent = norm_site_rpath(relpath(fpath, root))
                    content = try
                        read(fpath, String)
                    catch
                        ""
                    end
                    for m in eachmatch(Regex("\\\\includemd\\{([^}]+)\\}"), content)
                        register_md_dependency(rel_parent, strip(m.captures[1]))
                    end
                    for m in eachmatch(Regex("\\\\juliacode\\{([^}]+)\\}"), content)
                        register_md_dependency(rel_parent, strip(m.captures[1]))
                    end
                    for m in eachmatch(Regex("\\{\\{include_md\\s+([^}]+)\\}\\}"), content)
                        register_md_dependency(rel_parent, strip(m.captures[1]))
                    end
                end
            end
        end
    end

    scan_md_includes!()

    function get_parent_pages(child_rpath::AbstractString)
        norm_child = norm_site_rpath(child_rpath)
        parents = Set{String}()
        queue = [norm_child]
        visited = Set{String}([norm_child])

        while !isempty(queue)
            curr = popfirst!(queue)

            matching_parents = Set{String}()
            if haskey(MD_INCLUDED_IN, curr)
                union!(matching_parents, MD_INCLUDED_IN[curr])
            end
            for (k, ps) in MD_INCLUDED_IN
                if k == curr || normpath(k) == normpath(curr) || endswith(curr, "/" * k) || endswith(k, "/" * curr)
                    union!(matching_parents, ps)
                end
            end

            for p in matching_parents
                if p ∉ visited
                    push!(visited, p)
                    push!(queue, p)
                    _, pf = splitdir(p)
                    if !startswith(pf, "_")
                        push!(parents, p)
                    end
                end
            end
        end
        return parents
    end

    # Override process_file_err to suppress standalone rendering of _*.md fragments.
    # Without this, Franklin writes an HTML file for every fragment on disk, which
    # causes LiveServer to broadcast a reload to ALL connected browsers -- even those
    # viewing the parent page that includes the fragment.  Skipping the write means
    # only the parent page rebuild (done in the fd_loop wrapper below) triggers a reload.
    function process_file_err(case::Symbol, fpair::Pair{String, String},
                              head::AS="", pgfoot::AS="", foot::AS="",
                              t::Float64=0.0)::Nothing
        inp = joinpath(fpair...)

        layout_p = try path(:layout)   catch; "" end
        lit_p    = try path(:literate) catch; "" end
        rss_p    = try path(:rss)      catch; "" end

        if (!isempty(layout_p) && startswith(inp, layout_p)) ||
           (!isempty(lit_p)    && startswith(inp, lit_p))    ||
           (!isempty(rss_p)    && startswith(inp, rss_p))    ||
           endswith(inp, "config.md") ||
           endswith(inp, "utils.jl")  ||
           (case == :md && startswith(fpair.second, "_"))
            # Skip: infra / literate / rss / fragment files don't need a standalone page.
            @goto end_copyblock
        end

        outp = form_output_path(fpair.first, fpair.second, case)
        if case == :md
            FD_ENV[:SOURCE] = get_rpath(inp)
            convert_and_write(fpair..., head, pgfoot, foot, outp)
        elseif case == :html
            FD_ENV[:SOURCE] = get_rpath(inp)
            set_cur_rpath(joinpath(fpair...))
            set_page_env()
            raw_html  = read(inp, String)
            cond_add  = globvar(:generate_sitemap)::Bool && FD_ENV[:FULL_PASS]::Bool
            cond_add  && add_sitemap_item(html=true)
            proc_html = convert_html(raw_html) |> postprocess_page
            write(outp, proc_html)
        else
            if !isfile(outp) || (mtime(outp) < t && !filecmp(inp, outp))
                cp(inp, outp, force=true)
            end
            @label end_copyblock
        end
        FD_ENV[:SOURCE] = ""
        FD_ENV[:FULL_PASS] || FD_ENV[:SILENT_MODE] || rprint("\u2192 page updated [\u2713]")
        return nothing
    end

     # Wrap fd_loop once to rebuild parent pages when an included fragment changes.
    # Only wrap once; guard against re-wrapping on each process_utils reload.
    if !isdefined(Franklin, :_fd_loop_fragment_wrapped)
        const _fd_loop_fragment_wrapped = Ref(false)
    end
    if !isdefined(Franklin, :_fd_loop_hashes_initialized)
        const _fd_loop_hashes_initialized = Ref(false)
    end
    if !isdefined(Franklin, :_fd_loop_fragment_mtimes)
        const _fd_loop_fragment_mtimes = Dict{String,Float64}()
    end
    if !isdefined(Franklin, :_fd_loop_fragment_mtimes_initialized)
        const _fd_loop_fragment_mtimes_initialized = Ref(false)
    end
    if !isdefined(Franklin, :_fd_loop_native)
        const _fd_loop_native = Ref{Any}(nothing)
    end

    # This helper is redefined whenever Franklin reloads utils.jl. The stable
    # fd_loop wrapper below calls it through invokelatest, so watcher fixes can
    # be picked up on a full pass without stacking wrappers.
    function _fd_loop_with_md_fragments(cycle_counter::Int,
                                        fw::LiveServer.FileWatcher,
                                        watched_files::NamedTuple)::Nothing
        if !_fd_loop_fragment_mtimes_initialized[]
            _fd_loop_native[](cycle_counter, fw, watched_files)
            for fragment in keys(MD_INCLUDED_IN)
                fpath = joinpath(path(:folder), fragment)
                isfile(fpath) || continue
                _fd_loop_fragment_mtimes[fpath] = mtime(fpath)
            end
            _fd_loop_hashes_initialized[] = true
            _fd_loop_fragment_mtimes_initialized[] = true
            return nothing
        end

        # Refresh the include map before polling. This also covers fragments
        # discovered after Franklin's initial full pass or added while serving.
        scan_md_includes!()
        changed_frags = String[]
        for fragment in keys(MD_INCLUDED_IN)
            fpath = joinpath(path(:folder), fragment)
            isfile(fpath) || continue
            current_mtime = mtime(fpath)
            previous_mtime = get(_fd_loop_fragment_mtimes, fpath, current_mtime)
            if current_mtime > previous_mtime
                _fd_loop_fragment_mtimes[fpath] = current_mtime
                rp = try norm_site_rpath(relpath(fpath, path(:folder))) catch; "" end
                isempty(rp) || push!(changed_frags, rp)
            end
        end

        # Let Franklin handle regular pages and file touches as usual.
        _fd_loop_native[](cycle_counter, fw, watched_files)

        isempty(changed_frags) && return nothing

        layout = path(:layout)
        head   = read(joinpath(layout, "head.html"),      String)
        pgfoot = read(joinpath(layout, "page_foot.html"), String)
        foot   = read(joinpath(layout, "foot.html"),      String)

        parents_to_rebuild = Set{String}()
        scan_md_includes!()
        for frag_rp in changed_frags
            union!(parents_to_rebuild, get_parent_pages(frag_rp))
        end
        for parent_rp in sort!(collect(parents_to_rebuild))
            parent_full = joinpath(path(:folder), parent_rp)
            isfile(parent_full) || continue
            parent_dir, parent_file = splitdir(parent_full)
            parent_fpair = parent_dir => parent_file
            rprint("\u2192 rebuilding $(parent_rp) (included fragment changed)")
            process_file(:md, parent_fpair, head, pgfoot, foot, 0.0)
            if haskey(watched_files.md, parent_fpair)
                watched_files.md[parent_fpair] = mtime(parent_full)
            end
        end
        return nothing
    end

    if !_fd_loop_fragment_wrapped[]
        _fd_loop_fragment_wrapped[] = true
        let _orig_fd_loop = fd_loop
            _fd_loop_native[] = _orig_fd_loop
            function fd_loop(cycle_counter::Int,
                             fw::LiveServer.FileWatcher,
                             watched_files::NamedTuple)::Nothing
                Base.invokelatest(_fd_loop_with_md_fragments,
                                  cycle_counter, fw, watched_files)
            end
        end
    end
end)
