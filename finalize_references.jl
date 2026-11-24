using Base64

const PAGE_FOOT_TAG = r"<div\b(?=[^>]*\bclass\s*=\s*(?:\"[^\"]*\bpage-foot\b[^\"]*\"|'[^']*\bpage-foot\b[^']*'|[^\s>]*\bpage-foot\b[^\s>]*))[^>]*>"

reference_html_escape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;", "'" => "&#39;")

function html_attribute(attributes, name)
  pattern = Regex("\\b" * name * "\\s*=\\s*(?:\"([^\"]*)\"|'([^']*)'|([^\\s>]+))")
  matched = match(pattern, attributes)
  isnothing(matched) && return nothing
  return something(matched.captures...)
end

function insert_reference_list(html, section)
  if occursin(PAGE_FOOT_TAG, html)
    return replace(html, PAGE_FOOT_TAG => matched -> section * matched; count=1)
  elseif occursin("<!-- CONTENT ENDS HERE -->", html)
    return replace(html, "<!-- CONTENT ENDS HERE -->" => section * "<!-- CONTENT ENDS HERE -->"; count=1)
  elseif occursin(r"(?i)</article\s*>", html)
    return replace(html, r"(?i)</article\s*>" => section * "</article>"; count=1)
  elseif occursin(r"(?i)</body\s*>", html)
    return replace(html, r"(?i)</body\s*>" => section * "</body>"; count=1)
  end
  return html * section
end

function bibtex_source_entry(bib_path, key)
  source = read(bib_path, String)
  pattern = Regex("(?is)@([a-z]+)\\s*([({])\\s*" * replace(key, r"([\\.^$|?*+()\[\]{}])" => s"\\\1") * "\\s*,")
  matched = match(pattern, source)
  isnothing(matched) && error("Citation key '$key' was not found in $bib_path")
  start = matched.offset
  openchar = only(matched.captures[2])
  closechar = openchar == '{' ? '}' : ')'
  depth, quoted, escaped = 0, false, false
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
      depth == 0 && return strip(source[start:pos])
    end
  end
  error("Unclosed BibTeX entry '$key' in $bib_path")
end

function append_fallback_reference_list(html, html_path, site_dir)
  source_md = joinpath(dirname(site_dir), replace(relpath(html_path, site_dir), r"\.html$" => ".md"))
  isfile(source_md) || return nothing
  source = read(source_md, String)
  ref_match = match(r"(?m)^\s*@def\s+ref\s*=\s*[\"']([^\"']+)[\"']", source)
  isnothing(ref_match) && return nothing
  bib_path = normpath(joinpath(dirname(source_md), ref_match.captures[1]))
  isfile(bib_path) || return nothing

  citations = Dict{String,String}()
  for m in eachmatch(r"<a\b([^>]*)>", html)
    attrs = m.captures[1]
    occursin(r"(?:^|\s)class=(?:\"[^\"]*\bcitation\b[^\"]*\"|'[^']*\bcitation\b[^']*')", attrs) || continue
    href = something(html_attribute(attrs, "href"), "")
    target_match = match(r"^#(reference-([[:xdigit:]]+))$", href)
    isnothing(target_match) && continue
    target, key_hex = target_match.captures
    key = String(hex2bytes(key_hex))
    citations[target] = bibtex_source_entry(bib_path, key)
  end
  isempty(citations) && return nothing

  items = join("<li id=\"$(target)\" style=\"padding-left:2em;text-indent:-2em;\"><pre><code>" * reference_html_escape(entry) * "</code></pre></li>" for (target, entry) in sort(collect(citations); by=first))
  section = "<section class=\"reference-list\" aria-labelledby=\"reference-list-title\"><h1 id=\"reference-list-title\">Reference list</h1><ul style=\"list-style:none;margin:0;padding-left:0;\">$(items)</ul></section>"
  return insert_reference_list(html, section)
end

function append_link_label_reference_list(html)
  cited = Dict{String,String}()
  for m in eachmatch(r"<a\b([^>]*)>([\s\S]*?)</a>", html)
    attributes, label = m.captures
    occursin(r"(?:^|\s)class=(?:\"[^\"]*\bcitation\b[^\"]*\"|'[^']*\bcitation\b[^']*')", attributes) || continue
    href = something(html_attribute(attributes, "href"), "")
    startswith(href, "#reference-") || continue
    cited[href[2:end]] = label
  end
  isempty(cited) && return nothing

  items = join("<li id=\"$(target)\" style=\"padding-left:2em;text-indent:-2em;\">" * label * "</li>" for (target, label) in sort(collect(cited); by=first))
  section = "<section class=\"reference-list\" aria-labelledby=\"reference-list-title\"><h1 id=\"reference-list-title\">Reference list</h1><ul style=\"list-style:none;margin:0;padding-left:0;\">$(items)</ul></section>"
  return insert_reference_list(html, section)
end

function append_reference_lists(site_dir, output_names)
  citation_span = r"<span\b([^>]*)>([\s\S]*?)</span>"
  citation_anchor = r"<a\b([^>]*)>([^<]*)</a>"
  preview_script = r"<script>\(function\(\)\{\s*if \(window\.__bibReferenceListInstalled[\s\S]*?</script>"
  page_locator = r"(<a class=\"citation\"[^>]*href=\"#[^\"]+\">)([^<]*)(</a>\s*\))\s*(?:&#123;|\{)([^{}&]+)(?:&#125;|\})"
  processed = 0

  for (root, _, files) in walkdir(site_dir)
    for file in files
      endswith(file, ".html") || continue
      path = joinpath(root, file)
      original = read(path, String)
      html = replace(original, preview_script => "")

      for name in output_names
        output_file = joinpath(site_dir, "assets", "fidello", "code", "output", name * ".out")
        isfile(output_file) || continue
        error_span = "<span style=\"color:red;\">// `" * raw"\output{./code/" * name * "}`: could not find the relevant output file. //</span>"
        rendered_output = "<pre><code class=\"plaintext code-output\">" * reference_html_escape(read(output_file, String)) * "</code></pre>"
        html = replace(html, error_span => rendered_output)
      end

      has_reference_list = occursin(r"<section\b[^>]*\bclass=(?:\"[^\"]*\breference-list\b[^\"]*\"|'[^']*\breference-list\b[^']*')", html)
      if has_reference_list
        html = replace(html, "<h2 id=\"reference-list-title\">Reference list</h2>" => "<h1 id=\"reference-list-title\">Reference list</h1>")
        html = replace(html, r"(<section class=\"reference-list\"[^>]*>[\s\S]*?<h1 id=\"reference-list-title\">Reference list</h1>)<ol>" => s"\1<ul style=\"list-style:none;margin:0;padding-left:0;\">")
        html = replace(html, "</ol></section>" => "</ul></section>"; count=1)
        html = replace(html, r"<li id=\"(reference-[^\"]+)\">" => matched -> begin
          id = match(r"<li id=\"(reference-[^\"]+)\">", matched).captures[1]
          "<li id=\"$(id)\" style=\"padding-left:2em;text-indent:-2em;\">"
        end)
      end

      span_matches = [m for m in eachmatch(citation_span, html) if !isnothing(html_attribute(m.captures[1], "data-bib-key")) && !isnothing(html_attribute(m.captures[1], "data-bib-entry"))]
      anchor_matches = [m for m in eachmatch(citation_anchor, html) if occursin(r"(?:^|\s)class=(?:\"[^\"]*\bcitation\b[^\"]*\"|'[^']*\bcitation\b[^']*'|[^\s>]*\bcitation\b[^\s>]*)", m.captures[1]) && !isnothing(html_attribute(m.captures[1], "data-bib-key")) && !isnothing(html_attribute(m.captures[1], "data-bib-entry"))]
      if isempty(span_matches) && isempty(anchor_matches)
        if !has_reference_list
          fallback_html = append_fallback_reference_list(html, path, site_dir)
          if !isnothing(fallback_html)
            html = fallback_html
          else
            fallback_html = append_link_label_reference_list(html)
            !isnothing(fallback_html) && (html = fallback_html)
          end
        end
        if html != original
          write(path, html)
          processed += 1
        end
        continue
      end

      has_reference_list && (html = replace(html, r"<section\b[^>]*\bclass=(?:\"[^\"]*\breference-list\b[^\"]*\"|'[^']*\breference-list\b[^']*')[^>]*>[\s\S]*?</section>" => ""; count=1))
      reference_id(key) = "reference-" * bytes2hex(codeunits(key))
      cited = Dict{String,String}()
      for m in span_matches
        key = html_attribute(m.captures[1], "data-bib-key")
        entry = html_attribute(m.captures[1], "data-bib-entry")
        cited[String(base64decode(key))] = String(base64decode(entry))
      end
      for m in anchor_matches
        key = html_attribute(m.captures[1], "data-bib-key")
        entry = html_attribute(m.captures[1], "data-bib-entry")
        cited[String(base64decode(key))] = String(base64decode(entry))
      end

      html = replace(html, citation_span => matched -> begin
        captures = match(citation_span, matched).captures
        key_encoded = html_attribute(captures[1], "data-bib-key")
        entry_encoded = html_attribute(captures[1], "data-bib-entry")
        isnothing(key_encoded) && return matched
        isnothing(entry_encoded) && return matched
        key = String(base64decode(key_encoded))
        "<a class=\"citation\" data-bib-key=\"$key_encoded\" data-bib-entry=\"$entry_encoded\" href=\"#$(reference_id(key))\">$(captures[2])</a>"
      end)
      html = replace(html, citation_anchor => matched -> begin
        captures = match(citation_anchor, matched).captures
        key_encoded = html_attribute(captures[1], "data-bib-key")
        entry_encoded = html_attribute(captures[1], "data-bib-entry")
        isnothing(key_encoded) && return matched
        isnothing(entry_encoded) && return matched
        key = String(base64decode(key_encoded))
        "<a class=\"citation\" data-bib-key=\"$key_encoded\" data-bib-entry=\"$entry_encoded\" href=\"#$(reference_id(key))\">$(captures[2])</a>"
      end)
      html = replace(html, page_locator => matched -> begin
        captures = match(page_locator, matched).captures
        captures[1] * captures[2] * ", p. " * strip(captures[4]) * captures[3]
      end)

      sorted = sort(collect(cited); by=pair -> lowercase(strip(replace(last(pair), r"<[^>]+>" => ""))))
      items = join("<li id=\"$(reference_id(first(pair)))\" style=\"padding-left:2em;text-indent:-2em;\">" * last(pair) * "</li>" for pair in sorted)
      section = "<section class=\"reference-list\" aria-labelledby=\"reference-list-title\"><h1 id=\"reference-list-title\">Reference list</h1><ul style=\"list-style:none;margin:0;padding-left:0;\">$(items)</ul></section>"
      toc_link = "<li><a href=\"#reference-list-title\">Reference list</a></li>"
      if !occursin("href=\"#reference-list-title\"", html)
        html = replace(html, "</ol></div>" => toc_link * "</ol></div>"; count=1)
      end
      html = insert_reference_list(html, section)
      write(path, html)
      processed += 1
    end
  end
  return processed
end

function validate_reference_links(html_path)
  html = read(html_path, String)
  # Citation markup appears as string literals inside the browser preview
  # script. Validate the rendered HTML, not markup-shaped text in scripts or
  # comments.
  document_html = replace(html, r"<script\b[\s\S]*?</script>" => "", r"<!--[\s\S]*?-->" => "")
  citation_markers = [m for m in eachmatch(r"<span\b([^>]*)>", document_html)
                      if !isnothing(html_attribute(m.captures[1], "data-bib-key"))]
  isempty(citation_markers) || error("Unprocessed citation marker in $html_path")
  citation_links = collect(eachmatch(r"<a\b([^>]*)>", document_html))
  is_reference_citation(attributes) =
    occursin(r"(?:^|\s)class=(?:\"[^\"]*\bcitation\b[^\"]*\"|'[^']*\bcitation\b[^']*')", attributes) &&
    startswith(something(html_attribute(attributes, "href"), ""), "#reference-")
  reference_citations = [m for m in citation_links if is_reference_citation(m.captures[1])]
  has_citations = !isempty(reference_citations)
  has_reference_list = occursin(r"<section\b[^>]*\bclass=(?:\"[^\"]*\breference-list\b[^\"]*\"|'[^']*\breference-list\b[^']*')", document_html)
  has_citations && !has_reference_list &&
    error("Reference list missing from $html_path despite citation links")
  reference_list = match(r"<section\b[^>]*\bclass=(?:\"[^\"]*\breference-list\b[^\"]*\"|'[^']*\breference-list\b[^']*')", document_html)
  page_foot = match(PAGE_FOOT_TAG, document_html)
  if !isnothing(reference_list) && !isnothing(page_foot)
    reference_list.offset < page_foot.offset ||
      error("Reference list appears after the page footer in $html_path")
  end
  ids = Set(String(m.captures[1]) for m in eachmatch(r"id=\"([^\"]+)\"", document_html))
  for m in reference_citations
    target = something(html_attribute(m.captures[1], "href"), "")[2:end]
    target in ids || error("Citation target '#$target' is missing from $html_path")
  end
  return true
end

if abspath(PROGRAM_FILE) == @__FILE__
  site_dir = joinpath(@__DIR__, "__site")
  output_dir = joinpath(site_dir, "assets", "fidello", "code", "output")
  output_names = [first(splitext(file)) for file in readdir(output_dir) if endswith(file, ".out")]
  pages = append_reference_lists(site_dir, output_names)
  fidello_page = joinpath(site_dir, "fidello", "index.html")
  validate_reference_links(fidello_page)
  println("Finalized and verified citation links; processed $pages generated page(s).")
end
