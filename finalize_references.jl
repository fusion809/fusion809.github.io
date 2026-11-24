using Base64

reference_html_escape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;", "'" => "&#39;")

function html_attribute(attributes, name)
  pattern = Regex("\\b" * name * "\\s*=\\s*(?:\"([^\"]*)\"|'([^']*)'|([^\\s>]+))")
  matched = match(pattern, attributes)
  isnothing(matched) && return nothing
  return something(matched.captures...)
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

      has_reference_list = occursin("<section class=\"reference-list\"", html)
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
        if html != original
          write(path, html)
          processed += 1
        end
        continue
      end

      has_reference_list && (html = replace(html, r"<section class=\"reference-list\"[\s\S]*?</section>" => ""; count=1))
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

function validate_reference_links(html_path)
  html = read(html_path, String)
  occursin("<section class=\"reference-list\"", html) || error("Reference list missing from $html_path")
  any(m -> !isnothing(html_attribute(m.captures[1], "data-bib-key")), eachmatch(r"<span\b([^>]*)>", html)) && error("Unprocessed citation marker in $html_path")
  ids = Set(String(m.captures[1]) for m in eachmatch(r"id=\"([^\"]+)\"", html))
  for m in eachmatch(r"<a class=\"citation\"[^>]*href=\"#([^\"]+)\"", html)
    target = String(m.captures[1])
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
