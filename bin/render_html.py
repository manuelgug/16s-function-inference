#!/usr/bin/env python3
"""Render Markdown reports to self-contained HTML (images embedded), next to each .md file.

Usage: render_html.py report.md [more.md ...]

Each report.md becomes report.html in the same directory. Relative images are embedded as
data URIs, so the HTML file can be opened, mailed or archived on its own. Relative links to
other .md reports are rewritten to their .html twins. Blockquotes are styled as figure
explanations: reports put a "> **Figure N. ...**" block under every image.
"""
import base64
import html
import mimetypes
import re
import sys
from pathlib import Path

import markdown

CSS = """
:root { --bg:#fcfcfb; --surface:#ffffff; --ink:#1f1e1c; --ink2:#52514e; --rule:#e4e3df;
        --accent:#2a78d6; --note:#f3f2ee; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
        --bg:#1b1b1a; --surface:#232322; --ink:#ecebe7; --ink2:#b3b1ab; --rule:#3a3936;
        --accent:#6aa5ec; --note:#2b2a28; } }
:root[data-theme="dark"] { --bg:#1b1b1a; --surface:#232322; --ink:#ecebe7; --ink2:#b3b1ab; --rule:#3a3936;
        --accent:#6aa5ec; --note:#2b2a28; }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
       font:16px/1.6 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width: 1040px; margin: 0 auto; padding: 32px 16px 64px; }
h1 { font-size: 1.9rem; line-height: 1.25; margin: 0 0 1rem; }
h2 { font-size: 1.35rem; margin: 2.4rem 0 .8rem; padding-top: 1rem; border-top: 1px solid var(--rule); }
h3 { font-size: 1.1rem; margin: 1.6rem 0 .6rem; }
a { color: var(--accent); }
code { font: .9em ui-monospace, SFMono-Regular, Menlo, monospace; background: var(--note);
       padding: .1em .3em; border-radius: 4px; }
img { display:block; max-width:100%; height:auto; margin: 1.2rem auto .4rem; border-radius: 6px;
      border: 1px solid var(--rule); background:#fcfcfb; }
table { display:block; overflow-x:auto; border-collapse: collapse; margin: 1rem 0; font-size: .88rem;
        font-variant-numeric: tabular-nums; max-width: 100%; }
th, td { padding: .35rem .7rem; border-bottom: 1px solid var(--rule); text-align: left; white-space: nowrap; }
th { color: var(--ink2); font-weight: 600; border-bottom: 2px solid var(--rule); }
blockquote { margin: .4rem 0 1.6rem; padding: .8rem 1rem; background: var(--note); border-radius: 6px;
             border-left: 3px solid var(--accent); color: var(--ink); font-size: .95rem; }
blockquote p { margin: .35rem 0; }
.src { color: var(--ink2); font-size: .85rem; margin-top: 3rem; }
"""


def embed_images(text, base):
    def repl(m):
        src = m.group(2)
        if re.match(r"^(https?:|data:)", src):
            return m.group(0)
        p = base / src
        mime = mimetypes.guess_type(p.name)[0] or "application/octet-stream"
        return f'{m.group(1)}"data:{mime};base64,{base64.b64encode(p.read_bytes()).decode()}"'
    return re.sub(r'(<img[^>]*?src=)"([^"]+)"', repl, text)


def render(md_path):
    md_path = Path(md_path)
    text = md_path.read_text()
    m = re.search(r"^# (.+)$", text, re.M)
    title = m.group(1) if m else md_path.stem
    body = markdown.markdown(text, extensions=["tables", "fenced_code", "sane_lists"])
    body = embed_images(body, md_path.parent)
    body = re.sub(r'href="(?![a-z]+:)([^"#]+)\.md(#[^"]*)?"', r'href="\1.html\2"', body)
    page = (f'<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f"<title>{html.escape(title)}</title><style>{CSS}</style></head>\n<body><main>\n{body}\n"
            f'<p class="src">Rendered from <code>{html.escape(md_path.name)}</code>.</p>\n</main></body></html>\n')
    out = md_path.with_suffix(".html")
    out.write_text(page)
    return out


if __name__ == "__main__":
    for f in sys.argv[1:]:
        print(render(f))
