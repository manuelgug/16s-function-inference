import base64

from render_html import render

PNG = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR4nGP4z8DwHwAFAAH/iZk9HQAAAABJRU5ErkJggg==")


def test_render_embeds_images_and_rewrites_links(tmp_path):
    (tmp_path / "fig.png").write_bytes(PNG)
    md = tmp_path / "r.md"
    md.write_text("# Title\n\n![f](fig.png)\n\n> **Figure 1.** caption\n\n| a | b |\n|---|---|\n| 1 | 2 |\n\n"
                  "[other](sub/other.md#x) [web](https://example.org/a.md)\n")
    html = render(md).read_text()
    assert "<title>Title</title>" in html
    assert 'src="data:image/png;base64,' in html and "fig.png" not in html
    assert "<blockquote>" in html and "<table>" in html
    assert 'href="sub/other.html#x"' in html and 'href="https://example.org/a.md"' in html
