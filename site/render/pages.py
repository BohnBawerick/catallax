"""Shared chrome for the generated evidence pages."""
from html import escape

REPO = 'https://github.com/BohnBawerick/catallax'


def chrome(title, description, canonical, depth, current, body, source, scripts=''):
    """Wrap a generated evidence page in the shared bar and title block. `source` is (label, href, text)."""
    up = '../'*depth
    pages = [('two-kernels.html', 'Two kernels'), ('gallery.html', 'Gallery'), ('compare.html', 'Compare'), ('kernel.html', 'Kernel')]
    links = ''.join(f'\n  <a class="bar__link" href="{up}{href}"{" aria-current=\"page\"" if href == current else ""}>{name}</a>'
                    for href, name in pages)
    return f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)} | Catallax</title>
<meta name="description" content="{escape(description)}">
<link rel="canonical" href="https://kernel.zigrawbohnbawerick.org/{canonical}">
<meta name="color-scheme" content="dark">
<meta name="theme-color" content="#0b0c0d">
<link rel="icon" href="{up}brand/mark.svg" type="image/svg+xml">
<script>document.documentElement.classList.add('js')</script>
<link rel="stylesheet" href="{up}brand/fonts.css">
<link rel="stylesheet" href="{up}brand/tokens.css">
<link rel="stylesheet" href="{up}bare.css">
<link rel="stylesheet" href="{up}read.css">
<link rel="stylesheet" href="{up}evidence.css">
</head>
<body>
<svg width="0" height="0" style="position:absolute" aria-hidden="true">
<symbol id="s-projection" viewBox="0 0 48 24"><path d="M3 8.5 25 3v18L3 15.5Z"/><circle cx="37" cy="12" r="9"/><circle cx="37" cy="12" r="4"/></symbol>
</svg>
<a class="skip" href="#main">Skip to content</a>
<div class="frame" aria-hidden="true"></div>
<header class="bar">
  <a class="bar__brand" href="{up}index.html" aria-label="Catallax, home"><img src="{up}brand/mark.svg" alt="" width="28" height="28"><span>catallax</span></a>{links}
  <a class="bar__link" href="{REPO}" target="_blank" rel="noopener">Code</a>
</header>
{body}
<footer class="block mono">
  <div class="block__cell block__cell--sym"><svg class="sym sym--wide"><use href="#s-projection"/></svg></div>
  <div class="block__cell"><small>Title</small><span>{escape(title)}</span></div>
  <div class="block__cell"><small>{source[0]}</small><a href="{source[1]}" target="_blank" rel="noopener">{source[2]}</a></div>
  <div class="block__cell"><small>Drawn</small><a href="https://zigrawbohnbawerick.org/" target="_blank" rel="noopener">Zigrawbohnbawerick Pty Ltd</a></div>
  <div class="block__cell"><small>Licence</small><a href="{REPO}/blob/main/LICENSE" target="_blank" rel="noopener">MIT</a></div>
  <div class="block__cell"><small>Owner</small><a href="https://zigrawbohnbawerick.org/" target="_blank" rel="noopener">© 2026 Zigrawbohnbawerick Pty Ltd</a></div>
  <div class="block__cell"><small>Notes</small><span class="block__links"><a href="{up}why.html">Why</a><a href="{up}thanks.html">Thanks</a><a href="{up}THIRD-PARTY.txt">Third-party licences</a></span></div>
</footer>
<script src="{up}read.js"></script>{scripts}
</body>
</html>
'''


def heading(lines):
    spans = '\n'.join(f'      <span class="ln"><span>{line}</span></span>' for line in lines)
    return f'''    <h1>
{spans}
    </h1>
    <span class="dim js-only" aria-hidden="true"><i></i><b></b><i></i></span>'''
