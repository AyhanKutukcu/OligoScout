"""web/oligoscout_data.json + web/src/* -> web/index.html (tek dosya, çevrimdışı açılır).

Akış:  python web/build_web_data.py   (motor çıktılarını paketler)
       python web/build_site.py       (siteyi yeniden üretir)
"""
import json, os
W = os.path.dirname(os.path.abspath(__file__)); S = os.path.join(W, 'src')
rd = lambda *p: open(os.path.join(*p), encoding='utf-8').read()
raw = rd(W, 'oligoscout_data.json'); dig = json.loads(raw)['engine_params'].get('registry_digest') or ''
shell = rd(S, 'shell.html').replace('__LOGO__', rd(S, 'logo.svg')).replace('__DIGEST__', dig[:16] + '…')
js = ''.join(rd(S, f) + '\n' for f in ['engine.js', 'core.js', 'viz.js', 'views1.js', 'views2.js', 'views3.js', 'init.js'])
html = ('<!doctype html>\n<html lang="tr">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
        '<title>OligoScout</title>\n<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;700&family=Space+Grotesk:wght@500;600;700&display=swap">\n'
        '<style>' + rd(S, 'styles.css') + '</style>\n</head>\n<body>\n' + shell +
        '\n<script type="application/json" id="osdata">' + raw.replace('</', '<\\/') + '</script>\n<script>\n' + js + '</script>\n</body>\n</html>\n')
open(os.path.join(W, 'index.html'), 'w', encoding='utf-8').write(html)
print('web/index.html', len(html), 'bayt')
