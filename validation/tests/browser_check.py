"""Optional visual checks; run with an environment that has Playwright installed."""
import json
from pathlib import Path
from playwright.sync_api import sync_playwright
root=Path(__file__).resolve().parents[1]
with sync_playwright() as p:
 browser=p.chromium.launch(headless=True)
 page=browser.new_page(); network=[]
 page.on('request',lambda r: network.append(r.url) if r.url.startswith(('http:','https:')) else None)
 page.context.set_offline(True)
 page.goto((root/'report/index.html').as_uri())
 findings=[]
 for width,label in [(1440,'desktop'),(390,'mobile')]:
  page.set_viewport_size({'width':width,'height':1000})
  assert page.evaluate('document.documentElement.scrollWidth <= innerWidth'),label+' page clipping'
  page.screenshot(path=str(root/'report'/f'{label}.png'),full_page=False)
  findings.append(label+': no page overflow; screenshot captured')
 page.evaluate("document.querySelectorAll('details').forEach(d=>d.open=true)"); page.emulate_media(media='print'); page.pdf(path=str(root/'report/print-preview.pdf'),format='A4',print_background=True)
 assert not network,network
 browser.close()
 (root/'results/browser-verification.json').write_text(json.dumps({'offline':True,'runtime_network_requests':network,'checks':findings,'print':'PDF generated for inspection'},indent=2)+'\n')
 print('PASS: offline desktop/mobile layout and print PDF')
