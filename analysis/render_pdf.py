"""Render the standalone report offline with Chromium and verify its document structure."""
import hashlib,json,shutil
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'reports'
with sync_playwright() as p:
 browser=p.chromium.launch(headless=True,args=['--no-proxy-server']);context=browser.new_context(offline=True);page=context.new_page();page.goto((OUT/'experiment-report.html').as_uri());assert page.locator('svg').count()==5
 page.emulate_media(media='print');page.pdf(path=str(OUT/'experiment-report.pdf'),print_background=True,prefer_css_page_size=True,display_header_footer=True,header_template='<span></span>',footer_template='<div style="font-family:Arial;font-size:8px;color:#526579;width:100%;margin:0 15mm;display:flex;justify-content:space-between"><span>HG002 Chromosome 9 Variant Annotation · Educational benchmark analysis</span><span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>');browser.close()
shutil.copyfile(OUT/'experiment-report.pdf',ROOT/'report/experiment-report.pdf')
m=json.loads((OUT/'report-manifest.json').read_text());m['pdf_sha256']=hashlib.sha256((OUT/'experiment-report.pdf').read_bytes()).hexdigest();(OUT/'report-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print(OUT/'experiment-report.pdf')
