"""Verify guided explorer and portable report using actual saved data."""
import hashlib,json,time,shutil
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parent;URL='http://127.0.0.1:8000/analysis/report/index.html';SHOTS=ROOT/'report/screenshots';SHOTS.mkdir(exist_ok=True)
def main():
 data=json.loads((ROOT/'report/assets/dataset.json').read_text());checks=[];timings={};errors=[];remote=[]
 with sync_playwright() as p:
  browser=p.chromium.launch(headless=True,args=['--no-proxy-server']);context=browser.new_context(viewport={'width':1440,'height':1000});page=context.new_page()
  page.on('pageerror',lambda e:errors.append(str(e)));page.on('console',lambda m:errors.append(m.text) if m.type=='error' else None)
  context.route('**/*',lambda r:r.continue_() if r.request.url.startswith('http://127.0.0.1:8000/') else (remote.append(r.request.url),r.abort())[-1])
  started=time.perf_counter();page.goto(URL);page.wait_for_function("document.getElementById('filter-status').textContent.includes('1,307')");timings['initial_ready_ms']=round((time.perf_counter()-started)*1000,2)
  assert page.locator('#position-chart button').count()==40;assert page.locator('#scene-wrap').is_hidden();assert page.locator('#scene canvas').count()==0;checks.append('2D default with forty regions; 3D library loaded only on demand')
  def count(expected):page.wait_for_function('(n)=>document.getElementById("filter-status").textContent.includes("Filtered: "+n.toLocaleString("en-GB")+" alleles")',arg=expected)
  started=time.perf_counter();page.select_option('#type','indel');count(231);timings['immediate_type_filter_ms']=round((time.perf_counter()-started)*1000,2);assert 'short insertions or deletions' in page.locator('#filter-explanation').inner_text()
  # No Apply click is allowed in any of these regression checks.
  page.select_option('#term','intron_variant');count(168);page.select_option('#type','SNV');count(732)
  page.fill('#gene','coq');expected=sum(a['variant_type']=='SNV' and a['provider_most_severe']=='intron_variant' and 'COQ' in (a['genes']+';'+a['gene_ids']).upper() for a in data['alleles']);count(expected);assert 'COQ' in page.locator('#filter-explanation').inner_text();checks.append('Type, consequence and partial gene filters update automatically and combine correctly')
  page.click('#reset');count(1307);checks.append('Reset restores complete saved sample')
  page.locator('.region-controls summary').click();page.fill('#start','128300000');page.locator('#start').press('Tab');count(sum(a['position']>=128300000 for a in data['alleles']));page.click('#reset');count(1307);checks.append('Coordinate changes apply on change; reset restores boundaries')
  page.fill('#gene','NO_SUCH_GENE');count(0);assert 'No alleles match' in page.locator('#message').inner_text();page.click('#reset');count(1307);checks.append('Explicit empty state')
  page.locator('#position-chart button').first.focus();page.keyboard.press('Enter');page.wait_for_function('document.querySelector("#allele-detail table")!==null');assert 'In plain language' in page.locator('#allele-detail').inner_text();assert page.evaluate('document.activeElement.dataset.region');checks.append('2D bar selection, keyboard focus and plain-language allele explanation')
  opts=page.locator('#allele-select option').evaluate_all('(xs)=>xs.map(x=>x.value).filter(Boolean)');aid=opts[-1];expected=next(a['transcript_rows'] for a in data['alleles'] if a['allele_id']==aid);page.select_option('#allele-select',aid);page.wait_for_function('(a)=>document.querySelector("#allele-detail h3").textContent===a',arg=aid);page.wait_for_function('(n)=>document.querySelectorAll("#allele-detail table tbody tr").length===n',arg=expected);checks.append('Allele selection loads exact expected transcript rows')
  page.locator('#overview summary').click();page.locator('#bin-rows button').nth(2).click();page.wait_for_function('document.querySelector("#allele-detail table")!==null');checks.append('Exact coverage table provides equivalent region selection')
  page.click('#reset');count(1307);page.evaluate("window.scrollTo({top:0,behavior:'instant'})");page.screenshot(path=str(SHOTS/'desktop.png'));page.locator('#explorer').scroll_into_view_if_needed();page.screenshot(path=str(SHOTS/'explorer-controls.png'));page.locator('#position-chart').screenshot(path=str(SHOTS/'explorer-2d.png'))
  page.click('#view3d');page.wait_for_selector('#scene canvas');page.wait_for_timeout(1600);assert page.locator('#scene-wrap').is_visible();canvas=page.locator('#scene canvas');canvas.scroll_into_view_if_needed();box=canvas.bounding_box();hit=False
  for fx in [.35,.42,.49,.56,.63,.7]:
   for fy in [.3,.42,.54,.66,.78]:
    canvas.click(position={'x':box['width']*fx,'y':box['height']*fy});page.wait_for_timeout(20)
    if page.locator('#bin-detail').inner_text().startswith('chr9:'):hit=True;break
   if hit:break
  assert hit;checks.append('Optional 3D renders and geometry is selectable');page.wait_for_timeout(600);page.locator('#scene-wrap').screenshot(path=str(SHOTS/'explorer.png'))
  for control in ['#pause','#zoom-in','#zoom-out','#camera-reset']:page.click(control)
  page.select_option('#scale','sqrt');page.click('#view2d');assert page.locator('#scene-wrap').is_hidden();checks.append('Optional camera/scale controls and return to standard view')
  page.set_viewport_size({'width':390,'height':900});assert page.evaluate('document.documentElement.scrollWidth<=innerWidth');page.evaluate("window.scrollTo({top:0,behavior:'instant'})");page.screenshot(path=str(SHOTS/'mobile.png'));page.locator('#explorer').scroll_into_view_if_needed();page.screenshot(path=str(SHOTS/'mobile-explorer.png'));checks.append('390px layout has no page overflow')
  resources=page.evaluate('performance.getEntriesByType("resource").map(r=>({name:r.name,transferSize:r.transferSize}))')
  reduced=context.new_page();reduced.emulate_media(reduced_motion='reduce');reduced.goto(URL);reduced.wait_for_selector('#position-chart button');assert reduced.locator('#scene-wrap').is_hidden();reduced.close();checks.append('Reduced motion starts in standard view')
  disabled=p.chromium.launch(headless=True,args=['--no-proxy-server','--disable-webgl']);dp=disabled.new_page();dp.goto(URL);dp.wait_for_selector('#position-chart button');dp.click('#view3d');dp.wait_for_function("document.getElementById('message').textContent.includes('WebGL is unavailable')");dp.locator('#position-chart button').first.click();dp.wait_for_function('document.querySelector("#allele-detail table")!==null');disabled.close();checks.append('Actual WebGL failure retains full 2D selection and evidence')
  assert not errors,errors;assert not remote,remote
  # Shareable HTML must open completely offline and independently of project assets.
  portable=browser.new_context(offline=True,viewport={'width':1440,'height':1000});sp=portable.new_page();external=[];sp.on('request',lambda r:external.append(r.url) if r.url.startswith(('http:','https:')) else None);sp.goto((PROJECT/'reports/experiment-report.html').as_uri());assert sp.locator('svg').count()==5;assert sp.locator('h1').count()==1;assert not external;assert sp.evaluate('document.documentElement.scrollWidth<=innerWidth');sp.screenshot(path=str(SHOTS/'experiment-report.png'));sp.emulate_media(media='print');sp.pdf(path=str(PROJECT/'reports/experiment-report.pdf'),format='A4',print_background=True,display_header_footer=True,header_template='<span></span>',footer_template='<div style="font-family:Arial;font-size:8px;color:#526579;width:100%;margin:0 15mm;display:flex;justify-content:space-between"><span>HG002 Chromosome 9 Variant Annotation · Educational benchmark analysis</span><span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>',prefer_css_page_size=True);shutil.copyfile(PROJECT/'reports/experiment-report.pdf',ROOT/'report/experiment-report.pdf');portable.close();checks.append('Self-contained shareable HTML opens offline with five embedded scientific figures; PDF generated')
  result={'status':'passed','browser':browser.version,'viewport_checks':[1440,390],'console_errors':errors,'remote_requests':remote,'timings':timings,'resources':resources,'checks':checks};(ROOT/'results/browser-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='resources'},indent=2));browser.close()
if __name__=='__main__':main()
