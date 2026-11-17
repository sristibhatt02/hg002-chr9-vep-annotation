#!/usr/bin/env python3
"""Reproducible public Ensembl training annotation; Python standard library only."""
import argparse, csv, hashlib, html, json, os, re, sys, time, urllib.request, urllib.error
from pathlib import Path
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parent
INPUTS=[('var1',128328460,'TA','T'),('var2',128322349,'C','A'),('var3',128323079,'C','G'),('var4',128322917,'G','A'),('var5',128203516,'C','A')]
BASE='https://rest.ensembl.org'
OPTIONS={'canonical':1,'hgvs':1,'mane':1,'transcript_version':1}
SOURCES={'Training walkthrough':'https://training.ensembl.org/text_module/human_vep_walkthrough','Assembly walkthrough PDF':'https://training.ensembl.org/pdf_walkthrough.pdf','REST VEP documentation':BASE+'/documentation/info/vep_region_post','Sequence endpoint':BASE+'/documentation/info/sequence_region','Release metadata':BASE+'/documentation/info/data','Variant interpretation context':'https://www.genomicseducation.hee.nhs.uk/education/taught-courses/variant-interpretation/','Ensembl reuse terms':'https://www.ensembl.org/info/about/legal/disclaimer.html'}
def dump(p,x):
 p=ROOT/p; p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def sha(b): return hashlib.sha256(b).hexdigest()
def request(path,body=None):
 url=BASE+path
 print('Request '+url,flush=True)
 for attempt in range(3):
  try:
   req=urllib.request.Request(url,data=json.dumps(body).encode() if body else None,headers={'Content-Type':'application/json','Accept':'application/json','User-Agent':'public-training-vep-portfolio/1.0'})
   with urllib.request.urlopen(req,timeout=60) as r: return r.read()
  except urllib.error.HTTPError as e:
   if e.code not in (429,500,502,503,504) or attempt==2: raise
   delay=e.headers.get('Retry-After','2')
   try: delay=float(delay)
   except ValueError:
    from email.utils import parsedate_to_datetime
    delay=max(0,(parsedate_to_datetime(delay)-datetime.now(timezone.utc)).total_seconds())
   if delay>120: raise RuntimeError('Retry-After exceeds bounded wait; retry fetch later') from e
   time.sleep(delay)
  except (urllib.error.URLError,TimeoutError):
   if attempt==2: raise
   time.sleep(2**attempt)
def vcf():
 p=ROOT/'data/training-variants.grch38.vcf'
 expected='##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=9,length=138394717,assembly=GRCh38>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+''.join(f'9\t{pos}\t{i}\t{ref}\t{alt}\t.\t.\t.\n' for i,pos,ref,alt in INPUTS)
 if not p.exists(): p.write_text(expected)
 if p.read_text()!=expected: raise ValueError('VCF differs from the five fixed training records')
 return [line for line in expected.splitlines() if not line.startswith('#')]
def reconcile(raw):
 if not isinstance(raw,list): raise ValueError('VEP response must be a list')
 found={}; records={i:(pos,ref,alt) for i,pos,ref,alt in INPUTS}
 for obj in raw:
  tokens=obj.get('input','').split()
  if len(tokens)<5 or tokens[2] not in records: raise ValueError('Unmatched echoed input')
  i=tokens[2]; pos,ref,alt=records[i]
  if tokens[:5]!=['9',str(pos),i,ref,alt]: raise ValueError('Echoed identity mismatch')
  if i in found: raise ValueError('Duplicated response '+i)
  if obj.get('assembly_name')!='GRCh38' or str(obj.get('seq_region_name'))!='9': raise ValueError('Wrong assembly/chromosome')
  # REST trims the common anchor on the deletion. Echo identity is authoritative;
  # independently check this example's normalized forward-strand representation.
  expected=(pos+1,pos+1,'A/-') if i=='var1' else (pos,pos,ref+'/'+alt)
  if (obj.get('start'),obj.get('end'),obj.get('allele_string'))!=expected or obj.get('strand')!=1: raise ValueError('Unexpected normalized locus/allele '+i)
  found[i]=obj
 if set(found)!=set(records): raise ValueError('Missing inputs: '+str(set(records)-set(found)))
 return found
FIELDS=['row_id','variant_id','submitted_chromosome','submitted_position','ref','alt','returned_locus','returned_alleles','gene_symbol','gene_id','transcript_id','transcript_version','biotype','consequence_terms','mane_select','canonical','hgvsc','hgvsp']
def normalize(raw):
 found=reconcile(raw); rows=[]; other=[]; summaries=[]
 for i,pos,ref,alt in INPUTS:
  o=found[i]; ts=sorted(o.get('transcript_consequences',[]),key=lambda t:(t.get('transcript_id',''),json.dumps(t,sort_keys=True)))
  for n,t in enumerate(ts,1):
   vals=[f'{i}-t{n:03}',i,'9',pos,ref,alt,f"9:{o['start']}-{o['end']}",o['allele_string'],t.get('gene_symbol',''),t.get('gene_id',''),t.get('transcript_id',''),t.get('transcript_version',''),t.get('biotype',''),';'.join(t.get('consequence_terms',[])),t.get('mane_select',''),t.get('canonical',''),t.get('hgvsc',''),t.get('hgvsp','')]
   rows.append(dict(zip(FIELDS,vals)))
  for key in ('regulatory_feature_consequences','motif_feature_consequences','intergenic_consequences'):
   for feature in o.get(key,[]): other.append({'variant_id':i,'feature_type':key,'annotation':feature})
  summaries.append({'variant_id':i,'transcript_rows':len(ts),'distinct_transcripts':len({t['transcript_id'] for t in ts}),'consequence_terms':';'.join(sorted({c for t in ts for c in t.get('consequence_terms',[])})),'provider_most_severe':o.get('most_severe_consequence',''),'genes':';'.join(sorted({t.get('gene_symbol',t.get('gene_id','')) for t in ts}))})
 return rows,summaries,other

def fetch():
 (ROOT/'results/raw').mkdir(parents=True,exist_ok=True)
 lines=vcf(); checks=[]
 for i,pos,ref,alt in INPUTS:
  path=f'/sequence/region/human/9:{pos}..{pos+len(ref)-1}:1?content-type=application/json'
  b=request(path); o=json.loads(b)
  if o.get('seq','').upper()!=ref or 'GRCh38' not in o.get('id',''): raise ValueError('Reference mismatch or unverified build '+i+': '+str(o))
  checks.append({'variant_id':i,'expected':ref,'observed':o['seq'],'passed':True,'source':BASE+path,'response':o})
 body={'variants':lines,**OPTIONS}; endpoint='/vep/human/region'; raw=request(endpoint,body); normalize(json.loads(raw))
 metadata={}
 for key,path in [('data','/info/data'),('rest','/info/rest'),('software','/info/software'),('assembly','/info/assembly/homo_sapiens')]:
  try: metadata[key]=json.loads(request(path+'?content-type=application/json'))
  except Exception as e: metadata[key]={'not_reported':str(e)}
 stamp=datetime.now(timezone.utc).isoformat()
 manifest={'acquired_utc':stamp,'backend':'Ensembl REST batch VEP','endpoint':BASE+endpoint,'species':'Homo sapiens','assembly':'GRCh38','transcript_source':'Ensembl default; no transcript filtering','request':body,'input_sha256':sha((ROOT/'data/training-variants.grch38.vcf').read_bytes()),'output_sha256':sha(raw),'metadata':metadata,'vep_engine_version':'not reported','reference_checks':checks,'sources':SOURCES,'warnings':['Known variant IDs and clinical scores are not interpreted; co-location is not allele matching.'],'rights':'Public training examples, no patient metadata. Upstream data retain their own terms; review terms before redistribution.'}
 # Validate everything before touching the previous successful snapshot.
 (ROOT/'results/raw/vep-response.json').write_bytes(raw); dump('results/run-manifest.json',manifest); dump('results/raw/request.json',body)
 print('Fetched and validated all five inputs')
def tsv(path,rows,fields):
 with (ROOT/path).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
def load():
 m=json.loads((ROOT/'results/run-manifest.json').read_text()); raw=(ROOT/'results/raw/vep-response.json').read_bytes()
 if sha(raw)!=m['output_sha256'] or sha((ROOT/'data/training-variants.grch38.vcf').read_bytes())!=m['input_sha256']: raise ValueError('Snapshot hash mismatch')
 if len(m['reference_checks'])!=5 or not all(x['passed'] for x in m['reference_checks']): raise ValueError('Reference checks incomplete')
 return m,json.loads(raw)
def esc(v): return html.escape(str(v if v!='' else 'not reported'))
def table(headers,rows): return '<div class="scroll" tabindex="0" role="region" aria-label="Scrollable evidence table"><table><caption>'+esc(headers[0])+' and associated evidence</caption><thead><tr>'+''.join('<th scope="col">'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def choose(rows):
 def rank(r): return (not bool(r['mane_select']),not bool(r['canonical']),r['transcript_id'])
 selected=[sorted([r for r in rows if r['variant_id']==i],key=rank)[0] for i,*_ in INPUTS if any(r['variant_id']==i for r in rows)]
 if len(selected)<2: raise ValueError('Two transcript case studies unavailable')
 coding={'frameshift_variant','missense_variant','stop_gained','stop_lost','start_lost','inframe_deletion'}
 a=next((r for r in selected if coding.intersection(r['consequence_terms'].split(';'))),selected[0])
 b=next((r for r in selected if r['variant_id']!=a['variant_id'] and not coding.intersection(r['consequence_terms'].split(';'))),next(r for r in selected if r['variant_id']!=a['variant_id']))
 return [a,b]
def report():
 (ROOT/'report').mkdir(parents=True,exist_ok=True)
 vcf(); m,raw=load(); rows,summaries,other=normalize(raw)
 tsv('results/transcript-annotations.tsv',rows,FIELDS); tsv('results/variant-summary.tsv',summaries,list(summaries[0])); dump('results/other-features.json',other)
 cases=choose(rows); dump('results/case-studies.json',cases)
 explanations={'frameshift_variant':'The deletion changes the reading frame of this protein-coding transcript.','missense_variant':'The substitution changes a predicted amino acid in this transcript.','synonymous_variant':'The coding substitution leaves the predicted amino acid unchanged.','intron_variant':'The variant lies within an intron of this transcript; a splicing effect is not established by this label.','upstream_gene_variant':'The variant is upstream of this transcript; the label alone does not establish a regulatory effect.','downstream_gene_variant':'The variant is downstream of this transcript; the label alone does not establish a regulatory effect.'}
 sections=''
 for r in cases:
  alts=[x for x in rows if x['variant_id']==r['variant_id'] and x['consequence_terms']!=r['consequence_terms']]
  alt=alts[0] if alts else None
  reason='MANE Select annotation' if r['mane_select'] else ('Ensembl canonical annotation' if r['canonical'] else 'lexicographically first reported transcript; neither MANE Select nor canonical was available')
  sections+=f'<article><p class="eyebrow">Case study · {esc(r["variant_id"])}</p><h3>{esc(r["gene_symbol"])}: {esc(r["consequence_terms"])}</h3><p>{esc(" ".join(explanations.get(c,"The exact consequence term is reported below; its molecular effect depends on transcript context.") for c in r["consequence_terms"].split(";")))}</p>'+table(['Evidence field','Observed output'],[['Stable TSV row',r['row_id']],['GRCh38 submitted variant',f"9:{r['submitted_position']} {r['ref']} → {r['alt']}"],['Transcript',r['transcript_id']],['Version',r['transcript_version']],['Gene ID',r['gene_id']],['Exact consequences',r['consequence_terms']],['HGVSc',r['hgvsc']],['HGVSp',r['hgvsp']],['MANE Select',r['mane_select']],['Selection rationale',reason]])
  sections+='<p>'+ (f'Another transcript, <code>{esc(alt["transcript_id"])}</code> (TSV row {esc(alt["row_id"])}), reports <strong>{esc(alt["consequence_terms"])}</strong>.' if alt else 'No different consequence term set was reported across the other transcripts for this variant.')+'</p><aside><h4>Follow-up evidence needed</h4><p>No external clinical evidence lookup was performed. Ask: are there allele- and GRCh38-matched ClinVar submissions, with what review status and conflicts? What population frequency is appropriate to the proposed condition and inheritance? Is the gene–disease relationship established? Are functional assays, segregation and phenotype evidence available?</p></aside></article>'
 different=sum(len({r['consequence_terms'] for r in rows if r['variant_id']==s['variant_id']})>1 for s in summaries)
 content=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Five variants, multiple transcript consequences</title><style>
 :root{{color-scheme:light;--ink:#173236;--teal:#176b68}}*{{box-sizing:border-box}}body{{margin:0;background:#f4f7f5;color:var(--ink);font:16px/1.65 system-ui,sans-serif}}main{{max-width:1120px;margin:auto;padding:48px 28px}}header{{border-top:6px solid var(--teal);padding:28px 0}}h1{{font-size:clamp(2rem,5vw,3.6rem);line-height:1.1;max-width:850px;letter-spacing:-.035em}}h2{{font-size:1.7rem;margin-top:0}}h3{{font-size:1.35rem}}p{{max-width:78ch}}.eyebrow{{text-transform:uppercase;letter-spacing:.12em;font-size:.78rem;font-weight:750;color:var(--teal)}}section{{padding:32px 0;border-top:1px solid #c8d8d3}}article{{background:white;padding:26px;margin:24px 0;border:1px solid #d2dfda;border-radius:8px}}aside{{padding:12px 20px;background:#eef5f1;margin-top:24px}}a{{color:#075e69;text-underline-offset:3px}}a:focus-visible,.scroll:focus-visible,summary:focus-visible{{outline:3px solid #af530b;outline-offset:4px}}code{{font-family:ui-monospace,monospace;overflow-wrap:anywhere}}.scroll{{overflow:auto;max-width:100%}}table{{border-collapse:collapse;width:100%;font-size:.88rem}}caption{{text-align:left;color:#405d60;padding:10px 0}}th{{text-align:left;background:#e5eeea}}th,td{{padding:10px 12px;border-bottom:1px solid #d8e2de;vertical-align:top;overflow-wrap:anywhere}}details{{margin:18px 0}}details table{{min-width:680px}}summary{{cursor:pointer;font-weight:650}}.metrics{{font-size:1.1rem;color:var(--teal)}}@media(max-width:440px){{main{{padding:24px 16px}}article{{padding:16px}}th,td{{padding:8px}}}}@media print{{body{{background:white;font-size:10pt}}main{{padding:0;max-width:none}}.scroll{{overflow:visible}}table{{font-size:8pt;table-layout:fixed}}details table{{min-width:0}}article{{break-inside:avoid}}a{{color:inherit}}details{{display:block}}}}
 </style></head><body><main><header><p class="eyebrow">Educational analysis of public training variants</p><h1>Five variants, multiple transcript consequences</h1><p>How do predicted consequences differ across transcripts, and what evidence would be needed to assess clinical significance?</p><p>Homo sapiens · GRCh38 · {esc(m['backend'])}<br>Acquired {esc(m['acquired_utc'])}</p></header><section><h2>Findings</h2><p class="metrics">5 input variants · 5 accounted-for annotations · {len(rows)} transcript annotation rows</p><p>{different} of five inputs have more than one transcript consequence term set. The same genomic call can have different effects in different transcript models.</p><p>The selected examples report {esc(cases[0]['consequence_terms'])} for {esc(cases[0]['variant_id'])} and {esc(cases[1]['consequence_terms'])} for {esc(cases[1]['variant_id'])}. These are molecular predictions, not clinical classifications.</p></section><section><h2>Five-variant overview</h2>'''
 content+=table(['Variant','GRCh38 position','Submitted REF → ALT','Gene context','Transcript rows','Distinct transcripts','Provider most severe'],[[s['variant_id'],f'9:{pos}',f'{ref} → {alt}',s['genes'],s['transcript_rows'],s['distinct_transcripts'],s['provider_most_severe']] for s,(_,pos,ref,alt) in zip(summaries,INPUTS)])
 content+='<p>Counts measure transcript annotation rows for these five examples, not people or population prevalence. Multi-term annotations can contribute to more than one consequence term. Provider most-severe labels summarize across features and may differ from a selected transcript.</p></section><section><h2>Two transcript case studies</h2>'+sections+'</section><section><h2>Methods and quality control</h2><p>The fixed walkthrough records were expanded to eight-column VCF without changing coordinates, IDs or alleles. All five forward-strand REF intervals matched the official GRCh38 sequence, including both bases of var1. Batch REST annotations retain all returned transcripts; echoed input identity and normalized allele/locus are checked independently of response order. Non-transcript features are exported separately.</p>'
 content+=table(['Reference check','Submitted REF','Observed GRCh38 REF','Result'],[[c['variant_id'],c['expected'],c['observed'],'PASS'] for c in m['reference_checks']])
 content+='<p>Transcript selection prioritizes MANE Select, then Ensembl canonical, then stable transcript ID. Coding versus noncoding examples are preferred when present. Canonical/MANE flags select reporting context; they do not prove clinical significance.</p><p>Service/database metadata: <code>'+esc(json.dumps({k:({x:v for x,v in value.items() if x in ('assembly_name','assembly_accession','release','releases','not_reported')} if isinstance(value,dict) else value) for k,value in m['metadata'].items()},sort_keys=True))+'</code></p><p>VEP engine version: '+esc(m['vep_engine_version'])+'. Database releases and REST/software versions identify their own components, not a verified VEP engine version.</p><p>Run warnings: '+esc('; '.join(m['warnings']))+'</p></section><section><h2>Transcript evidence</h2><p><a href="../results/transcript-annotations.tsv">Complete transcript TSV</a> · <a href="../results/variant-summary.tsv">Variant summaries</a> · <a href="../results/other-features.json">Other feature annotations</a>. Missing optional fields are explicit; known co-located IDs are intentionally excluded from the transcript table.</p>'
 for i,*_ in INPUTS:
  subset=[r for r in rows if r['variant_id']==i]
  content+='<details><summary>'+esc(i)+f' · {len(subset)} transcript rows</summary>'+table(['Row ID','Transcript','Biotype','Consequence terms','HGVSc'],[[r[k] for k in ['row_id','transcript_id','biotype','consequence_terms','hgvsc']] for r in subset])+'</details>'
 content+='</section><section><h2>Limitations and sources</h2><p>This annotates supplied variant calls; it does not call variants from sequencing reads. Predicted molecular consequence alone cannot establish clinical significance. These five public teaching examples have no patient metadata and do not represent a population. No diagnostic, patient or NHS affiliation claim is made. Optional fields vary by backend and release. No ClinVar, population or gene–disease lookup was performed.</p><ul>'+''.join(f'<li><a href="{html.escape(url,quote=True)}">{esc(label)}</a></li>' for label,url in SOURCES.items())+'</ul><p>Local audit trail: <a href="../data/training-variants.grch38.vcf">VCF inputs</a> · <a href="../results/raw/vep-response.json">Untouched VEP response</a> · <a href="../results/raw/request.json">Exact request</a> · <a href="../results/run-manifest.json">Run manifest</a> · <a href="../README.md">README</a></p></section></main></body></html>'
 (ROOT/'report/index.html').write_text(content)
 print(f'Rendered {len(rows)} transcript rows; cases: '+', '.join(r['row_id'] for r in cases))
def verify():
 import unittest
 suite=unittest.defaultTestLoader.discover(str(ROOT/'tests')); result=unittest.TextTestRunner(verbosity=2).run(suite)
 if not result.wasSuccessful(): raise ValueError('Tests failed')
 m,raw=load(); rows,summary,_=normalize(raw)
 for p,expected in [('results/transcript-annotations.tsv',rows),('results/variant-summary.tsv',summary)]:
  with (ROOT/p).open() as f: actual=list(csv.DictReader(f,delimiter='\t'))
  if actual!=[{k:str(v) for k,v in r.items()} for r in expected]: raise ValueError('TSV reconciliation failed '+p)
 for r in json.loads((ROOT/'results/case-studies.json').read_text()):
  if r not in rows: raise ValueError('Case study not in raw-derived rows')
 before=sha((ROOT/'report/index.html').read_bytes()); report(); report()
 if sha((ROOT/'report/index.html').read_bytes())!=before: raise ValueError('Nondeterministic rebuild')
 from html.parser import HTMLParser
 class Check(HTMLParser):
  def __init__(self): super().__init__(); self.ids=[]
  def handle_starttag(self,tag,attrs):
   a=dict(attrs)
   if 'id' in a: self.ids.append(a['id'])
   if tag in ('script','iframe') or (tag=='link' and a.get('rel')=='stylesheet'): raise ValueError('Unexpected runtime dependency')
   if tag=='a' and not a.get('href','').startswith(('https:','#')):
    if not (ROOT/'report'/a['href']).resolve().is_file(): raise ValueError('Broken local link')
 c=Check(); c.feed((ROOT/'report/index.html').read_text()); assert len(c.ids)==len(set(c.ids))
 print('PASS: identity, hashes, TSV totals, cases, deterministic rebuild and local links')
def main():
 p=argparse.ArgumentParser(); p.add_argument('command',choices=['fetch','report','verify']); args=p.parse_args()
 try: globals()[args.command]()
 except Exception as e: print(f'FAILED: {e}',file=sys.stderr); sys.exit(1)
if __name__=='__main__': main()
