#!/usr/bin/env python3
"""GIAB HG002 bounded pilot, versioned exports and deterministic local report."""
import argparse, collections, csv, gzip, hashlib, html, json, os, platform, re, resource, shutil, subprocess, sys, tempfile, time, urllib.request, urllib.error
from pathlib import Path
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'cache'; RESULTS=ROOT/'results'; SCOPE=json.loads((ROOT/'scope.json').read_text())
GIAB='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/'
VCF_NAME='HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
BED_NAME='HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed'
FASTA_URL='https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz'
OPTIONS={'canonical':1,'hgvs':1,'mane':1,'transcript_version':1}
# Ensembl consequence ranking; primary source linked in README/data dictionary.
SEVERITY='transcript_ablation splice_acceptor_variant splice_donor_variant stop_gained frameshift_variant stop_lost start_lost transcript_amplification feature_elongation feature_truncation inframe_insertion inframe_deletion missense_variant protein_altering_variant splice_donor_5th_base_variant splice_region_variant splice_donor_region_variant splice_polypyrimidine_tract_variant incomplete_terminal_codon_variant start_retained_variant stop_retained_variant synonymous_variant coding_sequence_variant mature_miRNA_variant 5_prime_UTR_variant 3_prime_UTR_variant non_coding_transcript_exon_variant intron_variant NMD_transcript_variant non_coding_transcript_variant coding_transcript_variant upstream_gene_variant downstream_gene_variant TFBS_ablation TFBS_amplification TF_binding_site_variant regulatory_region_ablation regulatory_region_amplification regulatory_region_variant intergenic_variant sequence_variant'.split()
RANK={t:n for n,t in enumerate(SEVERITY)}
def consequence_class(term):
 if term in ('transcript_ablation','stop_gained','frameshift_variant','stop_lost','start_lost'): return 'Truncating / lost'
 if term.startswith('splice_'): return 'Splice related'
 if term in ('missense_variant','inframe_insertion','inframe_deletion','synonymous_variant','protein_altering_variant','coding_sequence_variant','start_retained_variant','stop_retained_variant','incomplete_terminal_codon_variant'): return 'Other coding'
 if 'UTR' in term: return 'UTR'
 if term in ('upstream_gene_variant','downstream_gene_variant'): return 'Flanking gene'
 if term=='intergenic_variant': return 'Intergenic'
 if 'regulatory' in term or 'TF' in term: return 'Regulatory / motif'
 return 'Other transcript'

def now(): return datetime.now(timezone.utc).isoformat()
def sha(data): return hashlib.sha256(data).hexdigest()
def filehash(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def save(p,x):
 p=Path(p); p.parent.mkdir(parents=True,exist_ok=True); p.write_text(json.dumps(x,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
def read(p): return json.loads(Path(p).read_text())
def tsv(p,rows,fields=None):
 rows=list(rows); fields=fields or list(rows[0])
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def request(url,body=None):
 for attempt in range(3):
  try:
   req=urllib.request.Request(url,data=json.dumps(body).encode() if body is not None else None,headers={'Content-Type':'application/json','Accept':'application/json','User-Agent':'HG002-bounded-pilot/1.0'})
   with urllib.request.urlopen(req,timeout=120) as r: return r.read()
  except urllib.error.HTTPError as e:
   if e.code not in (429,500,502,503,504) or attempt==2: raise
   delay=e.headers.get('Retry-After','3')
   try: delay=float(delay)
   except ValueError:
    from email.utils import parsedate_to_datetime
    delay=max(0,(parsedate_to_datetime(delay)-datetime.now(timezone.utc)).total_seconds())
   if delay>120: raise RuntimeError('Retry-After exceeds bounded retry; resume later') from e
   time.sleep(delay)
  except (urllib.error.URLError,TimeoutError):
   if attempt==2: raise
   time.sleep(2**attempt)
def download(url,p):
 p=Path(p)
 if p.exists(): return
 print('Download',url,flush=True)
 tmp=p.with_suffix(p.suffix+'.tmp')
 subprocess.run(['curl','--fail','--location','--retry','2','--max-time','300','--output',str(tmp),url],check=True)
 tmp.replace(p)
def acquire():
 import pysam
 CACHE.mkdir(exist_ok=True)
 for name,target in [(VCF_NAME+'.tbi','benchmark.vcf.gz.tbi'),(BED_NAME,'benchmark.bed')]: download(GIAB+name,CACHE/target)
 download(FASTA_URL,CACHE/'chr9.fa.gz')
 download(FASTA_URL.rsplit('/',1)[0]+'/CHECKSUMS',CACHE/'reference-CHECKSUMS')
 if not (CACHE/'chr9.fa').exists():
  with gzip.open(CACHE/'chr9.fa.gz','rb') as a,(CACHE/'chr9.fa').open('wb') as b: shutil.copyfileobj(a,b)
 pysam.faidx(str(CACHE/'chr9.fa'))
 if not (CACHE/'region.original.vcf').exists():
  f=pysam.TabixFile(GIAB+VCF_NAME,index=str(CACHE/'benchmark.vcf.gz.tbi'))
  with (CACHE/'region.original.vcf.tmp').open('w') as out:
   out.write('\n'.join(f.header)+'\n')
   for row in f.fetch('chr9',SCOPE['start_0based'],SCOPE['end_0based_exclusive']): out.write(row+'\n')
  (CACHE/'region.original.vcf.tmp').replace(CACHE/'region.original.vcf')
 provenance={'schema_version':'1.0','acquired_utc':now(),'release':'NIST v4.2.1','sample':'HG002','assembly':'GRCh38','chromosome_reference_accession':'NC_000009.12','reference_assembly_accession':'GCA_000001405.15 (GRCh38 base assembly; chromosome 9 unchanged in GRCh38.p14)','extraction':'HTSlib indexed remote Tabix query chr9:128000001-129000000; original VCF header retained','sources':{'vcf':GIAB+VCF_NAME,'index':GIAB+VCF_NAME+'.tbi','bed':GIAB+BED_NAME,'reference':FASTA_URL},'sha256':{p.name:filehash(p) for p in [CACHE/'benchmark.vcf.gz.tbi',CACHE/'benchmark.bed',CACHE/'chr9.fa.gz',CACHE/'chr9.fa',CACHE/'region.original.vcf']},'full_vcf_downloaded':False,'source_checksum_note':'SHA-256 recorded for acquired local bytes; full remote VCF not downloaded or claimed to have a full-file checksum. Official reference CHECKSUMS saved separately.'}
 if (CACHE/'source-manifest.json').exists():
  old=read(CACHE/'source-manifest.json')
  if old['sha256']==provenance['sha256']: return
 save(CACHE/'source-manifest.json',provenance)

def covered(start,end,bed): return any(a<=start and end<=b for a,b in bed)
def eligible(ref,alt): return bool(re.fullmatch('[ACGT]+',ref) and re.fullmatch('[ACGT]+',alt) and max(len(ref),len(alt))<=50 and ref!=alt)
def bed_regions():
 out=[]
 for line in (CACHE/'benchmark.bed').read_text().splitlines():
  c,a,b,*_=line.split('\t')
  if c not in ('9','chr9'): continue
  a,b=max(int(a),SCOPE['start_0based']),min(int(b),SCOPE['end_0based_exclusive'])
  if a<b: out.append((a,b))
 merged=[]
 for a,b in sorted(out):
  if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
  else: merged.append((a,b))
 return merged

def prepare():
 import pysam,pysam.bcftools
 start=time.perf_counter(); bed=bed_regions(); counts=collections.Counter(); excluded=[]; records=[]; origin={}
 fasta=pysam.FastaFile(str(CACHE/'chr9.fa')); a,b=SCOPE['start_0based'],SCOPE['end_0based_exclusive']
 def reject(key,reason,line): counts[reason]+=1; excluded.append({'original_allele_id':key,'reason':reason,'original_record':line})
 for num,line in enumerate((CACHE/'region.original.vcf').read_text().splitlines()):
  if line.startswith('#'): continue
  counts['original_records']+=1; p=line.split('\t'); chrom,pos,oid,ref,alts=p[:5]; pos=int(pos)
  gt=p[9].split(':')[p[8].split(':').index('GT')] if len(p)>9 and 'GT' in p[8].split(':') else None
  called=set(re.split(r'[/|]',gt)) if gt else set()
  for idx,alt in enumerate(alts.split(','),1):
   counts['ingested_alt_alleles']+=1; key=f'source-{num:06}-a{idx}'
   if p[6] not in ('PASS','.'): reject(key,'nonpass',line); continue
   if str(idx) not in called: reject(key,'not_present_in_sample_gt',line); continue
   if not eligible(ref,alt): reject(key,'unsupported_or_over_50bp',line); continue
   if not (a<=pos-1 and pos-1+len(ref)<=b): reject(key,'outside_full_ref_span',line); continue
   if not covered(pos-1,pos-1+len(ref),bed): reject(key,'outside_high_confidence',line); continue
   if fasta.fetch('9',pos-1,pos-1+len(ref)).upper()!=ref: reject(key,'reference_mismatch',line); continue
   counts['reference_matched_before_norm']+=1
   origin[key]={'original_allele_id':key,'original_record_id':oid,'original_alt_index':idx,'original_chromosome':chrom,'original_position':pos,'original_ref':ref,'original_alt':alt,'original_record':line}
   records.append(f'9\t{pos}\t{key}\t{ref}\t{alt}\t.\tPASS\t.')
 header='##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=9,length=138394717>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
 (CACHE/'split.vcf').write_text(header+'\n'.join(records)+'\n')
 args=['norm','--no-version','-f',str(CACHE/'chr9.fa'),'-c','e','-m','-any','-Ov',str(CACHE/'split.vcf')]
 normalized=pysam.bcftools.norm(*args[1:]); (CACHE/'bcftools-normalized.vcf').write_text(normalized)
 groups={}
 for line in normalized.splitlines():
  if line.startswith('#'): continue
  p=line.split('\t'); pos=int(p[1]); ref,alt=p[3:5]; key=p[2]; counts['normalized_alleles']+=1
  if not (a<=pos-1 and pos-1+len(ref)<=b) or not covered(pos-1,pos-1+len(ref),bed): reject(key,'normalized_outside_interval_or_hc',origin[key]['original_record']); continue
  if fasta.fetch('9',pos-1,pos-1+len(ref)).upper()!=ref: raise ValueError('Post-normalization reference mismatch')
  ident=f'9:{pos}:{ref}:{alt}'
  groups.setdefault(ident,{'allele_id':ident,'chromosome':'9','position':pos,'ref':ref,'alt':alt,'variant_type':'SNV' if len(ref)==len(alt)==1 else ('indel' if len(ref)!=len(alt) else 'MNV'),'origins':[]})['origins'].append(origin[key])
 alleles=sorted(groups.values(),key=lambda x:(x['position'],x['ref'],x['alt']))
 counts['deduplicated_eligible_alleles']=len(alleles); counts['duplicate_normalized_alleles']=counts['normalized_alleles']-counts['normalized_outside_interval_or_hc']-len(alleles)
 if len(alleles)>SCOPE['rest_max_alleles']: raise ValueError('Declared bounded REST budget exceeded; requires offline VEP')
 if counts['reference_mismatch']:
  tsv(CACHE/'failed-reference-exclusions.tsv',excluded,['original_allele_id','reason','original_record']); save(CACHE/'failed-reference-qc.json',dict(counts))
  raise ValueError('REF mismatches encountered; inspect cache/failed-reference-exclusions.tsv before acquisition')
 RESULTS.mkdir(exist_ok=True); tsv(RESULTS/'exclusions.tsv',excluded,['original_allele_id','reason','original_record']); save(CACHE/'alleles.json',alleles); save(CACHE/'prepare-qc.json',dict(counts))
 (CACHE/'normalized.vcf').write_text(header+'\n'.join(f"9\t{x['position']}\t{x['allele_id']}\t{x['ref']}\t{x['alt']}\t.\tPASS\t." for x in alleles)+'\n')
 provenance={'schema_version':'1.0','scope':SCOPE,'sources':read(CACHE/'source-manifest.json'),'tool_versions':{'python':platform.python_version(),'pysam':pysam.__version__,'embedded_bcftools':pysam.__samtools_version__},'normalization_argv':args,'reference_sha256':filehash(CACHE/'chr9.fa'),'normalized_vcf_sha256':filehash(CACHE/'normalized.vcf'),'elapsed_seconds':time.perf_counter()-start,'process_peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'platform':platform.platform(),'architecture':platform.machine(),'qc':dict(counts)}
 save(CACHE/'preparation-manifest.json',provenance); print('Prepared',dict(counts),flush=True)

def normalized_repr(x):
 pos=x['position']; ref=x['ref']; alt=x['alt']
 while ref and alt and ref[0]==alt[0]: pos+=1; ref=ref[1:]; alt=alt[1:]
 while ref and alt and ref[-1]==alt[-1]: ref=ref[:-1]; alt=alt[:-1]
 return pos,pos+len(ref)-1,(ref or '-')+'/'+(alt or '-')
def reconcile(raw,alleles):
 expected={x['allele_id']:x for x in alleles}; found={}
 if not isinstance(raw,list): raise ValueError('Not a list of VEP objects')
 for obj in raw:
  p=obj.get('input','').split()
  if len(p)<5 or p[2] not in expected: raise ValueError('Unmatched input')
  x=expected[p[2]]
  if p[:5]!=['9',str(x['position']),x['allele_id'],x['ref'],x['alt']]: raise ValueError('Echo identity mismatch')
  if x['allele_id'] in found: raise ValueError('Duplicate response')
  if obj.get('assembly_name')!='GRCh38' or str(obj.get('seq_region_name'))!='9' or obj.get('strand')!=1: raise ValueError('Assembly/strand mismatch')
  # A valid echo fixes the submitted identity; independently accept either trimmed
  # minimal representation or full VCF alleles if the service does not trim.
  actual=(obj.get('start'),obj.get('end'),obj.get('allele_string'))
  full=(x['position'],x['position']+len(x['ref'])-1,x['ref']+'/'+x['alt'])
  if actual not in (normalized_repr(x),full): raise ValueError(f'Unexpected returned representation {x["allele_id"]}: {actual}')
  found[x['allele_id']]=obj
 if set(found)!=set(expected): raise ValueError('Missing input annotations')
 return found

def annotate():
 live=json.loads(request('https://rest.ensembl.org/info/data?content-type=application/json'))
 if live!=read(ROOT/'docs/ensembl-release.json'): raise ValueError('REST database release changed; create a separately pinned snapshot instead of mixing releases')
 alleles=read(CACHE/'alleles.json'); prep=read(CACHE/'preparation-manifest.json'); started=time.perf_counter(); chunks=RESULTS/'chunks'; chunks.mkdir(exist_ok=True); index=[]
 fingerprint=sha(json.dumps({'normalized':prep['normalized_vcf_sha256'],'reference':prep['reference_sha256'],'options':OPTIONS,'engine':'REST human GRCh38','scope':SCOPE},sort_keys=True).encode())
 for n,offset in enumerate(range(0,len(alleles),100)):
  subset=alleles[offset:offset+100]; body={'variants':[f"9 {x['position']} {x['allele_id']} {x['ref']} {x['alt']} . . ." for x in subset],**OPTIONS}
  key=sha(json.dumps({'body':body,'fingerprint':fingerprint},sort_keys=True).encode()); p=chunks/(key+'.json'); meta=chunks/(key+'.manifest.json')
  if p.exists() and meta.exists() and read(meta)['sha256']==filehash(p):
   reconcile(read(p),subset); print('Resume valid chunk',n+1,flush=True)
  else:
   print('Annotating chunk',n+1,'/',(len(alleles)+99)//100,flush=True)
   b=request('https://rest.ensembl.org/vep/human/region',body); reconcile(json.loads(b),subset)
   tmp=p.with_suffix('.tmp'); tmp.write_bytes(b); tmp.replace(p)
   save(meta,{'fingerprint':fingerprint,'request':body,'acquired_utc':now(),'sha256':sha(b)})
  index.append({'file':p.name,'manifest':meta.name,'sha256':filehash(p)})
 acquired=[read(chunks/x['manifest'])['acquired_utc'] for x in index]
 inputs=RESULTS/'inputs'/fingerprint; inputs.mkdir(parents=True,exist_ok=True)
 for name in ['alleles.json','normalized.vcf','region.original.vcf','prepare-qc.json']:
  shutil.copyfile(CACHE/name,inputs/name)
 shutil.copyfile(RESULTS/'exclusions.tsv',inputs/'exclusions.tsv'); save(inputs/'bed.json',bed_regions())
 manifest={'input_directory':'inputs/'+fingerprint,'schema_version':'1.0','status':'bounded_pilot_complete; chromosome_run_blocked','backend':'Ensembl REST batch VEP (bounded pilot only)','endpoint':'https://rest.ensembl.org/vep/human/region','options':OPTIONS,'transcript_source':'Ensembl default, all reported transcripts','acquired_utc':max(acquired),'first_chunk_utc':min(acquired),'fingerprint':fingerprint,'chunks':index,'preparation':prep,'scope':SCOPE,'service_metadata':read(ROOT/'docs/ensembl-release.json'),'vep_engine_version':'not reported','alleles':len(alleles),'annotation_elapsed_seconds_this_invocation':time.perf_counter()-started,'process_peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'clinical_overlay':SCOPE['clinvar_status']}
 # A single atomic pointer promotes only complete, validated chunk collections.
 tmp=RESULTS/'snapshot.pending.json'; save(tmp,manifest); tmp.replace(RESULTS/'snapshot.json'); print('Promoted complete snapshot',len(alleles),flush=True)

def severity(terms):
 terms=list(terms)
 unknown=set(terms)-set(RANK)
 if unknown: raise ValueError('Unrecognized severity terms; update documented ranking: '+str(unknown))
 return min(terms,key=lambda t:RANK[t]) if terms else ''
def representative(rows):
 ranked=sorted(rows,key=lambda r:(not bool(r['mane_select']),not bool(r['canonical']),r['transcript_id'],r['row_id']))
 r=ranked[0]; tier='MANE Select' if r['mane_select'] else ('canonical' if r['canonical'] else 'fallback: transcript accession sort')
 ties=sum((bool(x['mane_select'])==bool(r['mane_select']) and (bool(x['canonical'])==bool(r['canonical']) if not r['mane_select'] else True)) for x in rows)
 return r,tier,ties

def current_raw():
 m=read(RESULTS/'snapshot.json'); inputs=RESULTS/m['input_directory']; alleles=read(inputs/'alleles.json')
 if filehash(inputs/'normalized.vcf')!=m['preparation']['normalized_vcf_sha256']: raise ValueError('Prepared input differs from snapshot')
 objs=[]
 for x in m['chunks']:
  p=RESULTS/'chunks'/x['file']
  if filehash(p)!=x['sha256']: raise ValueError('Raw hash mismatch')
  objs.extend(read(p))
 return m,alleles,reconcile(objs,alleles)

def analyze():
 m,alleles,raw=current_raw(); clinical=read(RESULTS/'clinvar-overlay.json') if (RESULTS/'clinvar-overlay.json').exists() else None; out=Path(tempfile.mkdtemp(prefix='.staging-exports-',dir=RESULTS)); transcript=[]; other=[]; summary=[]; pair_rows=[]; selected={}; missing=collections.Counter()
 for x in alleles:
  aid=x['allele_id']; obj=raw[aid]; rows=[]
  for n,t in enumerate(sorted(obj.get('transcript_consequences',[]),key=lambda t:(t.get('gene_id',''),t.get('transcript_id',''),json.dumps(t,sort_keys=True))),1):
   tid=t.get('transcript_id',''); version=t.get('transcript_version','') or (tid.rsplit('.',1)[1] if re.search(r'\.\d+$',tid) else '')
   r={'row_id':f'{aid}/t{n:04}','allele_id':aid,'gene_symbol':t.get('gene_symbol',''),'gene_id':t.get('gene_id',''),'transcript_id':tid,'transcript_version':version,'biotype':t.get('biotype',''),'consequence_terms':';'.join(t.get('consequence_terms',[])),'most_severe_transcript_term':severity(t.get('consequence_terms',[])),'mane_select':t.get('mane_select',''),'canonical':t.get('canonical',0),'hgvsc':t.get('hgvsc',''),'hgvsp':t.get('hgvsp',''),'variant_allele':t.get('variant_allele','')}
   for k in ['gene_symbol','transcript_id','mane_select','hgvsc','hgvsp']:
    if not r[k]: missing[k]+=1
   rows.append(r); transcript.append(r)
  for key in ('regulatory_feature_consequences','motif_feature_consequences','intergenic_consequences'):
   for feature in obj.get(key,[]): other.append({'allele_id':aid,'feature_type':key,'annotation':feature})
  genes=collections.defaultdict(list)
  for r in rows: genes[r['gene_id'] or '(not reported)'].append(r)
  selected[aid]=[]
  for gid,rs in sorted(genes.items()):
   pick,tier,ties=representative(rs); most=severity(t for r in rs for t in r['consequence_terms'].split(';') if t); selected[aid].append(pick)
   pair_rows.append({'allele_id':aid,'gene_id':gid,'gene_symbol':pick['gene_symbol'],'transcript_rows':len(rs),'distinct_term_sets':len({tuple(sorted(r['consequence_terms'].split(';'))) for r in rs}),'different_across_transcripts':int(len({tuple(sorted(r['consequence_terms'].split(';'))) for r in rs})>1),'selected_row_id':pick['row_id'],'selection_tier':tier,'selection_ties':ties,'selected_consequences':pick['consequence_terms'],'selected_most_severe':pick['most_severe_transcript_term'],'all_transcript_most_severe':most,'selected_misses_more_severe':int(bool(most) and RANK[most]<RANK[pick['most_severe_transcript_term']])})
  fields={k:v for k,v in x.items() if k!='origins'}
  summary.append({**fields,'provider_most_severe':obj.get('most_severe_consequence',''),'consequence_class':consequence_class(obj.get('most_severe_consequence','')),'transcript_rows':len(rows),'genes':';'.join(sorted({r['gene_symbol'] or r['gene_id'] for r in rows})),'gene_ids':';'.join(sorted({r['gene_id'] for r in rows})),'transcript_terms':';'.join(sorted({t for r in rows for t in r['consequence_terms'].split(';')})),'clinvar_status':('exact match in ClinVar 2026-11-08' if clinical['hits'].get(aid) else 'no exact eligible match in ClinVar 2026-11-08') if clinical else m['clinical_overlay']})
 spectrum=[]
 for unit,items,field in [('alleles',summary,'provider_most_severe'),('transcript_rows',transcript,'most_severe_transcript_term')]:
  c=collections.Counter((r.get('variant_type','all'),r[field]) for r in items)
  for (typ,term),count in sorted(c.items()): spectrum.append({'unit':unit,'variant_type':typ,'term':term,'count':count,'denominator':len(items),'type_denominator':sum(r.get('variant_type','all')==typ for r in items)})
 multi=collections.Counter(t for r in transcript for t in r['consequence_terms'].split(';'))
 tsv(out/'transcript-term-counts.tsv',[{'term':t,'transcript_rows_containing_term':v,'denominator_transcript_rows':len(transcript),'multi_counted':'yes'} for t,v in sorted(multi.items())])
 bed=read(RESULTS/m['input_directory']/'bed.json'); bins=[]; width=SCOPE['bin_width_bp']
 for a in range(SCOPE['start_0based'],SCOPE['end_0based_exclusive'],width):
  b=min(a+width,SCOPE['end_0based_exclusive']); cov=sum(max(0,min(b,y)-max(a,x)) for x,y in bed); members=[x for x in summary if a<=x['position']-1<b]
  c=collections.Counter(x['provider_most_severe'] for x in members)
  bins.append({'bin_id':f'9:{a+1}-{b}','start_0based':a,'end_0based_exclusive':b,'high_confidence_bases':cov,'alleles':len(members),'alleles_per_million_hc_bases':round(len(members)*1e6/cov,4) if cov else None,'terms':dict(sorted(c.items()))})
 gene_rows=[]
 for gid in sorted({p['gene_id'] for p in pair_rows}):
  ps=[p for p in pair_rows if p['gene_id']==gid]; gene_rows.append({'gene_id':gid,'gene_symbol':ps[0]['gene_symbol'],'distinct_alleles':len({p['allele_id'] for p in ps}),'transcript_rows':sum(p['transcript_rows'] for p in ps),'pairs_differing':sum(p['different_across_transcripts'] for p in ps),'pairs_selected_misses_more_severe':sum(p['selected_misses_more_severe'] for p in ps)})
 candidates=[]; used=set()
 def pick_case(label,predicate):
  found=next((x for x in summary if x['allele_id'] not in used and predicate(x)),None)
  if found:
   aid=found['allele_id']; used.add(aid); rs=[r for r in transcript if r['allele_id']==aid]; picks=selected[aid]
   pair=next((p for p in pair_rows if p['allele_id']==aid and p['selected_misses_more_severe']),None) or next((p for p in pair_rows if p['allele_id']==aid and p['different_across_transcripts']),None)
   r=next((r for r in picks if pair and r['row_id']==pair['selected_row_id']),picks[0] if picks else None)
   if label=='Coding SNV': r=next((v for v in picks if v['most_severe_transcript_term'] in ('missense_variant','stop_gained','synonymous_variant','stop_lost','start_lost')),r)
   if label=='Indel' and picks: r=min(picks,key=lambda v:RANK[v['most_severe_transcript_term']])
   if label=='Missing HGVS evidence': r=next((v for v in picks if not v['hgvsc']),r)
   alt=next((r2 for r2 in sorted(rs,key=lambda v:RANK[v['most_severe_transcript_term']]) if r and r2['gene_id']==r['gene_id'] and set(r2['consequence_terms'].split(';'))!=set(r['consequence_terms'].split(';'))),None)
   candidates.append({'label':label,'allele':found,'selected':r,'alternative':alt,'selection':next((p for p in pair_rows if r and p['selected_row_id']==r['row_id']),None),'clinvar_matches':clinical['hits'].get(aid,[]) if clinical else [],'source_raw_chunk':next(z['file'] for z in m['chunks'] if any(o['input'].split()[2]==aid for o in read(RESULTS/'chunks'/z['file'])))})
 coding_pool=[x for x in summary if x['variant_type']=='SNV' and x['provider_most_severe']=='missense_variant']
  # Actual matched evidence is preferred for a teaching case, never for cohort totals.
 matched_coding={x['allele_id'] for x in coding_pool if clinical and clinical['hits'].get(x['allele_id'])}
 coding_ids=matched_coding or {x['allele_id'] for x in coding_pool}
 pick_case('Coding SNV',lambda x:x['allele_id'] in coding_ids)
 coding_indels={x['allele_id'] for x in summary if x['variant_type']=='indel' and x['provider_most_severe'] in ('frameshift_variant','inframe_deletion','inframe_insertion')}
 pick_case('Indel',lambda x:x['allele_id'] in coding_indels if coding_indels else x['variant_type']=='indel')
 missed_pairs=[p for p in pair_rows if p['selected_misses_more_severe']]
 best_rank=min((RANK[p['all_transcript_most_severe']] for p in missed_pairs),default=999)
 differing={p['allele_id'] for p in missed_pairs if RANK[p['all_transcript_most_severe']]==best_rank}
 if not differing: differing={p['allele_id'] for p in pair_rows if p['different_across_transcripts']}
 pick_case('Transcript sensitivity',lambda x:x['allele_id'] in differing)
 no_hgvs={r['allele_id'] for group in selected.values() for r in group if not r['hgvsc']}
 pick_case('Missing HGVS evidence',lambda x:x['allele_id'] in no_hgvs)
 tsv(out/'allele-summary.tsv',summary); tsv(out/'transcripts.tsv',transcript); tsv(out/'transcript-selection.tsv',pair_rows); tsv(out/'gene-summary.tsv',gene_rows); tsv(out/'consequence-spectrum.tsv',spectrum)
 tsv(out/'bins.tsv',[{**{k:v for k,v in b.items() if k!='terms'},'term_counts_json':json.dumps(b['terms'],sort_keys=True)} for b in bins]); save(out/'other-features.json',other); save(out/'cases.json',candidates)
 tsv(out/'origin-map.tsv',[{'allele_id':x['allele_id'],**o} for x in alleles for o in x['origins']])
 inputs=RESULTS/m['input_directory']; shutil.copyfile(inputs/'normalized.vcf',out/'normalized.vcf'); shutil.copyfile(inputs/'region.original.vcf',out/'original-region.vcf'); shutil.copyfile(inputs/'exclusions.tsv',out/'exclusions.tsv')
 (out/'high-confidence.bed').write_text(''.join(f'chr9\t{a}\t{b}\n' for a,b in bed))
 metrics={'input_alleles':len(alleles),'annotated_alleles':len(raw),'transcript_rows':len(transcript),'unannotated_alleles':len(alleles)-len(raw),'no_transcript_alleles':sum(not x['transcript_rows'] for x in summary),'allele_gene_pairs':len(pair_rows),'pairs_differing':sum(p['different_across_transcripts'] for p in pair_rows),'pairs_selected_misses_more_severe':sum(p['selected_misses_more_severe'] for p in pair_rows),'fallback_pairs':sum(p['selection_tier'].startswith('fallback') for p in pair_rows),'tied_pairs':sum(p['selection_ties']>1 for p in pair_rows),'high_confidence_bases':sum(b-a for a,b in bed),'missing_optional_transcript_fields':dict(missing),'variant_types':dict(collections.Counter(x['variant_type'] for x in summary))}
 if clinical:
  metrics['clinvar']=clinical['qc']; shutil.copyfile(RESULTS/'clinvar-matches.tsv',out/'clinvar-matches.tsv'); save(out/'clinvar-overlay.json',clinical)
 save(out/'metrics.json',metrics)
 # Compact overview plus per-allele transcript details in separate local assets.
 details=Path(tempfile.mkdtemp(prefix='.staging-details-',dir=RESULTS))
 for aid in raw:
  save(details/(sha(aid.encode())[:20]+'.json'),{'allele':next(x for x in summary if x['allele_id']==aid),'transcripts':[r for r in transcript if r['allele_id']==aid],'clinvar_matches':clinical['hits'].get(aid,[]) if clinical else [],'selected':selected[aid]})
 compact={'schema_version':'1.0','term_classes':{t:consequence_class(t) for t in SEVERITY},'scope':{**SCOPE,'clinvar_status':'available: exact normalized overlap with ClinVar 2026-11-08' if clinical else SCOPE['clinvar_status']},'clinical_overlay':{'source':clinical['source'],'qc':clinical['qc']} if clinical else None,'manifest':{k:m[k] for k in ('acquired_utc','backend','service_metadata','vep_engine_version')},'metrics':metrics,'bins':bins,'spectrum':spectrum,'cases':candidates,'genes':gene_rows,'alleles':[{**x,'detail_file':'assets/details/'+sha(x['allele_id'].encode())[:20]+'.json'} for x in summary]}
 save(out/'report-dataset.json',compact)
 save(out/'export-manifest.json',{'schema_version':'1.0','snapshot_fingerprint':m['fingerprint'],'files':{p.name:filehash(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='export-manifest.json'},'metrics':metrics})
 destination=RESULTS/'exports'; backup=RESULTS/'exports.previous'
 if backup.exists(): shutil.rmtree(backup)
 if destination.exists(): destination.rename(backup)
 try: out.rename(destination)
 except Exception:
  if backup.exists(): backup.rename(destination)
  raise
 target=ROOT/'report/assets/details'; previous=ROOT/'report/assets/details.previous'
 if previous.exists(): shutil.rmtree(previous)
 if target.exists(): target.rename(previous)
 details.rename(target); shutil.copyfile(destination/'report-dataset.json',ROOT/'report/assets/dataset.json')
 out=destination
 save(RESULTS/'run-manifest.json',{'schema_version':'1.0','annotation_snapshot':m,'clinical_overlay':{'source':clinical['source'],'qc':clinical['qc']} if clinical else None,'exports_sha256':filehash(out/'export-manifest.json'),'scope_status':'completed 1 Mb pilot; full chromosome blocked'})
 print('Analyzed',metrics,flush=True)


def render():
 data=read(ROOT/'report/assets/dataset.json'); m=data['metrics']; s=data['scope']; e=html.escape
 def table(headers,rows): return '<div class="table-scroll" tabindex="0"><table><thead><tr>'+''.join('<th scope="col">'+e(str(x))+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+e(str(v))+'</td>' for v in r)+'</tr>' for r in rows)+'</tbody></table></div>'
 facts=f"<li><strong>{m['input_alleles']:,} DNA changes were analyzed.</strong> These are existing calls from one public benchmark sample, within a fixed one-million-base region. They are not patient diagnoses.</li><li><strong>The same change can have different predicted effects.</strong> A gene can produce different RNA versions, called transcripts. Consequence sets differed across those versions in {m['pairs_differing']:,} of {m['allele_gene_pairs']:,} variant–gene comparisons.</li><li><strong>A standard transcript does not show every annotation.</strong> In {m['pairs_selected_misses_more_severe']:,} comparisons, another transcript had a more severe molecular term than the representative. This is why the full transcript evidence is retained; it does not establish clinical risk.</li>"
 cases=''
 for case in data['cases']:
  x=case['allele']; r=case['selected']; alt=case['alternative']; selection=case['selection']
  fields=[['Allele (GRCh38)',x['allele_id']],['Provider severity summary',x['provider_most_severe']],['Gene',(r['gene_symbol'] or r['gene_id']) if r else 'No transcript'],['Selected transcript',r['transcript_id'] if r else 'Not reported'],['Selected consequence',r['consequence_terms'] if r else 'Not reported'],['Selection rationale',selection['selection_tier'] if selection else 'No eligible transcript'],['HGVSc',r['hgvsc'] or 'Not reported' if r else 'Not reported'],['HGVSp',r['hgvsp'] or 'Not reported' if r else 'Not reported'],['Selected source row',r['row_id'] if r else 'Not reported'],['Alternative transcript',alt['transcript_id'] if alt else 'No different same-gene term set'],['Alternative consequence',alt['consequence_terms'] if alt else 'Not applicable'],['Alternative source row',alt['row_id'] if alt else 'Not applicable']]
  clinical_rows=[[h['clinvar_variation_id'],h['classification_summary'],h['review_status'],h['conditions']] for h in case.get('clinvar_matches',[])]
  clinical_html=('<p>Exact normalized ClinVar match, release 2026-11-08:</p>'+table(['Variation ID','Submitted classification summary','Review status','Conditions'],clinical_rows)+''.join('<a href="'+e(h['url'])+'">ClinVar '+e(h['clinvar_variation_id'])+'</a> ' for h in case.get('clinvar_matches',[]))) if clinical_rows else '<p>No exact eligible match in ClinVar release 2026-11-08. This is not evidence of benignity.</p>'
  cases+='<article class="case"><p class="kicker">'+e(case['label'])+'</p><h3>'+e(x['allele_id'])+'</h3>'+table(['Evidence','Observed output'],fields)+clinical_html+'<p><a href="../results/chunks/'+e(case['source_raw_chunk'])+'">Raw VEP chunk</a> · <a href="../results/exports/transcripts.tsv">Full transcript rows</a></p><p class="note">Clinical significance is not established. Follow-up requires allele/build-matched assertions and review status, condition-appropriate frequency, gene–disease validity, phenotype, functional and segregation evidence.</p></article>'
 qc=dict(read(RESULTS/'snapshot.json')['preparation']['qc']); qc.update({'annotation_response_objects':m['annotated_alleles'],'export_transcript_rows':m['transcript_rows'],'alleles_without_transcript_overlap':m['no_transcript_alleles']})
 exports=['clinvar-matches.tsv','clinvar-overlay.json','allele-summary.tsv','transcripts.tsv','transcript-selection.tsv','gene-summary.tsv','consequence-spectrum.tsv','transcript-term-counts.tsv','bins.tsv','origin-map.tsv','original-region.vcf','normalized.vcf','exclusions.tsv','cases.json','metrics.json','export-manifest.json']
 spectrum_counts=collections.Counter()
 for x in data['spectrum']:
  if x['unit']=='alleles': spectrum_counts[x['term']]+=x['count']
 max_count=max(spectrum_counts.values(),default=1)
 bars='<figure><figcaption>Distinct alleles by provider most-severe consequence; denominator '+str(m['input_alleles'])+' HG002 alleles. One category per allele.</figcaption>'+''.join('<div class="spectrum-row"><span>'+e(term)+'</span><div class="bar-track"><div class="bar-fill" style="width:'+str(count*100/max_count)+'%"></div></div><strong>'+str(count)+'</strong></div>' for term,count in spectrum_counts.most_common())+'</figure>'
 overview='<h2>Observed consequences</h2>'+bars+table(['Unit','Type','Most-severe term','Count','All-allele denominator','Type denominator'],[[x['unit'],x['variant_type'],x['term'],x['count'],x['denominator'],x['type_denominator']] for x in data['spectrum'] if x['unit']=='alleles'])+'<details><summary>Transcript-level spectrum: '+str(m['transcript_rows'])+' annotation rows</summary><p>One most-severe term per transcript row; these rows are distinct from allele counts. The separate multi-term export counts each row in every applicable term.</p>'+table(['Transcript most-severe term','Rows','Transcript-row denominator'],[[x['term'],x['count'],x['denominator']] for x in data['spectrum'] if x['unit']=='transcript_rows'])+'</details>'
 methods=f'<h2>Methods, QC and reproducibility</h2><p>HG002 is one public benchmark sample, not a population cohort. GIAB NIST v4.2.1 GRCh38 VCF and matching noinconsistent BED were queried with a tabix index. The 1 Mb interval was fixed before examining consequences. Full anchored REF spans must be contained in the interval and a high-confidence BED region before and after normalization. VCF positions are 1-based; BED and bin boundaries are 0-based half-open.</p><p>ALT alleles present in the sample genotype were split, reference-checked and normalized with embedded bcftools 1.21 (pysam 0.23.3) and the official Ensembl GRCh38 chromosome 9 FASTA. SNVs, indels and MNVs with alleles ≤50 bp were eligible. Duplicate normalized chromosome/position/REF/ALT identities were merged with all original allele links retained.</p><p>Annotations use genuine bounded REST batches, all transcripts retained, with HGVS, MANE, canonical and transcript versions requested. Chunks resume only with matching input/settings/reference fingerprints and raw hashes. Transcript selection is per allele–gene pair: MANE Select, then canonical, then transcript accession order; ties are resolved deterministically and counted. Ensembl database release is '+e(str(data['manifest']['service_metadata']))+'. VEP engine version: not separately reported.</p>'+table(['Processing step or exclusion','Count'],sorted(qc.items()))
 limitations='<h2>Scope and limitations</h2><p class="note">Partial chromosome coverage: chr9:128,000,001–129,000,000 (GRCh38). The rest of chromosome 9 was not annotated. This pilot must not be described as a chromosome-wide or whole-genome analysis.</p><p>Offline VEP is blocked: the standard release 116 cache is listed at 26G compressed versus approximately 22 GiB free disk before extraction, and Docker has no running daemon. REST is used only for this bounded pilot.</p><p>ClinVar overlay uses the fixed 2026-11-08 GRCh38 release. It found exact normalized matches for '+str(m.get('clinvar',{}).get('hg002_exact_matched_alleles',0))+' of '+str(m['input_alleles'])+' alleles; '+str(m.get('clinvar',{}).get('hg002_unmatched_alleles',0))+' have no exact eligible match in this release. All chromosome 9 ClinVar alleles were queried, with the same ACGT/50 bp eligibility and reference normalization. Classification strings, review status, conditions and conflicts are retained as source assertions. No match is not evidence of benignity; co-located VEP IDs are not exact ClinVar matches. Predicted molecular consequence does not establish clinical significance. No disease hotspots, population frequencies, patient outcomes, variant-calling accuracy or clinical classifications are inferred.</p>'
 future='<h2>Future scope</h2><ul><li>Whole chromosome and genome: provision a compatible pinned offline engine/cache/reference and sufficient storage, then measure chromosome runtime and QC before parallelization.</li><li>Additional public samples: use matched annotations and shared high-confidence regions; population claims require an appropriate cohort design.</li><li>Formal calling benchmark: compare an independent query callset to matched truth VCF/BED using a representation-aware comparator before estimating precision/recall.</li><li>Clinical evidence: normalize a dated ClinVar release and population dataset consistently; retain exact allele matches, conditions, conflicts and review status.</li><li>Release drift: re-annotate a fixed sentinel panel across pinned releases and audit transcript and consequence changes.</li></ul>'
 sources='<h2>Sources and audit trail</h2><ul><li><a href="https://www.nist.gov/programs-projects/genome-bottle">NIST GIAB</a></li><li><a href="'+GIAB+VCF_NAME+'">Fixed GIAB v4.2.1 VCF</a> · <a href="'+GIAB+BED_NAME+'">Matching BED</a></li><li><a href="'+FASTA_URL+'">Ensembl GRCh38 chromosome 9 reference</a></li><li><a href="https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html">Ensembl consequence severity order</a></li><li><a href="https://rest.ensembl.org/documentation/info/vep_region_post">REST VEP</a> · <a href="https://static.ensembl.org/info/docs/tools/vep/script/vep_options.html">Offline VEP options</a></li><li><a href="https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/">ClinVar downloads</a> · <a href="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20261108.vcf.gz">Fixed 2026-11-08 GRCh38 release</a></li><li><a href="../results/snapshot.json">Snapshot manifest</a> · <a href="../../README.md">README</a> · <a href="../DATA_DICTIONARY.md">Data dictionary</a></li></ul><p>Exports: '+ ' · '.join('<a href="../results/exports/'+x+'">'+x+'</a>' for x in exports)+'</p>'
 continuity='<h2>Training validation continuity</h2><p>The separate five teaching variants and 216 transcript rows still pass regression. On COQ4 ENST00000300452.8, var3 reports missense_variant (c.134C&gt;G; p.Ser45Cys), whereas var1 reports intron_variant (c.402+2580del). Those examples are excluded from all HG002 totals. <a href="../../validation/report/index.html">Open Training validation report</a>.</p>'
 body=f'<header><p class="kicker">HG002 Chromosome 9 Variant Annotation · educational benchmark analysis</p><h1>One DNA change can have different transcript effects.</h1><p class="lede">We studied existing DNA calls from HG002, a public benchmark sample, to ask how predicted effects change across versions of a gene. Completed scope: GRCh38 chr9:128,000,001–129,000,000. The full chromosome run remains blocked.</p><p class="meta">GIAB v4.2.1 · Ensembl REST · acquired {e(data["manifest"]["acquired_utc"])}</p><ul class="findings">{facts}</ul></header>'
 explorer=(ROOT/'templates/explorer.html').read_text()
 sections='<section id="consequences">'+overview+f'<h3>Transcript selection sensitivity</h3><p>Unit: allele–gene pairs with transcript rows. {m["pairs_differing"]:,}/{m["allele_gene_pairs"]:,} pairs differ across transcripts; {m["pairs_selected_misses_more_severe"]:,} select a less severe term than another same-gene transcript. {m["fallback_pairs"]:,} pairs require accession-order fallback; {m["tied_pairs"]:,} pairs have ties at the preferred tier. These comparisons use molecular severity, not clinical pathogenicity.</p><p><a href="../results/exports/transcript-selection.tsv">Every selection and comparison</a> · <a href="../results/exports/transcript-term-counts.tsv">Multi-counted transcript term spectrum</a></p></section><section id="cases"><h2>Traceable cases</h2><div class="cases">'+cases+'</div></section><section>'+continuity+'</section><section id="methods">'+methods+'</section><section>'+limitations+'</section><section>'+future+'</section><section>'+sources+'</section>'
 shell=lambda inner,scripts: '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HG002 · HG002 Chromosome 9 Variant Annotation</title><link rel="stylesheet" href="assets/style.css?v=3"></head><body><a class="skip" href="#main">Skip to content</a><nav aria-label="Main navigation"><a href="#main">HG002 / GRCh38</a><a href="#explorer">Explorer</a><a href="experiment-report.html">Full report</a><a href="experiment-report.pdf">PDF</a></nav><main id="main">'+inner+'</main>'+scripts+'</body></html>'
 (ROOT/'report/index.html').write_text(shell(body+'<div class="report-links"><a href="experiment-report.html">Read the complete experiment report</a><a href="experiment-report.pdf">Download the report as PDF</a></div>'+explorer+'<details class="technical-details"><summary>Technical results, case-study source rows and reproducibility</summary>'+sections+'</details>','<script type="module" src="assets/app.js?v=3"></script>'))
 static_bins='<section id="explorer"><h2>Static genomic distribution</h2>'+table(['GRCh38 bin','Distinct alleles','HC bases','Alleles / million HC bases'],[[b['bin_id'],b['alleles'],b['high_confidence_bases'],b['alleles_per_million_hc_bases'] if b['high_confidence_bases'] else 'Not available'] for b in data['bins']])+'</section>'
 (ROOT/'report/summary.html').write_text(shell(body+static_bins+sections,''))
 print('Rendered interactive report and standalone summary')

def verify():
 m,alleles,raw=current_raw(); clinical=read(RESULTS/'clinvar-overlay.json') if (RESULTS/'clinvar-overlay.json').exists() else None; out=RESULTS/'exports'; metrics=read(out/'metrics.json'); dataset=read(ROOT/'report/assets/dataset.json'); exp=read(out/'export-manifest.json')
 for name,digest in exp['files'].items():
  if filehash(out/name)!=digest: raise ValueError('Export changed: '+name)
 assert metrics==dataset['metrics']; assert metrics['input_alleles']==len(raw)==len(alleles)
 with (out/'transcripts.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
 assert len(rows)==metrics['transcript_rows']; ids={r['row_id'] for r in rows}; assert len(ids)==len(rows)
 assert sum(b['alleles'] for b in dataset['bins'])==len(alleles)
 assert sum(b['high_confidence_bases'] for b in dataset['bins'])==metrics['high_confidence_bases']
 import pysam
 fasta=pysam.FastaFile(str(CACHE/'chr9.fa'))
 assert filehash(CACHE/'chr9.fa')==m['preparation']['reference_sha256']
 for a in alleles: assert fasta.fetch('9',a['position']-1,a['position']-1+len(a['ref'])).upper()==a['ref']
 if (RESULTS/'clinvar-overlay.json').exists():
  from clinvar import verify_overlay
  verify_overlay()
 # Check a deterministic stratified sample, plus every selected/alternative case row.
 check=[rows[n] for n in range(0,len(rows),max(1,len(rows)//30))]
 for c in dataset['cases']:
  for r in [c['selected'],c['alternative']]:
   if r: check.append(r)
 for r in check:
  matches=[t for t in raw[r['allele_id']].get('transcript_consequences',[]) if t.get('transcript_id','')==r['transcript_id'] and t.get('gene_id','')==r['gene_id']]
  assert any(';'.join(t['consequence_terms'])==r['consequence_terms'] and t.get('hgvsc','')==r['hgvsc'] and t.get('hgvsp','')==r['hgvsp'] for t in matches),r['row_id']
 # Rebuild renderer twice without networking; analytical export hashes stay fixed.
 paths=[ROOT/'report/index.html',ROOT/'report/summary.html']; before={str(p):filehash(p) for p in paths}; render(); render(); assert before=={str(p):filehash(p) for p in paths}
 from html.parser import HTMLParser
 class Links(HTMLParser):
  def __init__(self,p): super().__init__(); self.path=p; self.ids=[]
  def handle_starttag(self,tag,attrs):
   a=dict(attrs)
   if a.get('id'): self.ids.append(a['id'])
   target=a.get('href',a.get('src',''))
   if target and not target.startswith(('https://','#','http://')): assert (self.path.parent/target.split('?',1)[0]).resolve().is_file(),target
 for p in paths:
  parser=Links(p); parser.feed(p.read_text()); assert len(parser.ids)==len(set(parser.ids))
 save(RESULTS/'verification.json',{'status':'passed','alleles':len(alleles),'transcript_rows':len(rows),'stratified_and_case_rows_checked':len(check),'deterministic_report_rebuild':True,'raw_hashes':True,'export_hashes':True,'bin_and_coverage_totals':True,'local_links':True})
 print('PASS: snapshot, exports, case/source rows, stratified sample, bin totals and deterministic offline rendering')

def main():
 p=argparse.ArgumentParser(); p.add_argument('command',choices=['acquire','prepare','annotate','analyze','render','verify']); args=p.parse_args()
 try: globals()[args.command]()
 except Exception as e: print('FAILED:',e,file=sys.stderr); raise
if __name__=='__main__': main()
