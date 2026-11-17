"""Dated, exact normalized ClinVar overlap; no clinical inference."""
import collections,json,sys,time
from pathlib import Path
import pysam,pysam.bcftools
from pipeline import ROOT,CACHE,RESULTS,download,save,read,filehash,now,eligible
URL='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20261108.vcf.gz'
def acquire():
 download(URL+'.tbi',CACHE/'clinvar_20261108.vcf.gz.tbi')
 if not (CACHE/'clinvar-chr9.original.vcf').exists():
  print('Query dated ClinVar chromosome 9',flush=True)
  f=pysam.TabixFile(URL,index=str(CACHE/'clinvar_20261108.vcf.gz.tbi'))
  chrom='9' if '9' in f.contigs else 'chr9'
  p=CACHE/'clinvar-chr9.original.vcf.tmp'
  with p.open('w') as out:
   out.write('\n'.join(f.header)+'\n')
   for line in f.fetch(chrom): out.write(line+'\n')
  p.replace(CACHE/'clinvar-chr9.original.vcf')
 save(CACHE/'clinvar-source.json',{'release':'2026-11-08','acquired_utc':now(),'url':URL,'index_url':URL+'.tbi','extraction':'all chromosome 9 records via index; full remote VCF not downloaded','source_region_sha256':filehash(CACHE/'clinvar-chr9.original.vcf'),'index_sha256':filehash(CACHE/'clinvar_20261108.vcf.gz.tbi')})
def info_map(s):
 return dict(item.split('=',1) if '=' in item else (item,True) for item in s.split(';'))
def overlay():
 original={}; split=[]; qc=collections.Counter(); fasta=pysam.FastaFile(str(CACHE/'chr9.fa'))
 for num,line in enumerate((CACHE/'clinvar-chr9.original.vcf').read_text().splitlines()):
  if line.startswith('#'): continue
  p=line.split('\t'); pos=int(p[1]); ref=p[3]; qc['source_chr9_records']+=1
  for n,alt in enumerate(p[4].split(','),1):
   qc['source_chr9_alt_alleles']+=1
   if not eligible(ref,alt): qc['unsupported_or_over50bp']+=1; continue
   if fasta.fetch('9',pos-1,pos-1+len(ref)).upper()!=ref: qc['reference_mismatch']+=1; continue
   key=f'cv{num}a{n}';original[key]={'clinvar_variation_id':p[2],'original_position':pos,'original_ref':ref,'original_alt':alt,'original_alt_index':n,'info':info_map(p[7]),'original_record':line}
   split.append(f'9\t{pos}\t{key}\t{ref}\t{alt}\t.\t.\t.')
 header='##fileformat=VCFv4.2\n##contig=<ID=9,length=138394717>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
 (CACHE/'clinvar-split.vcf').write_text(header+'\n'.join(split)+'\n')
 norm=pysam.bcftools.norm('--no-version','-f',str(CACHE/'chr9.fa'),'-c','e','-m','-any','-Ov',str(CACHE/'clinvar-split.vcf'))
 (CACHE/'clinvar-normalized.vcf').write_text(norm)
 lookup=collections.defaultdict(list)
 for line in norm.splitlines():
  if line.startswith('#'): continue
  p=line.split('\t');key=f'9:{p[1]}:{p[3]}:{p[4]}'
  lookup[key].append(original[p[2]])
 hits={}; rows=[]
 for a in read(CACHE/'alleles.json'):
  aid=a['allele_id']; entries=lookup.get(aid,[]); hits[aid]=[]
  for entry in entries:
   info=entry['info']; varid=entry['clinvar_variation_id']
   hit={'allele_id':aid,'clinvar_variation_id':varid,'classification_summary':info.get('CLNSIG','not reported'),'review_status':info.get('CLNREVSTAT','not reported'),'conditions':info.get('CLNDN','not reported'),'condition_identifiers':info.get('CLNDISDB','not reported'),'conflicting_classifications':info.get('CLNSIGCONF','not reported'),'somatic_clinical_impact':info.get('SCI','not reported'),'oncogenicity':info.get('ONC','not reported'),'allele_id_clinvar':info.get('ALLELEID','not reported'),'url':'https://www.ncbi.nlm.nih.gov/clinvar/variation/'+varid+'/','release':'2026-11-08','match_method':'exact chromosome/position/REF/ALT after bcftools 1.21 normalization','original_record':entry['original_record']}
   rows.append(hit); hits[aid].append(hit)
 from pipeline import tsv
 tsv(RESULTS/'clinvar-matches.tsv',rows,['allele_id','clinvar_variation_id','classification_summary','review_status','conditions','condition_identifiers','conflicting_classifications','somatic_clinical_impact','oncogenicity','allele_id_clinvar','url','release','match_method','original_record'])
 qc['hg002_alleles']=len(hits);qc['hg002_exact_matched_alleles']=sum(bool(v) for v in hits.values());qc['hg002_unmatched_alleles']=sum(not v for v in hits.values());qc['matched_clinvar_records']=len(rows)
 save(RESULTS/'clinvar-overlay.json',{'schema_version':'1.0','source':read(CACHE/'clinvar-source.json'),'qc':dict(qc),'hits':hits,'scope_note':'All chr9 ClinVar small ACGT alleles <=50 bp were independently reference-checked and normalized. Unmatched means no exact eligible match in this fixed release, not benign or absent from all ClinVar representations. Original INFO strings retain assertion categories, conflicts and condition context; no clinical classification is assigned by this project.'})
 print('ClinVar',dict(qc),flush=True)
if __name__=='__main__':
 if sys.argv[1]=='acquire': acquire()
 elif sys.argv[1]=='overlay': overlay()

def verify_overlay():
 """Re-normalize matched source records and check source assertions independently."""
 import tempfile
 from pipeline import tsv
 overlay=read(RESULTS/'clinvar-overlay.json');source=overlay['source']
 if filehash(CACHE/'clinvar-chr9.original.vcf')!=source['source_region_sha256']: raise ValueError('ClinVar source hash mismatch')
 expected=[];lines=[]
 for aid,entries in overlay['hits'].items():
  for h in entries:
   p=h['original_record'].split('\t'); info=info_map(p[7]); key='match'+str(len(expected))
   assert p[2]==h['clinvar_variation_id']
   assert info.get('CLNSIG','not reported')==h['classification_summary']
   assert info.get('CLNREVSTAT','not reported')==h['review_status']
   assert info.get('CLNDN','not reported')==h['conditions']
   lines.append(f'9\t{p[1]}\t{key}\t{p[3]}\t{p[4]}\t.\t.\t.');expected.append(aid)
 with tempfile.TemporaryDirectory() as tmp:
  p=Path(tmp)/'matches.vcf';p.write_text('##fileformat=VCFv4.2\n##contig=<ID=9,length=138394717>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+'\n'.join(lines)+'\n')
  norm=pysam.bcftools.norm('--no-version','-f',str(CACHE/'chr9.fa'),'-c','e','-Ov',str(p))
  for line in norm.splitlines():
   if line.startswith('#'):continue
   fields=line.split('\t');assert f'9:{fields[1]}:{fields[3]}:{fields[4]}'==expected[int(fields[2][5:])]
 assert sum(bool(v) for v in overlay['hits'].values())==overlay['qc']['hg002_exact_matched_alleles']
 print('PASS: every matched ClinVar source identity, normalization and promoted assertion fields')
