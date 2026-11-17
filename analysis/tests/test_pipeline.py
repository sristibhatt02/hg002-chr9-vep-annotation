"""Synthetic test cases; no fixtures contribute to the scientific report."""
import copy,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pipeline as p
import pysam,pysam.bcftools
class TestAnnotation(unittest.TestCase):
 def allele(self,pos=100,ref='A',alt='G'):
  return {'allele_id':f'9:{pos}:{ref}:{alt}','position':pos,'ref':ref,'alt':alt}
 def obj(self,x):
  start,end,allele=p.normalized_repr(x)
  return {'input':f"9 {x['position']} {x['allele_id']} {x['ref']} {x['alt']} . . .",'assembly_name':'GRCh38','seq_region_name':'9','strand':1,'start':start,'end':end,'allele_string':allele}
 def test_same_position_different_alleles_and_reorder(self):
  xs=[self.allele(100,'A','G'),self.allele(100,'A','T')]
  self.assertEqual(set(p.reconcile([self.obj(x) for x in xs][::-1],xs)),{x['allele_id'] for x in xs})
 def test_wrong_identity_missing_duplicate_build(self):
  x=self.allele();o=self.obj(x)
  for bad in [[],[o,o],[{**o,'assembly_name':'GRCh37'}],[{**o,'input':'9 99 unknown A G'}]]:
   with self.assertRaises(ValueError):p.reconcile(bad,[x])
 def test_deletion_insertion_minimal(self):
  self.assertEqual(p.normalized_repr(self.allele(100,'TA','T')),(101,101,'A/-'))
  self.assertEqual(p.normalized_repr(self.allele(100,'T','TA')),(101,100,'-/A'))
 def test_interval_edges_full_ref_rule(self):
  bed=[(100,200)]
  self.assertTrue(p.covered(100,101,bed));self.assertTrue(p.covered(198,200,bed))
  self.assertFalse(p.covered(99,101,bed));self.assertFalse(p.covered(199,201,bed))
  self.assertFalse(p.covered(100,102,[(100,101),(102,103)]))
 def test_supported_types(self):
  for ref,alt in [('A','*'),('A','<DEL>'),('N','A'),('A','A'),('A','G'*51)]:self.assertFalse(p.eligible(ref,alt))
  self.assertTrue(p.eligible('AA','A'))
 def test_missing_transcripts_valid(self):
  x=self.allele();self.assertEqual(p.reconcile([self.obj(x)],[x])[x['allele_id']].get('transcript_consequences',[]),[])
 def test_selection_and_severity(self):
  rows=[{'row_id':'a','transcript_id':'T1','mane_select':'','canonical':1},{'row_id':'b','transcript_id':'T2','mane_select':'NM_1','canonical':0}]
  r,tier,ties=p.representative(rows);self.assertEqual(r['row_id'],'b');self.assertEqual(tier,'MANE Select')
  self.assertEqual(p.severity(['intron_variant','missense_variant']),'missense_variant')
  with self.assertRaises(ValueError):p.severity(['unknown_term'])
 def test_bcftools_split_left_align_and_ref_failure(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);fa=root/'ref.fa';fa.write_text('>9\nCAAAAAGTCG\n');pysam.faidx(str(fa))
   vcf=root/'in.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=9,length=10>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n9\t4\ttest\tAA\tA,AAA\t.\tPASS\t.\n')
   result=pysam.bcftools.norm('--no-version','-f',str(fa),'-m','-any','-c','e',str(vcf))
   records=[l.split('\t') for l in result.splitlines() if not l.startswith('#')]
   self.assertEqual(len(records),2);self.assertTrue(all(int(r[1])<4 for r in records));self.assertTrue(all(r[2]=='test' for r in records))
   vcf.write_text(vcf.read_text().replace('AA\tA,AAA','GG\tG'))
   with self.assertRaises(Exception):pysam.bcftools.norm('-f',str(fa),'-c','e',str(vcf))
 def test_real_dataset_identity_and_totals(self):
  self.assertEqual(p.SCOPE['sample'],'HG002');self.assertEqual(p.SCOPE['release'],'GIAB NIST v4.2.1')
  if (p.CACHE/'prepare-qc.json').exists():
   qc=p.read(p.CACHE/'prepare-qc.json');self.assertEqual(qc['ingested_alt_alleles'],qc['reference_matched_before_norm']+qc['outside_high_confidence']);self.assertEqual(len(p.read(p.CACHE/'alleles.json')),qc['deduplicated_eligible_alleles'])
 def test_production_prepare_multiallelic_gt_provenance_and_edge(self):
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);cache=root/'cache';results=root/'results';cache.mkdir();results.mkdir()
   (cache/'chr9.fa').write_text('>9\nCAAAAAGTCG\n');pysam.faidx(str(cache/'chr9.fa'))
   (cache/'benchmark.bed').write_text('chr9\t0\t10\n')
   (cache/'source-manifest.json').write_text('{}')
   header='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'
   rows=['chr9\t7\tmulti\tG\tA,T\t.\tPASS\t.\tGT\t1|2','chr9\t4\tdel\tAA\tA\t.\tPASS\t.\tGT\t0/1','chr9\t10\tedge\tGAA\tG\t.\tPASS\t.\tGT\t0/1','chr9\t8\tnotcalled\tT\tC\t.\tPASS\t.\tGT\t0/0']
   (cache/'region.original.vcf').write_text(header+'\n'.join(rows)+'\n')
   with patch.object(p,'CACHE',cache),patch.object(p,'RESULTS',results),patch.object(p,'SCOPE',{**p.SCOPE,'start_0based':0,'end_0based_exclusive':10}):p.prepare()
   alleles=p.read(cache/'alleles.json');qc=p.read(cache/'prepare-qc.json')
   self.assertEqual(len(alleles),3);self.assertEqual(qc['outside_full_ref_span'],1);self.assertEqual(qc['not_present_in_sample_gt'],1)
   multi=[a for a in alleles if a['position']==7];self.assertEqual({a['alt'] for a in multi},{'A','T'});self.assertEqual({a['origins'][0]['original_alt_index'] for a in multi},{1,2})
   deletion=next(a for a in alleles if a['variant_type']=='indel');self.assertEqual(deletion['position'],1);self.assertEqual(deletion['origins'][0]['original_position'],4)
 def test_clinvar_exact_allele_not_coordinate(self):
  # Identical position with different ALT never matches the normalized key.
  lookup={'9:100:A:G':['example']}
  self.assertEqual(lookup.get('9:100:A:T',[]),[])
  self.assertNotEqual(p.normalized_repr(self.allele(100,'A','G')),p.normalized_repr(self.allele(100,'A','T')))
if __name__=='__main__':unittest.main()
