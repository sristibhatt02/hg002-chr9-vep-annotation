"""Labelled synthetic fixtures exercise normalization only; never scientific output."""
import copy
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import run

def fixture():
 output=[]
 for i,pos,ref,alt in run.INPUTS:
  start=pos+1 if i=='var1' else pos
  output.append({'input':f'9 {pos} {i} {ref} {alt} . . .','assembly_name':'GRCh38','seq_region_name':'9','strand':1,'start':start,'end':start,'allele_string':'A/-' if i=='var1' else ref+'/'+alt,'transcript_consequences':[{'transcript_id':'SYNTHETIC_A','consequence_terms':['intron_variant']},{'transcript_id':'SYNTHETIC_B','consequence_terms':['splice_region_variant','intron_variant']}]})
 return output
class Normalization(unittest.TestCase):
 def test_reordered_deletion_multiple_missing_optional(self):
  raw=fixture(); a=run.normalize(raw); b=run.normalize(raw[::-1]); self.assertEqual(a,b); self.assertEqual(len(a[0]),10); self.assertEqual(a[0][0]['returned_alleles'],'A/-'); self.assertEqual(a[0][0]['hgvsc'],'')
 def test_missing_input(self):
  with self.assertRaisesRegex(ValueError,'Missing'): run.normalize(fixture()[:-1])
 def test_duplicate(self):
  f=fixture()
  with self.assertRaisesRegex(ValueError,'Duplicated'): run.normalize(f+[f[0]])
 def test_wrong_identity_build_allele(self):
  for key,value in [('input','9 1 var1 TA T'),('assembly_name','GRCh37'),('allele_string','T/A')]:
   f=fixture(); f[0][key]=value
   with self.assertRaises(ValueError): run.normalize(f)
 def test_html_escaping(self): self.assertEqual(run.esc('<script>'),'&lt;script&gt;')
