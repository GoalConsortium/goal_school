#!/bin/env python

import sys
import json
import argparse
import pathlib
import re
import math

def get_args():
    '''Define arguments.'''
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-j', '--jsonfile',
		        help='JSON Report file',
		        type=argparse.FileType('r'),
                        required=True)

    args = parser.parse_args()
    return args

def update_vclass(varclass):
  varclass = re.sub('Missense.+', 'Missense_Mutation',varclass)
  varclass = re.sub('Splice.+', 'Splice_Region',varclass)
  varclass = re.sub('Stop g.+', 'Nonsense_Mutation',varclass)
  varclass = re.sub('Stop l.+', 'Nonstop_Mutation',varclass)
  varclass = re.sub('Frameshift.*', 'Frame_Shift',varclass)
  varclass = re.sub('Inframe del.*', 'In_Frame_Del',varclass)
  varclass = re.sub('Inframe ins.*', 'In_Frame_Ins',varclass)
  varclass = re.sub('disruptive_inframe_ins.*', 'In_Frame_Ins',varclass)
  varclass = re.sub('disruptive_inframe_del.*', 'In_Frame_Del',varclass)
  varclass = re.sub('Start loss.*','Translation_Start_Site',varclass)
  return varclass
  

args = get_args()
jsonref = json.load(args.jsonfile)
patient = jsonref['patient']
orderid = jsonref['order']['accessionId']
reportid = jsonref['report']['reportId']

if 'emrID' in patient:
  mrn = patient['emrID']
elif 'emrId' in patient:
  mrn = patient['emrId']
elif 'emr_id' in patient:
  mrn = patient['emr_id']


results = jsonref['results']
for samp in jsonref['specimens']:
  if samp['sampleCategory'] == 'tumor':
    tumorsite = samp['sampleSite']
  elif samp['sampleCategory'] == 'cfDNA specimen':
    tumorsite = samp['sampleSite']
            
reportedVars = results['somaticPotentiallyActionableMutations'] + results['somaticBiologicallyRelevantVariants'] + results['somaticVariantsOfUnknownSignificance']

pname = orderid + '.patient.txt'
pout = open(pname,'w')
pinfo = [mrn,patient['tempusId'],patient['diagnosis']]
pout.write('\t'.join(pinfo) + '\n')

sname = orderid + '.sample.txt'
sout = open(sname,'w')
sinfo = [patient['tempusId'],orderid,tumorsite]
sout.write('\t'.join(sinfo) + '\n')

mafname = orderid + '.maf'
maf = open(mafname,'w')
fusname = orderid + '.fusion.list'
fusout = open(fusname,'w')

maf.write('\t'.join(['Hugo_Symbol','Entrez_Gene_Id','Variant_Classification','Tumor_Sample_Barcode','HGVSp_Short','t_alt_count','t_ref_count','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2']) + '\n')
for var in reportedVars:
    if var.get('nucleotideAlteration'):
      chrpos = var['nucleotideAlteration']
      geno = chrpos.split(':')
      vclass = var['variantDescription']
      altct = round(float(var['allelicFraction'])*1.5)
      mafinfo = [var['gene'],var['entrezId'],update_vclass(vclass),orderid,var['HGVS.p'],str(altct),str(150)] + geno[0:]
      maf.write('\t'.join(mafinfo) + '\n')
    elif var.get('fusionType'):
        fusout.write(var['gene'] + '\n')
    elif 'variants' in var:
      for genevar in var['variants']:
          chrpos = genevar['nucleotideAlteration']
          geno = chrpos.split(':')
          vclass = genevar['variantDescription']
          altct = round(float(genevar['allelicFraction'])*1.5)
          mafinfo = [var['gene'],var['entrezId'],update_vclass(vclass),orderid,genevar['HGVS.p'],str(altct),str(150)] + geno[0:]
          maf.write('\t'.join(mafinfo) + '\n')

#Frame_Shift_Del
#Frame_Shift_Ins

