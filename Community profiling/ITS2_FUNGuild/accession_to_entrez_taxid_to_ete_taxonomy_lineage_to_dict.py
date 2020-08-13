# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 10:33:04 2017

@author: dl923
"""
from Bio import Entrez
import time
from ete3 import NCBITaxa
import pickle

ncbi = NCBITaxa()

ncbi.update_taxonomy_database()


##conda install -c etetoolkit ete3
##conda install -c anaconda biopython
##IF THIS SCRIPT FAILS, RUN SPYDER WITH ADMINSTRATOR PRIVLEDGES

save_pre_taxid_dict='yes'
save_post_taxid_dict='yes'

with open('ISDE_blastp_taxonomy_entez_ete3_PRE_taxid.pickle', 'rb') as handle:
     ORF_taxa_dict=pickle.load(handle)
         

taxid_file='blastP\BlastP_result_files\ISDE_cat_blastp_out.csv'

f=open(taxid_file)

ORF_taxa_dict={}
missing_taxids=[]
taxid_list=[]

while True:
    line=f.readline().strip('\n')
    if line=='':
        break
    if line.count(',')==8:
         ORF=line.split(',',1)[0]
         Axsn=line.split(',',3)[2]
         Annotation=line.split(',', 4)[3]
         e_val=line.split(',', 5)[4]
         if ORF not in ORF_taxa_dict:
             ORF_taxa_dict[ORF]={}
        # if ORF in ORF_taxa_dict:
        #     print('duplicate')
         ORF_taxa_dict[ORF]['Accession']=Axsn
         ORF_taxa_dict[ORF]['Tophit']=Annotation 
         ORF_taxa_dict[ORF]['Tophit_eval']=e_val
    elif line.count(',')>=9:
         diff=line.count(',')-8
         ORF=line.split(',',1)[0]
         Axsn=line.split(',',3)[2]
         Annotation=line.split(',', 4+diff)[3+diff]
         e_val=line.split(',', 5+diff)[4+diff]
         if ORF not in ORF_taxa_dict:
             ORF_taxa_dict[ORF]={}
        # if ORF in ORF_taxa_dict:
         #    print('duplicate')
         ORF_taxa_dict[ORF]['Accession']=Axsn
         ORF_taxa_dict[ORF]['Tophit']=Annotation 
         ORF_taxa_dict[ORF]['Tophit_eval']=e_val

###ncbi efetch taxids from accessions####
Entrez.email = 'dl923@york.ac.uk'
records={}
loopcount=0
failed_accessions={}

for x in ORF_taxa_dict:
    loopcount+=1
    if loopcount%100==0:
        print('\n#\nParsed ' + str(loopcount) + ' taxids (of ' + str(len(ORF_taxa_dict)) + ') succesfully\n#\n')
    try:
        handle=Entrez.esummary(db="protein", id=str(ORF_taxa_dict[x]['Accession']), retmode="xml")
        record=Entrez.read(handle)
        for rec in record:
            ORF_taxa_dict[x]['Taxid']=str(rec['TaxId'])
        print('ORF: ' + x + ' Taxid: ' + str(ORF_taxa_dict[x]['Taxid']) + ' Accession: ' + str(ORF_taxa_dict[x]['Accession']))
        records[x]=record
        time.sleep(0.35)  
    except:
        print('Failed Accession: ' + str(ORF_taxa_dict[x]['Accession']))
        failed_accessions['x']='NA'
        ORF_taxa_dict[x]['Taxid']='NA'
        continue
    

if save_pre_taxid_dict.lower()=='yes':
    with open('ISDE_blastp_taxonomy_entez_ete3_PRE_taxid_NEW.pickle', 'wb') as handle:
         pickle.dump(ORF_taxa_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


desired_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
global no_ids
no_ids=[]

def get_desired_ranks_output_dict(accession, taxid, desired_ranks):
    taxdict={}
    try:
         lineage = ncbi.get_lineage(taxid)
         lineage2ranks = ncbi.get_rank(lineage)
         lineagetranslated=ncbi.get_taxid_translator(lineage)
         ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
         taxdict={}
         for x in desired_ranks:
              if x in list(ranks2lineage.keys()):
                    taxdict[x]={}
              else:
                    taxdict[x]='NA'    
         for taxrank, taxid in ranks2lineage.items():
               if taxrank in desired_ranks:
                    taxdict[taxrank]=lineagetranslated[taxid]
               else:
                    taxdict[taxrank]='NA'
         taxdict['lineage']=lineagetranslated
         return taxdict
    except:
         print('Failed to gather taxid for id: ' + str(taxid))
         for x in desired_ranks:
             taxdict[x]='NA'
         taxdict['lineage']='NA'
         no_ids.append(accession)
         return taxdict

loopcount=0
for x in ORF_taxa_dict:
    loopcount+=1
    if loopcount%250==0:
        print('\n#\nParsed ' + str(loopcount) + ' taxids (of ' + str(len(ORF_taxa_dict)) + ') succesfully\n#\n')
    ORF_taxa_dict[x]['Taxonomy']=get_desired_ranks_output_dict(x, ORF_taxa_dict[x]['Taxid'], desired_ranks)
    print(ORF_taxa_dict[x]['Taxonomy']['phylum'])

print('Taxids not identifiable for: ' + str(len(no_ids)) + ' hits')
  
if save_post_taxid_dict.lower()=='yes':          
     with open('ISDE_blastp_taxonomy_entez_ete3_NEW.pickle', 'wb') as handle:
         pickle.dump(ORF_taxa_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

