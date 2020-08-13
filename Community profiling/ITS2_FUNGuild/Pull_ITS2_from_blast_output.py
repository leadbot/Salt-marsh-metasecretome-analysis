# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:58:06 2020

@author: Dan
"""

########    Read blast .txt. output with 10 hits per query

import re
from Bio import Entrez
import time
import pickle
import pandas as pd

blast_dict={}

blast_file='BLAST_output_Block_A.txt'
max_OTU='OTU_517'

blast_file='BLAST_output_Block_B.txt'
max_OTU='OTU_999'

blast_file='BLAST_output_Block_C.txt'
max_OTU='OTU_1180'




f=open(blast_file, 'r')
linenum=0
originalpos=0
while True:
    linenum+=1
    f.seek(originalpos)
    line=f.readline().strip('\n')
    originalpos=f.tell()
    if line.startswith('Query #'):
        OTU=line.split(' ', 3)[2]
        blast_dict[OTU]={}
        print(OTU)
        #NB there are 14 lines of interest after the
        headcount=0
        for x in list(range(1,250)):
            line=f.readline().strip('\n')
            #read current position
            currentpos=f.tell()
            line2=f.readline().strip('\n')
            #return to previous position to not offset line reads if uneven numbers (e.g. line1 not containing > or skipping query)
            f.seek(currentpos)
            if line.startswith('Query #'):
                break
            if line.startswith('>'):
                headcount+=1
                print(line, line2, headcount)
                if 'Hit_'+str(headcount) not in blast_dict[OTU]:
                     blast_dict[OTU]['Hit_'+str(headcount)]={}
                blast_dict[OTU]['Hit_'+str(headcount)]['Description']=line
                blast_dict[OTU]['Hit_'+str(headcount)]['Accession']=line2.split(':')[1].split(' ')[1]                       
    if line == '' and OTU == max_OTU:
        print('end',OTU, line)
        break
f.close()



regex_expressions=['internal', 'its2']
ITS2_OTUs=[]
accessions_list=[]
tophit_ITS2_hit_and_accession={}
for otu in blast_dict:
    for hit in blast_dict[otu]:
        for reg in regex_expressions:
             if bool(re.search(reg, blast_dict[otu][hit]['Description'].lower()))==True:
                 if otu not in ITS2_OTUs: 
                      ITS2_OTUs.append(otu)  
                 if otu not in tophit_ITS2_hit_and_accession:
                     tophit_ITS2_hit_and_accession[otu]={}
                     tophit_ITS2_hit_and_accession[otu]['Hit']=hit
                     tophit_ITS2_hit_and_accession[otu]['Description']=blast_dict[otu][hit]['Description']  
                     tophit_ITS2_hit_and_accession[otu]['Accession']=blast_dict[otu][hit]['Accession'] 
                     accessions_list.append(blast_dict[otu][hit]['Accession'])    

##Check all ITS2 blast OTUs are tophit
for x in ITS2_OTUs:
         if bool(re.search('internal', blast_dict[x]['Hit_1']['Description'].lower()))==False and\
                bool(re.search('its2', blast_dict[x]['Hit_1']['Description'].lower()))==False:
                    print(x + 'not hit1')
    
   
c=0              
for otu in blast_dict:
    if otu not in ITS2_OTUs:
        c+=1
        print(otu)
        print(str(c)+'^')
        if 'Hit_1' in blast_dict[otu]:
             print(blast_dict[otu]['Hit_1']['Description'])
             print(blast_dict[otu]['Hit_2']['Description'])
             time.sleep(1)
             
             
###ncbi efetch taxids from accessions####
famcount=0
Entrez.email = 'dl923@york.ac.uk'
records={}
loopcount=0
failed_accessions={}
desired_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
#db="protein"
db="nuccore"
famcount=0

for x in tophit_ITS2_hit_and_accession:
    loopcount+=1
    if loopcount%100==0:
        print('\n#\nParsed ' + str(loopcount) + ' taxids (of ' + str(len(tophit_ITS2_hit_and_accession)) + ') succesfully\n#\n')
    try:
        handle=Entrez.esummary(db=db, id=str(tophit_ITS2_hit_and_accession[x]['Accession']), retmode="xml")
        record=Entrez.read(handle)
        for rec in record:
            tophit_ITS2_hit_and_accession[x]['Taxid']=str(rec['TaxId'])
        print('ORF: ' + x + ' Taxid: ' + str(tophit_ITS2_hit_and_accession[x]['Taxid']) + ' Accession: ' + str(tophit_ITS2_hit_and_accession[x]['Accession']))
        records[x]=record
        time.sleep(0.35)
        taxhandle=Entrez.efetch(db="Taxonomy", id=str(tophit_ITS2_hit_and_accession[x]['Taxid']), retmode='xml')
        taxrecords = Entrez.read(taxhandle)
        tophit_ITS2_hit_and_accession[x]['Lineage']=taxrecords[0]["Lineage"]
        #### pull individual ranks
        tophit_ITS2_hit_and_accession[x]['Taxonomy']={}
        returned_ranks=[]
        dicnum=0
        for lin in taxrecords[0]['LineageEx']:
           returned_ranks.append(tuple((lin['Rank'], lin['ScientificName'], lin['TaxId'], dicnum)))
           dicnum+=1 
        matched_ranks=[]                     
        for tup in returned_ranks:
            if tup[0] in desired_ranks:
                tophit_ITS2_hit_and_accession[x]['Taxonomy'][tup[0]]=str(tup[1])
                matched_ranks.append(tup[0])
                if tup[0]=='family':
                    famcount+=1
        for rnk in desired_ranks:
            if rnk not in matched_ranks:
                tophit_ITS2_hit_and_accession[x]['Taxonomy'][rnk]='NA'            
    except:
        print('Failed Accession: ' + str(tophit_ITS2_hit_and_accession[x]['Accession']))
        failed_accessions['x']='NA'
        continue             
   
       
with open('Blast_ITS2_hits_only_393.pickle', 'wb') as handle:
    pickle.dump(tophit_ITS2_hit_and_accession, handle, protocol=pickle.HIGHEST_PROTOCOL)    
    
with open('Blast_all_unknown_hits_610.pickle', 'wb') as handle:
    pickle.dump(blast_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)        
    
             
###check how many contain fungi
c=0
for x in tophit_ITS2_hit_and_accession:
         if bool(re.search('fungi', tophit_ITS2_hit_and_accession[x]['Lineage'].lower()))==True:
                    print(x + 'fungi')   
                    c+=1

####pull into a single dataframe
dicttodf={}
for otu in tophit_ITS2_hit_and_accession:
    dicttodf[otu]=tophit_ITS2_hit_and_accession[otu]['Taxonomy']

its2_blast_df=pd.DataFrame.from_dict(dicttodf)


fungi_dict={}
off_target_sequence_dict={}
for otu in tophit_ITS2_hit_and_accession:
    if bool(re.search('fungi', tophit_ITS2_hit_and_accession[otu]['Lineage'].lower()))==True:
         fungi_dict[otu]=tophit_ITS2_hit_and_accession[otu]['Taxonomy']
    else:
         off_target_sequence_dict[otu]=tophit_ITS2_hit_and_accession[otu]['Taxonomy']
        
with open('Blast_fungi_hits_only_351.pickle', 'wb') as handle:
    pickle.dump(fungi_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('Blast_offtarget_hits_only_42.pickle', 'wb') as handle:
    pickle.dump(off_target_sequence_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)    