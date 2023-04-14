# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:30:02 2017

@author: dl923
"""
import xlrd 
import time
from copy import deepcopy
import re 
import pickle
import matplotlib
import matplotlib.patches as patches
import math
import matplotlib.path as mpath
#from matplotlib_venn import venn3
import turtle 
#import canvasvg
#import cairosvg
#import math
#from reportlab.pdfgen import canvas
#from reportlab.lib.utils import ImageReader

#changes

#from pylab import *

#NB dependencies for cairosvg... pip install weasyprint, conda install cairo...
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
#matplotlib.use('TkAgg')


#CAZy homologs (≥ 1e-10) represented only 0.717-0.986 mol% of the total metasecretome in line with previous in vitro reports [alessi 2018); however, as a proportion of identified proteins (2.19% ± 0.27) was greatly reduced in contrast to in vitro studies which report 5-15% CAZy metasecretome representation [xxx] 
###df_master_mean_named[:252][df_master_mean_named.iloc[:252]['Ten']>0]
#num_cazy_per_timepoint=np.array([66, 78, 89,89])
#num_PMO_per_timepoint=np.array([3816, 3503, 3709,3720])
#pct_cazy_total=(num_cazy_per_timepoint/num_PMO_per_timepoint)*100
#pct_mean=pct_cazy_total.mean()
#pct_stdev=pct_cazy_total.std()
###########DATA FILTERING PARAMETERS###############
##default e value is 1e-10!!!
e_value_filter=1e-10
sig_seq_filter=1
phylum_percent_filter=2
class_percent_filter=2
manual_target_list=['AA', 'CE']
manual_entry_of_CBM_only_ORFs='No'
enzyme_class_filter=0.1
extract_top_n_class=100
vector_length=900
###################################################

with open('blastP/BlastP_result_files/Taxonomy_from_blastP/ISDE_blastp_taxonomy_entez_ete3.pickle', 'rb') as fp:
    taxa_master_dict = pickle.load(fp)    
master_nested_dict={}

with open('Old_analysis/old_analysis_annotation_dict.pickle', 'rb') as op:
    old_analysis_annotation_dict = pickle.load(op) 

##########################################################     INPUTS     ##########################################################
with open('blastP/BlastP_result_files/Taxonomy_from_blastP/ISDE_blastp_taxonomy_entez_ete3.pickle', 'rb') as fp:
    taxa_master_dict = pickle.load(fp) 

with open('CBM_activity_table_pickles/cbm_pickle_cellulose.pickle', 'rb') as handle:
    cbm_pickle_cellulose=pickle.load(handle)
    
with open('CBM_activity_table_pickles/cbm_pickle_hemicellulose.pickle', 'rb') as handle:
    cbm_pickle_hemicellulose=pickle.load(handle)
    
with open('CBM_activity_table_pickles/cbm_pickle_chitin.pickle', 'rb') as handle:
    cbm_pickle_chitin=pickle.load(handle)
    
file_master='mascot_100k_comp_searches_with_db_0pt05_dbCAN/DB0pt05_vs_0pt1SIG/ISDE_Master_mascot_COMPOSITE_100k_search_0pt05DB_with_0pt1sig_11268U.xlsx'
file_dbCAN='mascot_100k_comp_searches_with_db_0pt05_dbCAN/DB0pt05_vs_0pt1SIG/dbCAN_OUTPUT_DB0pt05_vs_0pt1SIG.xlsx'
###
read_and_write='no'
all_peptide_matching_ORFs_file='mascot_100k_comp_searches_with_db_0pt05_dbCAN/DB0pt05_vs_0pt1SIG/ISDE_100K_master_0pt05DB_with_0pt1SIGtranslated_total_block_results_11268U.fna'
all_peptide_matching_ORFs_file_nuc='mascot_100k_comp_searches_with_db_0pt05_dbCAN/DB0pt05_vs_0pt1SIG/ISDE_100K_master_0pt05DB_with_0pt1SIGtotal_block_results_11268U.fna'
####
#blastp_tophit='blastP/BlastP_result_files/ISDE_cat_blastp_out.xlsx'
pfam_cat='Pfam/Pfam_results_out/ISDE_PFam_cat_master.xlsx'
####################################################################################################################################

##########################################################     OUTPUTS     #########################################################
excel_out_molarpct='Outputs/output_ORFs_mol_pctage_' + time.strftime('%H%M_%S') + '_.csv'
#output_file='output_dbCANN_annotated_sequences_only_gt300_gt500_ ' + time.strftime('%H%M%S') + '_.fna'
output_file='Outputs/output_translated_dbCANN_annotated_sequences_only_gt300_gt500_ ' + time.strftime('%H%M%S') + '_.fna'
output_file_nuc='Outputs/output_NUCLEOTIDE_dbCANN_annotated_sequences_only_gt300_gt500_ ' + time.strftime('%H%M%S') + '_.fna'
pfam_compatible_out='Outputs/output_translated_dbCANN_annotated_sequences_only_gt300_gt500_pFAM_compatible ' + time.strftime('%H%M%S') + '_.fna'
####################################################################################################################################

workbook_master=xlrd.open_workbook(file_master)
workbook_pmc=xlrd.open_workbook(file_dbCAN)
#Write create a dictionary of sheetnames (key) and indices (create two lists then merge into dictionary using dict(zip(a+b)))
sheetlist=workbook_master.sheet_names()
indexlist=list(range(0,workbook_master.nsheets))
xl_dict=dict(zip(sheetlist, indexlist))
#filter only the sheets containing contigs and empai scores using dict comprehension
filter_xl_dict= {k: v for k, v in xl_dict.items() if k.startswith('ISDE')}
#####Write peptide matching contigs to list then dictionary

peptidewb=workbook_pmc.sheet_by_index(0)
peptide_matching_ORF_list=(list(filter(None, (set(peptidewb.col_values(colx=0, start_rowx=0, end_rowx=None))))))

###########################################################################################################################
#Create dict with dbcann annotations in a list and ORFS as key, to allow easy sorting afterwards                            
annotation_dict={}
sum_activity_dict={}
e_value_filter_keys={}
rowcount=0
for rownum in range(peptidewb.nrows):
    rowcount+=1
    if peptidewb.cell_value(rownum, 1)[:-4].split('_',1)[0] not in sum_activity_dict:
        sum_activity_dict[peptidewb.cell_value(rownum, 1)[:-4].split('_',1)[0]]={}
    if (peptidewb.cell_value(rownum, 0)) not in annotation_dict:
        temp={}
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]={}
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]['e_val']=peptidewb.cell_value(rownum, 2)
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]['start']=peptidewb.cell_value(rownum, 5)
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]['end']=peptidewb.cell_value(rownum, 6)
        if peptidewb.cell_value(rownum, 2) <= e_value_filter:
                  annotation_dict[(peptidewb.cell_value(rownum, 0))]={}
                  annotation_dict[(peptidewb.cell_value(rownum, 0))].update(temp)
        elif peptidewb.cell_value(rownum, 2) > e_value_filter:
                  e_value_filter_keys[(peptidewb.cell_value(rownum, 0))]=temp     
    elif (peptidewb.cell_value(rownum, 0)) in annotation_dict:
        temp={}
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]={}
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]['e_val']=peptidewb.cell_value(rownum, 2)
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]['start']=peptidewb.cell_value(rownum, 5)
        temp[peptidewb.cell_value(rownum, 1)[:-4]+'_'+str(rowcount)]['end']=peptidewb.cell_value(rownum, 6)
        if peptidewb.cell_value(rownum, 2) <= e_value_filter:
                  annotation_dict[(peptidewb.cell_value(rownum, 0))].update(temp)
        elif peptidewb.cell_value(rownum, 2) > e_value_filter:
            if peptidewb.cell_value(rownum, 0) in e_value_filter_keys:
                  e_value_filter_keys[(peptidewb.cell_value(rownum, 0))].update(temp)
            elif peptidewb.cell_value(rownum, 0) not in e_value_filter_keys:
                  e_value_filter_keys[(peptidewb.cell_value(rownum, 0))]=temp                                     

sum_activity_dict_skeleton_copy=deepcopy(sum_activity_dict)

ORFS=0
Domains=0
for x in e_value_filter_keys:
    ORFS+=1
    for y in e_value_filter_keys[x]:
        Domains+=1
print('There are ' + str(ORFS) +' ORFs containing domains below threshold - NB THIS DOES NOT MEAN THE ORF WILL BE FILTERED - it may containg another domain')   
print('There are ' + str(Domains)+ ' domains below the threshold')                                                                                                   
###############################################Pfam annotation dict########################################################
#workbook_pfam=xlrd.open_workbook(pfam_cat)
#peptidepfam=workbook_pfam.sheet_by_index(0)
#pfam_annotation_dict={}
#pfam_master_dict={}
#for rownum in range(peptidepfam.nrows):
#    rowcount+=1
#    ORF=peptidepfam.cell_value(rownum, 0)[1:].split(' ',1)[0]
#    if ORF in annotation_dict and (peptidepfam.cell_value(rownum, 0)) not in pfam_annotation_dict:
#        temp={}
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]={}
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]['e_val']=peptidepfam.cell_value(rownum, 12)
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]['start']=peptidepfam.cell_value(rownum, 3)
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]['end']=peptidepfam.cell_value(rownum, 4)
#        pfam_annotation_dict[ORF]={}
#        pfam_annotation_dict[ORF].update(temp)
#    elif (ORF) in annotation_dict and (peptidepfam.cell_value(rownum, 0)) in pfam_annotation_dict:
#        temp={}
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]={}
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]['e_val']=peptidepfam.cell_value(rownum, 12)
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]['start']=peptidepfam.cell_value(rownum, 3)
#        temp[peptidepfam.cell_value(rownum, 6)+'_'+str(rowcount)]['end']=peptidepfam.cell_value(rownum, 4)
#        pfam_annotation_dict[ORF].update(temp)  

workbook_pfam=xlrd.open_workbook(pfam_cat)
peptidepfam=workbook_pfam.sheet_by_index(0)
pfam_annotation_dict={}
pfam_master_dict={}
rowcount=0
for rownum in range(peptidepfam.nrows):
    rowcount+=1
    ORF=peptidepfam.cell_value(rownum, 0)[1:].split(' ',1)[0]
    temp={}
    temp[str(peptidepfam.cell_value(rownum, 6))+'_'+str(rowcount)]={}
    temp[str(peptidepfam.cell_value(rownum, 6))+'_'+str(rowcount)]['e_val']=peptidepfam.cell_value(rownum, 12)
    temp[str(peptidepfam.cell_value(rownum, 6))+'_'+str(rowcount)]['start']=peptidepfam.cell_value(rownum, 3)
    temp[str(peptidepfam.cell_value(rownum, 6))+'_'+str(rowcount)]['end']=peptidepfam.cell_value(rownum, 4)
    if ORF not in pfam_master_dict:
        pfam_master_dict[ORF]={}
        pfam_master_dict[ORF].update(temp)
    if ORF in pfam_master_dict:
        pfam_master_dict[ORF].update(temp)          
    if ORF in annotation_dict and (peptidepfam.cell_value(rownum, 0)) not in pfam_annotation_dict:
        pfam_annotation_dict[ORF]={}
        pfam_annotation_dict[ORF].update(temp)
    elif (ORF) in annotation_dict and (peptidepfam.cell_value(rownum, 0)) in pfam_annotation_dict:
        pfam_annotation_dict[ORF].update(temp)  

pfam_DUF_dict={} 
pfam_list_search=['DUF','peroxidase','glyco_hydro']
for rownum in range(peptidepfam.nrows):
    for x in pfam_list_search:
         if bool(re.search(x,str(peptidepfam.cell_value(rownum, 6))))==True:
             pfam_DUF_dict[peptidepfam.cell_value(rownum, 0)]=peptidepfam.cell_value(rownum, 6)                        
#Use the discarded e value filtered dictionary to map low e value domains onto ORFs containing high e value domains
# E.g. an ORF with 2 domains, GH6 e-5 and CBM10 e-100, with GH6 would be filtered, but the confidence of both is high
#Therefore re-map the GH6 (which would be filtered above), back onto the ORF in annotation dict
supplanted=0
for orfname in e_value_filter_keys:
    if orfname in annotation_dict:
        supplanted+=1
        annotation_dict[orfname].update(e_value_filter_keys[orfname])
print('There have been ' + str(supplanted) + ' low confidence domains supplanted back into ORFs which have a high confidence domains')

for pmc in peptide_matching_ORF_list: 
    if str(pmc) == 'Query': peptide_matching_ORF_list.remove('Query')
print("There are " + str(len(peptide_matching_ORF_list)) + " ORFs with CAZy annotated domains present")

##################################################################################################
def create_dict(input):
    # fw=open(filename[0:-4]+'_fungal_seq_only.fna', 'a')
     #dictname=(str((k.split('_', 1))[1:])).strip("'[]'")#[:-5])+'_dict'
     samplename=str(k.split('_',1)[1])
     PMCstat_dict={}
     for rownum in range(input.nrows)[25:]:   
     #my_dict[key]=value when assigning to dictionary
     #Check if the ORF key is already in dict, if not create an empty dict
     #Check if sample is a key, within the ORF dict, if not, create an empty dict
     #Then append using update function the empai, sigseq and molarpct as the value to the empty dictionary keys above
     #Delete the empty key
                    if wb.cell_value(rownum, 3) not in master_nested_dict:
                        master_nested_dict[input.cell_value(rownum, 3)]={}
                    if samplename not in master_nested_dict[input.cell_value(rownum, 3)]:
                        master_nested_dict[input.cell_value(rownum, 3)][samplename]={}
                        PMCstat_dict['sigseq'] = (input.cell_value(rownum, 9))
                        PMCstat_dict['empai'] = (input.cell_value(rownum, 10))
                        PMCstat_dict['molarpct'] = (input.cell_value(rownum, 11))
                        master_nested_dict[input.cell_value(rownum, 3)][samplename].update(PMCstat_dict)
     #This will result in a dict like{ ORFname { W1_c201 {sigseq:1,empai:1,molarpct:1} } }
     if '' in master_nested_dict:
          del master_nested_dict['']
     print('Deleting empty key from ' + samplename)
     return
##################################################################################################
def annotate_dictionary_with_activity_filter_e_value(input_dict, annotation_dict, output_dict):
    for key in input_dict:
       ORF=str(key)
       keylist=[]
       if key in annotation_dict:
           keylist = [x.split('_', 1)[0] for x in list(annotation_dict[key].keys())]
           dictnam=', '.join(list(reversed(sorted(keylist))))+' (' + str(ORF)+ ')'
           output_dict[dictnam]=input_dict[key]
    return (output_dict)
##################################################################################################
def annotate_dictionary_with_activity_WITHOUT_filter(input_dict, annotation_dict, output_dict):
    for key in input_dict:
       ORF=str(key)
       keylist=[]
       if key in annotation_dict:
           keylist = [x.split('_', 1)[0] for x in list(annotation_dict[key].keys())]
           dictnam=', '.join(list(reversed(sorted(keylist))))+' (' + str(ORF)+ ')'
           output_dict[dictnam]=input_dict[key]
       else:
           output_dict[key]=input_dict[key]
    return (output_dict)
###################################################################################################
print("Reading spreadsheet " + str(file_master) + "\n")
print("Opening and extracting protein matching ORFs from spreadsheet")
for k in filter_xl_dict:
    wb=workbook_master.sheet_by_index(filter_xl_dict[k])
    create_dict(wb)    
###################CROSS-REFERENCING#############################################################
crossref_nested_dict={}
#NB/ REMEMBER the emboss getorf used to translate the sequences adds an additional _1 to the end of the header name, so use [:-2]
print("Consctructed new dictionary via cross referencing CAZy annotated peptide matching ORFs")

###################This script will cross reference and pull ORFs into a crossref dictionary from dbCAN output######################
for key, value in master_nested_dict.items():
    temp_dict={}
    for listitem in peptide_matching_ORF_list:
            if key == listitem:
                crossref_nested_dict[key]=value
print("Crossref complete")
######################################
#####Insert empty keys into the DEEPCOPY of master nested dict, so they are all the same length with 0 values######
copy_crossref_nested_dict=deepcopy(crossref_nested_dict)
for key in copy_crossref_nested_dict:
     sample_dict={}
     temp_dict={}
     temp_dict['sigseq'] = 0
     temp_dict['empai'] = 0
     temp_dict['molarpct'] = 0
     for samplename in filter_xl_dict:
          if samplename.split('_',1)[1] not in copy_crossref_nested_dict[key]:
                sample_dict[(samplename.split('_',1)[1])]=temp_dict
     copy_crossref_nested_dict[key].update(sample_dict)
################Build master_nested_dict######################
for key in master_nested_dict:
     for samplename in filter_xl_dict:
          if samplename.split('_',1)[1] in master_nested_dict[key]:
               master_nested_dict[key][(samplename.split('_',1)[1])].pop('sigseq')
               master_nested_dict[key][(samplename.split('_',1)[1])].pop('empai')
               placeholder=master_nested_dict[key][(samplename.split('_',1)[1])]['molarpct']
               master_nested_dict[key][(samplename.split('_',1)[1])]=placeholder
          elif samplename.split('_',1)[1] not in master_nested_dict[key]:
                master_nested_dict[key][(samplename.split('_',1)[1])]=0
########################################ALL DATA########################
df_master=pd.DataFrame.from_dict(master_nested_dict, orient='index')
df_master_mean=pd.DataFrame(index=df_master.index)
df_master_mean['One']=(df_master['W1_C201'] + df_master['W1_3_C429'] + df_master['W1_4_C429'])/3
df_master_mean['Three']=(df_master['W3_C201'] + df_master['W3_2_C429'] + df_master['W3_1_C429'])/3
df_master_mean['Five']=(df_master['W5_C201'] + df_master['W5_25_C429'] + df_master['W5_100_C429'])/3
df_master_mean['Ten']=(df_master['W10_C201'] + df_master['W10_100_i_C429'] + df_master['W10_100_ii_C429'])/3
              
df_master_mean['sum']=(df_master_mean['One'] + df_master_mean ['Three'] + df_master_mean ['Five'] + df_master_mean ['Ten'])             
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_master_mean.values, df_master_mean.values.T))
x = sns.clustermap(df_master_mean, row_linkage=row_linkage, col_linkage=col_linkage, figsize=(9,28), standard_scale=0, cmap='Blues', yticklabels=False)

df_master_mean['sum']=(df_master_mean['One'] + df_master_mean ['Three'] + df_master_mean ['Five'] + df_master_mean ['Ten']) 
df_master_mean['rank']=df_master_mean['sum'].rank(ascending=False)
df_master_mean[df_master_mean['rank']<200]

######WRITE_all_dbCANN_annotated_seqs_to_file###########


if read_and_write == 'yes':
     seqs={}
     for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
         seqs[seq.id]=seq
     target_ORFs=[]
     ORFcount=0
     for contigname in peptide_matching_ORF_list:
         ORFname=contigname
         if ORFname in seqs:
             target_ORFs.append(seqs[ORFname])
             ORFcount+=1
     print(target_ORFs)
     SeqIO.write(target_ORFs, output_file, 'fasta')
     print("Written " + str(ORFcount) + " sequences in total\nWriting complete\n") 



   
######WRITE_all_dbCANN_annotated_seqs_to_file_NUCLEOTIDE###########
if read_and_write == 'yes':
     seqs_nuc={}
     for seq in SeqIO.parse(all_peptide_matching_ORFs_file_nuc, 'fasta'):
         seqs_nuc[seq.id]=seq
     target_ORFs=[]
     ORFcount=0
     for contigname in peptide_matching_ORF_list:
         ORFname=contigname
         if ORFname in seqs_nuc:
             target_ORFs.append(seqs_nuc[ORFname])
             ORFcount+=1
     print(target_ORFs)
     SeqIO.write(target_ORFs, output_file_nuc, 'fasta')
     print("Written " + str(ORFcount) + " sequences in total\nWriting complete\n")




  
###########################################
#####TEST####
#for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
#    seq.description=''
#    print(seq.format('fasta'))
####TEST#####
#####Anotation_compatible_fna#########################
seqs={}
target_ORFs=[]
ORFcount=0
for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
    seq.description=''
    seqs[seq.id]=seq
    target_ORFs.append(seq)
    ORFcount+=1
print(str(len(target_ORFs)))
print(str(ORFcount))
SeqIO.write(target_ORFs, pfam_compatible_out, 'fasta')
#SeqIO.write(target_ORFs, 'composite_translated_matched_contigs_12sample_Pfam_compatible.fna', 'fasta')
print("Written " + str(ORFcount) + " sequences in total\nWriting complete\n")
#######################################################

###############Write_excel################
#Create samplename list
#As dictionaries are unordered, we need to loop through the list to get the order correct  
samplenamelist=[]
for x in filter_xl_dict:
     if x.startswith('Composite'):
         samplenamelist.append(x.split('_',1)[1])

print('Annotating ORF IDs with enzyme classes')

taxa_master_dict_named={}
taxa_master_dict_named=annotate_dictionary_with_activity_filter_e_value(taxa_master_dict, annotation_dict, taxa_master_dict_named)

copy_crossref_nested_dict_named={}
copy_crossref_nested_dict_named=annotate_dictionary_with_activity_filter_e_value(copy_crossref_nested_dict, annotation_dict, copy_crossref_nested_dict_named)

e_value_filter_keys_named={}
e_value_filter_keys_named=annotate_dictionary_with_activity_filter_e_value(e_value_filter_keys, e_value_filter_keys, e_value_filter_keys_named)

###################################ANALYSIS BEFORE FILTERING ETC#############################
master_nested_dict_named={}
master_nested_dict_named=annotate_dictionary_with_activity_WITHOUT_filter(master_nested_dict, annotation_dict, master_nested_dict_named)

df_master_named=pd.DataFrame.from_dict(master_nested_dict_named, orient='index')
df_master_mean_named=pd.DataFrame(index=df_master_named.index)
df_master_mean_named['One']=(df_master_named['W1_C201'] + df_master_named['W1_3_C429'] + df_master_named['W1_4_C429'])/3
df_master_mean_named['Three']=(df_master_named['W3_C201'] + df_master_named['W3_2_C429'] + df_master_named['W3_1_C429'])/3
df_master_mean_named['Five']=(df_master_named['W5_C201'] + df_master_named['W5_25_C429'] + df_master_named['W5_100_C429'])/3
df_master_mean_named['Ten']=(df_master_named['W10_C201'] + df_master_named['W10_100_i_C429'] + df_master_named['W10_100_ii_C429'])/3

#df_master_mean_named['sum']=(df_master_mean_named['One'] + df_master_mean_named ['Three'] + df_master_mean_named ['Five'] + df_master_mean_named ['Ten'])             
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_master_mean_named.values, df_master_mean_named.values.T))
x = sns.clustermap(df_master_mean_named, row_linkage=row_linkage, col_linkage=col_linkage, figsize=(9,28), standard_scale=0, cmap='Blues', yticklabels=False)

############################PIE CHART WITH UNKNOWN MOLAR PERCENTAGE PROTEINS###############################
GHw1=0
GHw3=0
GHw5=0
GHw10=0
AAw1=0
AAw3=0
AAw5=0
AAw10=0
CEw1=0
CEw3=0
CEw5=0
CEw10=0
PLw1=0
PLw3=0
PLw5=0
PLw10=0
dockw1=0
dockw3=0
dockw5=0
dockw10=0
SLHw1=0
SLHw3=0
SLHw5=0
SLHw10=0
ukwnw1=0
ukwnw3=0
ukwnw5=0
ukwnw10=0
cohw1=0
cohw3=0
cohw5=0
cohw10=0
for x in df_master_mean_named.index:
     if bool(re.search('GH', x))==True:
        GHw1+=df_master_mean_named['One'][x]
        GHw3+=df_master_mean_named['Three'][x]
        GHw5+=df_master_mean_named['Five'][x]
        GHw10+=df_master_mean_named['Ten'][x]
     elif bool(re.search('AA', x))==True:
        AAw1+=df_master_mean_named['One'][x]
        AAw3+=df_master_mean_named['Three'][x]
        AAw5+=df_master_mean_named['Five'][x]
        AAw10+=df_master_mean_named['Ten'][x]
     elif bool(re.search('CE', x))==True:
        CEw1+=df_master_mean_named['One'][x]
        CEw3+=df_master_mean_named['Three'][x]
        CEw5+=df_master_mean_named['Five'][x]
        CEw10+=df_master_mean_named['Ten'][x]
     elif bool(re.search('PL', x))==True:
       PLw1+=df_master_mean_named['One'][x]
       PLw3+=df_master_mean_named['Three'][x]
       PLw5+=df_master_mean_named['Five'][x]
       PLw10+=df_master_mean_named['Ten'][x]
     elif bool(re.search('dockerin', x))==True:
       dockw1+=df_master_mean_named['One'][x]
       dockw3+=df_master_mean_named['Three'][x]
       dockw5+=df_master_mean_named['Five'][x]
       dockw10+=df_master_mean_named['Ten'][x]
     elif bool(re.search('SLH', x))==True:
       SLHw1+=df_master_mean_named['One'][x]
       SLHw3+=df_master_mean_named['Three'][x]
       SLHw5+=df_master_mean_named['Five'][x]
       SLHw10+=df_master_mean_named['Ten'][x]
     elif bool(re.search('cohesin', x))==True:
       cohw1+=df_master_mean_named['One'][x]
       cohw3+=df_master_mean_named['Three'][x]
       cohw5+=df_master_mean_named['Five'][x]
       cohw10+=df_master_mean_named['Ten'][x]
     else:
       ukwnw1+=df_master_mean_named['One'][x]  
       ukwnw3+=df_master_mean_named['Three'][x]
       ukwnw5+=df_master_mean_named['Five'][x]
       ukwnw10+=df_master_mean_named['Ten'][x]   
    
pie_chart_array=np.array([['Index', 'One', 'Three', 'Five', 'Ten'], ['GH', GHw1, GHw3, GHw5, GHw10], ['AA', AAw1, AAw3, AAw5, AAw10], ['CE', CEw1, CEw3, CEw5, CEw10], ['PL', PLw1, PLw3, PLw5, PLw10], ['Dockerin', dockw1, dockw3, dockw5, dockw10], ['SLH', SLHw1, SLHw3, SLHw5, SLHw10], ['cohensin', cohw1, cohw3, cohw5, cohw10], ['Non-CAZy', ukwnw1, ukwnw3, ukwnw5, ukwnw10]])       
dfpiechart=pd.DataFrame(data=pie_chart_array[1:,1:], index=pie_chart_array[1:, 0], columns=pie_chart_array[0,1:])
#change values to floats
for x in dfpiechart:
    dfpiechart[x]=dfpiechart[x].astype(float)
###calculate percentages
for x in dfpiechart:
    newstr=x+'_pct'
    dfpiechart[newstr]=dfpiechart[x]/dfpiechart[x].sum()*100
fontsize=14
plt.rc('font', size=fontsize)
fig_pie, ((ax_pie1, ax_pie3), (ax_pie5, ax_pie10)) = plt.subplots(2, 2, figsize=(12,12))
labels=list(dfpiechart.index)
piecolors=sns.xkcd_palette(['seafoam', 'bubble gum pink', 'pastel yellow', 'sky', 'orchid', 'pinky red', 'pale lilac', 'light grey'])
pie_handle_list=[]
colct=0
ukwnpctw1=dfpiechart['One']['Non-CAZy']
ukwnpctw3=dfpiechart['Three']['Non-CAZy']
ukwnpctw5=dfpiechart['Five']['Non-CAZy']
ukwnpctw10=dfpiechart['Ten']['Non-CAZy']
for x in labels:
    pie_handle_list.append(patches.Patch(color=piecolors[colct], label=x))
    colct+=1   
ax_pie1.pie(dfpiechart['One_pct'], colors=piecolors, startangle=90)
ax_pie1.axis('equal')
plt.text(-2.9,2.5, str(ukwnpctw1)[:6]+'%', fontsize=20)
plt.text(-.2,2.5, str(ukwnpctw3)[:6]+'%', fontsize=20)
plt.text(-2.9,-.1, str(ukwnpctw5)[:6]+'%', fontsize=20)
plt.text(-.2,-.1, str(ukwnpctw10)[:6]+'%', fontsize=20)
ax_pie1.set_title('Week one', size=20)
ax_pie3.pie(dfpiechart['Three_pct'], colors=piecolors, startangle=90)
ax_pie3.axis('equal')
ax_pie3.set_title('Week three', size=20)
ax_pie5.pie(dfpiechart['Five_pct'], colors=piecolors, startangle=90)
ax_pie5.axis('equal')
ax_pie5.set_title('Week five', size=20)
ax_pie10.pie(dfpiechart['Ten_pct'], colors=piecolors, startangle=90)
ax_pie10.axis('equal')
ax_pie10.set_title('Week ten', size=20)
legend=plt.legend(handles=pie_handle_list, loc=(-1.15,-.3), title='Enzyme class', ncol=4, fontsize=20)
plt.setp(legend.get_title(),fontsize=20)
plt.text(-3.75,3.75,'a)', fontsize=28)
plt.show()                    
fig_pie.savefig('OUT_pie/CAZy_vs_total_proteome.png')   
plt.clf()
plt.close(fig_pie)               
#####################FILTER_BY_NUMBER_OF_SIGNIFICANT_SEQUENCES_AND_BUILD_DICTIONARY_FOR_DATAFRAME(MOLARPCT)#######################
crossref_molarpct={}
for k in copy_crossref_nested_dict_named:
    if k not in crossref_molarpct:
         crossref_molarpct[k]={}
    for k2 in copy_crossref_nested_dict_named[k]:
          if k2 not in crossref_molarpct[k]:
               crossref_molarpct[k][k2]={}
          if copy_crossref_nested_dict_named[k][k2]['sigseq']>=sig_seq_filter:
              crossref_molarpct[k][k2]=copy_crossref_nested_dict_named[k][k2]['molarpct']
          else:
              crossref_molarpct[k][k2]=0
        
print("Filtered proteins of less than " + str(sig_seq_filter) + " significant peptide matches")
#########################################################################################   
keys_to_filter_list=['GT']
#########################################################################################
def filter_by_activity(input_dict, list_to_filter):
     del_list=[]
     countpop=0
     for x in input_dict:
         for filt in list_to_filter:
             if bool(re.search(filt, x))==True:
                 del_list.append(x)
     for x in del_list:
         input_dict.pop(x)
         countpop+=1
     print('Popped ' + str(countpop) + ' keys by filtering items in list: ' + str(list_to_filter))
     return input_dict            
#########################################################################################################################################
#Filter with the definition above
crossref_molarpct=filter_by_activity(crossref_molarpct, keys_to_filter_list)
print('\n#\n#\nThere are ' +str(len(crossref_molarpct)) + ' ORFs taken forward for analysis in total post filtering\n#\n#')
###filtering complete##########################filtering complete##########################filtering complete##########################filtering complete#################################

#################################CONSTRUCT LISTS OF XXX ONLY CONTAINING ORFs############################
CBM_only_ORF_list=[]
for x in crossref_molarpct:
    if (re.search(r'CBM', x)) and not (re.search(r'GH',x)) and not (re.search(r'CE',x)) and not (re.search(r'AA',x)) and not (re.search(r'GT',x)) and not (re.search(r'PL',x)):
        CBM_only_ORF_list.append(x) 
        
AA_only_ORF_list=[]
for x in crossref_molarpct:
    if (re.search(r'AA', x)) and not (re.search(r'GH',x)) and not (re.search(r'CE',x)) and not (re.search(r'GT',x)) and not (re.search(r'PL',x)) and not (re.search(r'CBM',x)):
        AA_only_ORF_list.append(x)
        
CE_only_ORF_list=[]
for x in crossref_molarpct:
    if (re.search(r'CE', x)) and not (re.search(r'GH',x)) and not (re.search(r'AA',x)) and not (re.search(r'GT',x)) and not (re.search(r'PL',x)) and not (re.search(r'CBM',x)):
        CE_only_ORF_list.append(x)
        
Cellulosome_associated_only_ORF_list=[]        
for x in crossref_molarpct:
    if not (re.search(r'CE', x)) and not (re.search(r'GH',x)) and not (re.search(r'AA',x)) and not (re.search(r'GT',x)) and not (re.search(r'PL',x)) and not (re.search(r'CBM',x)):
        Cellulosome_associated_only_ORF_list.append(x)
        
PL_only_ORF_list=[]        
for x in crossref_molarpct:
    if not (re.search(r'PL', x)) and not (re.search(r'GH',x)) and not (re.search(r'AA',x)) and not (re.search(r'GT',x)) and not (re.search(r'CE',x)) and not (re.search(r'CBM',x)):
        PL_only_ORF_list.append(x) 

#######################################Top XXX most abudant ORFS  - bar plot######################################        
        
CAZy_ORF_mean_dict=deepcopy(crossref_molarpct)
df_CAZy_ORF=pd.DataFrame.from_dict(CAZy_ORF_mean_dict, orient='index')
df_mean_CAZy_ORF=pd.DataFrame(index=df_CAZy_ORF.index)
df_mean_CAZy_ORF['One']=(df_CAZy_ORF['W1_C201'] + df_CAZy_ORF['W1_3_C429'] + df_CAZy_ORF['W1_4_C429'])/3
df_mean_CAZy_ORF['Three']=(df_CAZy_ORF['W3_C201'] + df_CAZy_ORF['W3_2_C429'] + df_CAZy_ORF['W3_1_C429'])/3
df_mean_CAZy_ORF['Five']=(df_CAZy_ORF['W5_C201'] + df_CAZy_ORF['W5_25_C429'] + df_CAZy_ORF['W5_100_C429'])/3
df_mean_CAZy_ORF['Ten']=(df_CAZy_ORF['W10_C201'] + df_CAZy_ORF['W10_100_i_C429'] + df_CAZy_ORF['W10_100_ii_C429'])/3
df_mean_CAZy_ORF['sum']=(df_mean_CAZy_ORF['One'] + df_mean_CAZy_ORF['Three'] + df_mean_CAZy_ORF['Five'] + df_mean_CAZy_ORF['Ten'])
df_mean_CAZy_ORF['sum_rank']=df_mean_CAZy_ORF['sum'].rank(ascending=False)
top_x_ORFs=20
#ORF_top_20_overall=df_mean_CAZy_ORF[df_mean_CAZy_ORF['sum_rank']<=top_x_ORFs].sort(columns='sum_rank').index.tolist()        
ORF_top_20_overall=df_mean_CAZy_ORF[df_mean_CAZy_ORF['sum_rank']<=top_x_ORFs].sort_values(by=['sum_rank']).index.tolist()

###stats stuff
nums=df_CAZy_ORF.sum(axis=0)
cola=df_CAZy_ORF.columns
w1=[nums[8], nums[9], nums[10]]
w3=[nums[0], nums[1], nums[11]]
w5=[nums[2], nums[3], nums[4]]
w10=[nums[5], nums[6], nums[7]]


print('Plotting mean molar percentage for top ' + str(top_x_ORFs) + ' ORFs')

f, axarr = plt.subplots(4, 5, figsize=(15, 12))
sns.set(style="white")
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)
colour_barplot=sns.xkcd_palette(['clear blue'])
#axarr[row, col]
print('Top 20 ORFs')
for ORF in list(range(len(ORF_top_20_overall))):
     xvals=list(df_mean_CAZy_ORF.loc[ORF_top_20_overall[ORF]][:4])
     width=0.65
     ind=[0, 1, 2, 3]
     yerrx=list(df_mean_CAZy_ORF.loc[ORF_top_20_overall[ORF]][4:8])
     chartname=ORF_top_20_overall[ORF].split('(', 1)[0]+'\n'+ORF_top_20_overall[ORF].split('(',1)[1].strip(')')
     if ORF<=4:
          axarr[0, ORF].bar(ind, xvals, width, color=colour_barplot)# yerr=yerrx)
          axarr[0, ORF].set_title(chartname, size=10.5)
          axarr[0, ORF].set_ylim([0, 0.4])
          axarr[0, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
     elif ORF<=9:
          axarr[1, ORF-5].bar(ind, xvals, width, color=colour_barplot)# yerr=yerrx)
          axarr[1, ORF-5].set_title(chartname, size=10.5)
          axarr[1, ORF-5].set_ylim([0, 0.4])
          axarr[1, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
     elif ORF<=14:
          axarr[2, ORF-10].bar(ind, xvals, width, color=colour_barplot)# yerr=yerrx)
          axarr[2, ORF-10].set_title(chartname, size=10.5)
          axarr[2, ORF-10].set_ylim([0, 0.4])
          axarr[2, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
     elif ORF<=19:
          axarr[3, ORF-15].bar(ind, xvals, width, color=colour_barplot)# yerr=yerrx)
          axarr[3, ORF-15].set_xlabel('Week', fontsize=16)
          axarr[3, ORF-15].set_ylim([0, 0.4])
          axarr[3, ORF-15].set_title(chartname, size=10.5)
          axarr[3, ORF-15].set_xticks([0, 1, 2, 3])          
          axarr[3, ORF-15].set_xticklabels(['One', 'Three', 'Five', 'Ten'], size=16)
          axarr[3, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
plt.setp([a.get_xticklabels() for a in axarr[2, :]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 3]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 4]], visible=False)
sns.despine(top=True, right=True)
f.savefig('OUT_grids/top_20_ORFs.png')
plt.show()
plt.clf()
plt.close(f)  
                
######################################CREATE DF FOR AA, CE, PL, Cellulosome associated only#########################################
df_CBM_only_ORFs=deepcopy(df_mean_CAZy_ORF)
df_CBM_only_ORFs=df_mean_CAZy_ORF.loc[CBM_only_ORF_list]
df_CBM_only_ORFs['rank_rank']=df_CBM_only_ORFs['sum_rank'].rank(ascending=True)

#################################################Update crossref_molarpct dictionary with taxonomy################################## 
for k in crossref_molarpct:
          if k in taxa_master_dict_named:
                    crossref_molarpct[k]['Taxonomy']=taxa_master_dict_named[k]['Taxonomy']
  
df=pd.DataFrame.from_dict(crossref_molarpct, orient='index')

#create empty dataframe with indexes and name the index
df_mean=pd.DataFrame(index=df.index)
df_mean.index.names=['Normalised molar percentage']
#create new empty dataframe, then take the mean of the replicates and put them into that dataframe
df_mean['One']=(df['W1_C201'] + df['W1_3_C429'] + df['W1_4_C429'])/3
df_mean['Three']=(df['W3_C201'] + df['W3_2_C429'] + df['W3_1_C429'])/3
df_mean['Five']=(df['W5_C201'] + df['W5_25_C429'] + df['W5_100_C429'])/3
df_mean['Ten']=(df['W10_C201'] + df['W10_100_i_C429'] + df['W10_100_ii_C429'])/3
       
###Class grouping pie charts#####
#                               #
#################################
GHw1=0
GHw3=0
GHw5=0
GHw10=0
AAw1=0
AAw3=0
AAw5=0
AAw10=0
CEw1=0
CEw3=0
CEw5=0
CEw10=0
PLw1=0
PLw3=0
PLw5=0
PLw10=0
dockw1=0
dockw3=0
dockw5=0
dockw10=0
SLHw1=0
SLHw3=0
SLHw5=0
SLHw10=0
cohw1=0
cohw3=0
cohw5=0
cohw10=0
for x in df_mean.index:
     if bool(re.search('GH', x))==True:
        GHw1+=df_mean['One'][x]
        GHw3+=df_mean['Three'][x]
        GHw5+=df_mean['Five'][x]
        GHw10+=df_mean['Ten'][x]
     if bool(re.search('AA', x))==True:
        AAw1+=df_mean['One'][x]
        AAw3+=df_mean['Three'][x]
        AAw5+=df_mean['Five'][x]
        AAw10+=df_mean['Ten'][x]
     if bool(re.search('CE', x))==True:
        CEw1+=df_mean['One'][x]
        CEw3+=df_mean['Three'][x]
        CEw5+=df_mean['Five'][x]
        CEw10+=df_mean['Ten'][x]
     if bool(re.search('PL', x))==True:
       PLw1+=df_mean['One'][x]
       PLw3+=df_mean['Three'][x]
       PLw5+=df_mean['Five'][x]
       PLw10+=df_mean['Ten'][x]
     if bool(re.search('dockerin', x))==True:
       dockw1+=df_mean['One'][x]
       dockw3+=df_mean['Three'][x]
       dockw5+=df_mean['Five'][x]
       dockw10+=df_mean['Ten'][x]
     if bool(re.search('SLH', x))==True:
       SLHw1+=df_mean['One'][x]
       SLHw3+=df_mean['Three'][x]
       SLHw5+=df_mean['Five'][x]
       SLHw10+=df_mean['Ten'][x]
     if bool(re.search('cohesin', x))==True:
       cohw1+=df_mean['One'][x]
       cohw3+=df_mean['Three'][x]
       cohw5+=df_mean['Five'][x]
       cohw10+=df_mean['Ten'][x]
      
pie_chart_array=np.array([['Index', 'One', 'Three', 'Five', 'Ten'], ['GH', GHw1, GHw3, GHw5, GHw10], ['AA', AAw1, AAw3, AAw5, AAw10], ['CE', CEw1, CEw3, CEw5, CEw10], ['PL', PLw1, PLw3, PLw5, PLw10], ['Dockerin', dockw1, dockw3, dockw5, dockw10], ['SLH', SLHw1, SLHw3, SLHw5, SLHw10], ['cohesin', cohw1, cohw3, cohw5, cohw10]])       
dfpiechart2=pd.DataFrame(data=pie_chart_array[1:,1:], index=pie_chart_array[1:, 0], columns=pie_chart_array[0,1:])
#change values to floats
for x in dfpiechart2:
    dfpiechart2[x]=dfpiechart2[x].astype(float)
###calculate percentages
for x in dfpiechart2:
    newstr=x+'_pct'
    dfpiechart2[newstr]=dfpiechart2[x]/dfpiechart2[x].sum()*100
fontsize=14
plt.rc('font', size=fontsize)
fig_pie, ((ax_pie1, ax_pie3), (ax_pie5, ax_pie10)) = plt.subplots(2, 2, figsize=(12,12))
labels=list(dfpiechart2.index)
piecolors=sns.xkcd_palette(['seafoam', 'bubble gum pink', 'pastel yellow', 'sky', 'orchid', 'pinky red', 'pale lilac'])
pie_handle_list=[]
colct=0
for x in labels:
    pie_handle_list.append(patches.Patch(color=piecolors[colct], label=x))
    colct+=1   
ax_pie1.pie(dfpiechart2['One_pct'], colors=piecolors, startangle=90)
ax_pie1.axis('equal')
ax_pie1.set_title('Week one', size=20)
ax_pie3.pie(dfpiechart2['Three_pct'], colors=piecolors, startangle=90)
ax_pie3.axis('equal')
ax_pie3.set_title('Week three', size=20)
ax_pie5.pie(dfpiechart2['Five_pct'], colors=piecolors, startangle=90)
ax_pie5.axis('equal')
ax_pie5.set_title('Week five', size=20)
ax_pie10.pie(dfpiechart2['Ten_pct'], colors=piecolors, startangle=90)
ax_pie10.axis('equal')
ax_pie10.set_title('Week ten', size=20)
legend=plt.legend(handles=pie_handle_list, loc=(-1.15,-.3), title='Enzyme class', ncol=4, fontsize=20)
plt.setp(legend.get_title(),fontsize=20)
plt.text(-3.75,3.75,'b)', fontsize=28)
plt.show()
fig_pie.savefig('OUT_pie/enzyme_class_by_week.png', dpi=500)  
plt.clf()
plt.close(fig_pie)  
##########################################KINGDOM_PHYLUM############################################################################################
print('Creating databases for taxonomic represntation of enzyme classes')
poplist=[]
sum_activity_taxa_dict_total={}
sum_activity_taxa_dict_w1={} 
sum_activity_taxa_dict_w3={} 
sum_activity_taxa_dict_w5={} 
sum_activity_taxa_dict_w10={}
sum_activity_taxa_dict_molpct_w1={} 
sum_activity_taxa_dict_molpct_w3={} 
sum_activity_taxa_dict_molpct_w5={} 
sum_activity_taxa_dict_molpct_w10={} 
sum_activity_taxa_dict_molpct_error_w1={} 
sum_activity_taxa_dict_molpct_error_w3={} 
sum_activity_taxa_dict_molpct_error_w5={} 
sum_activity_taxa_dict_molpct_error_w10={}
for k in sum_activity_dict:
    sum_activity_taxa_dict_total[k]=[]
    sum_activity_taxa_dict_w1[k]=[]
    sum_activity_taxa_dict_w3[k]=[] 
    sum_activity_taxa_dict_w5[k]=[] 
    sum_activity_taxa_dict_w10[k]=[]
    sum_activity_taxa_dict_molpct_w1[k]={}
    sum_activity_taxa_dict_molpct_w3[k]={} 
    sum_activity_taxa_dict_molpct_w5[k]={} 
    sum_activity_taxa_dict_molpct_w10[k]={}
    sum_activity_taxa_dict_molpct_error_w1[k]={}
    sum_activity_taxa_dict_molpct_error_w3[k]={} 
    sum_activity_taxa_dict_molpct_error_w5[k]={} 
    sum_activity_taxa_dict_molpct_error_w10[k]={}
    
w1list=['W1_C201','W1_3_C429','W1_4_C429']    
w3list=['W3_C201','W3_2_C429','W3_1_C429']
w5list=['W5_C201','W5_25_C429','W5_100_C429']
w10list=['W10_C201','W10_100_i_C429','W10_100_ii_C429']
#create a list of taxa that is present within ALL samples, so that they can be assigned colours to ensure##
##homogeneity across graphs - a taxonomic colour index of sorts####
taxonomic_level='phylum'
total_relevant_taxa_list=[]
####kingdom_phylum#####
for x in sum_activity_dict:
    for k in crossref_molarpct:
        x1=str(x)+' '
        x2=str(x)+','
        if bool(re.search(x1, k))==True or bool(re.search(x2, k))==True:
             sum_activity_taxa_dict_total[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in total_relevant_taxa_list:
                 total_relevant_taxa_list.append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             w1temp_1=crossref_molarpct[k][w1list[0]]
             w1temp_2=crossref_molarpct[k][w1list[1]]
             w1temp_3=crossref_molarpct[k][w1list[2]]
             w1temp_avg=(w1temp_1+w1temp_2+w1temp_3)/3
             if not w1temp_1 == 0 or not w1temp_2 == 0 or not w1temp_3 == 0:
                sum_activity_taxa_dict_w1[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w1temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] in sum_activity_taxa_dict_molpct_w1[x]: 
                     sum_activity_taxa_dict_molpct_w1[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w1temp_avg+sum_activity_taxa_dict_molpct_w1[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in sum_activity_taxa_dict_molpct_w1[x]:
                     sum_activity_taxa_dict_molpct_w1[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w1temp_avg                    
             w3temp_1=crossref_molarpct[k][w3list[0]]
             w3temp_2=crossref_molarpct[k][w3list[1]]
             w3temp_3=crossref_molarpct[k][w3list[2]]
             w3temp_avg=(w3temp_1+w3temp_2+w3temp_3)/3
             if not w3temp_1 == 0 or not w3temp_2 == 0 or not w3temp_3 == 0:
                sum_activity_taxa_dict_w3[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w3temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] in sum_activity_taxa_dict_molpct_w3[x]: 
                     sum_activity_taxa_dict_molpct_w3[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w3temp_avg+sum_activity_taxa_dict_molpct_w3[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in sum_activity_taxa_dict_molpct_w3[x]:
                     sum_activity_taxa_dict_molpct_w3[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w3temp_avg
             w5temp_1=crossref_molarpct[k][w5list[0]]
             w5temp_2=crossref_molarpct[k][w5list[1]]
             w5temp_3=crossref_molarpct[k][w5list[2]]
             w5temp_avg=(w5temp_1+w5temp_2+w5temp_3)/3
             if not w5temp_1 == 0 or not w5temp_2 == 0 or not w5temp_3 == 0:
                sum_activity_taxa_dict_w5[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w5temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] in sum_activity_taxa_dict_molpct_w5[x]: 
                     sum_activity_taxa_dict_molpct_w5[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w5temp_avg+sum_activity_taxa_dict_molpct_w5[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]not in sum_activity_taxa_dict_molpct_w5[x]:
                     sum_activity_taxa_dict_molpct_w5[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w5temp_avg                     
             w10temp_1=crossref_molarpct[k][w10list[0]]
             w10temp_2=crossref_molarpct[k][w10list[1]]
             w10temp_3=crossref_molarpct[k][w10list[2]]
             w10temp_avg=(w10temp_1+w10temp_2+w10temp_3)/3
             if not w10temp_1 == 0 or not w10temp_2 == 0 or not w10temp_3 == 0:
                sum_activity_taxa_dict_w10[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w10temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] in sum_activity_taxa_dict_molpct_w10[x]: 
                     sum_activity_taxa_dict_molpct_w10[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w10temp_avg+sum_activity_taxa_dict_molpct_w10[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]        
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in sum_activity_taxa_dict_molpct_w10[x]:
                     sum_activity_taxa_dict_molpct_w10[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w10temp_avg
             for k2 in crossref_molarpct[k]:
                 if not k2 == 'Taxonomy':
                     if k2 in sum_activity_dict[x]:
                         sum_activity_dict[x][k2]=(crossref_molarpct[k][k2]+sum_activity_dict[x][k2])
                     elif k2 not in sum_activity_dict[x]:
                         sum_activity_dict[x][k2]=crossref_molarpct[k][k2]
##########################################CLASS############################################################################################
                                                                #
                                                                #
                                                                #
                                                 ####Whole_proteome_by_phylum####
taxonomic_level='phylum'
proteome_by_phylum_one={}
proteome_by_phylum_three={}
proteome_by_phylum_five={}
proteome_by_phylum_ten={}
proteome_by_phylum_one['NA']=0
proteome_by_phylum_three['NA']=0
proteome_by_phylum_five['NA']=0
proteome_by_phylum_ten['NA']=0
for x in df_master_mean.index:
    if x in taxa_master_dict:
         if taxa_master_dict[x]['Taxonomy'][taxonomic_level] not in proteome_by_phylum_one:
             proteome_by_phylum_one[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=0
             proteome_by_phylum_three[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=0
             proteome_by_phylum_five[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=0
             proteome_by_phylum_ten[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=0
         else: 
             proteome_by_phylum_one[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=proteome_by_phylum_one[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]+df_master_mean.loc[x]['One']
             proteome_by_phylum_three[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=proteome_by_phylum_three[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]+df_master_mean.loc[x]['Three']
             proteome_by_phylum_five[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=proteome_by_phylum_five[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]+df_master_mean.loc[x]['Five']
             proteome_by_phylum_ten[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]=proteome_by_phylum_ten[taxa_master_dict[x]['Taxonomy'][taxonomic_level]]+df_master_mean.loc[x]['Ten']
    elif x not in taxa_master_dict:
        proteome_by_phylum_one['NA']=proteome_by_phylum_one['NA']+df_master_mean.loc[x]['One']
        proteome_by_phylum_three['NA']=proteome_by_phylum_three['NA']+df_master_mean.loc[x]['Three']
        proteome_by_phylum_five['NA']=proteome_by_phylum_five['NA']+df_master_mean.loc[x]['Five']
        proteome_by_phylum_ten['NA']=proteome_by_phylum_ten['NA']+df_master_mean.loc[x]['Ten']

onedel=[]
threedel=[]
fivedel=[]
tendel=[]

for x in proteome_by_phylum_one:
    if proteome_by_phylum_one[x] < 0.5:
         onedel.append(x)
for x in onedel:
    del proteome_by_phylum_one[x]
    
for x in proteome_by_phylum_three:
    if proteome_by_phylum_three[x] < 0.5:
         threedel.append(x)
for x in threedel:
    del proteome_by_phylum_three[x]

for x in proteome_by_phylum_five:
    if proteome_by_phylum_five[x] < 0.5:
         fivedel.append(x)
for x in fivedel:
    del proteome_by_phylum_five[x]

for x in proteome_by_phylum_ten:
    if proteome_by_phylum_ten[x] < 0.5:
         tendel.append(x)
for x in tendel:
    del proteome_by_phylum_ten[x] 

#insert empty keys
keylist=list(proteome_by_phylum_one.keys())+list(proteome_by_phylum_three.keys())+list(proteome_by_phylum_five.keys())+list(proteome_by_phylum_ten.keys())
for x in set(keylist):
    if x not in proteome_by_phylum_one:
        proteome_by_phylum_one[x]=0
    if x not in proteome_by_phylum_three:
        proteome_by_phylum_three[x]=0
    if x not in proteome_by_phylum_five:
        proteome_by_phylum_five[x]=0
    if x not in proteome_by_phylum_ten:
        proteome_by_phylum_ten[x]=0

df_proteome_taxa=pd.DataFrame([proteome_by_phylum_one, proteome_by_phylum_three, proteome_by_phylum_five, proteome_by_phylum_ten], index=(['One', 'Three', 'Five', 'Ten'])).T
for x in df_proteome_taxa:
    df_proteome_taxa[x]=df_proteome_taxa[x].astype(float)
###calculate percentages
for x in df_proteome_taxa:
    newstr=x+'_pct'
    df_proteome_taxa[newstr]=df_proteome_taxa[x]/df_proteome_taxa[x].sum()*100
fontsize=14
plt.rc('font', size=fontsize)
fig_pie, ((ax_pie1, ax_pie3), (ax_pie5, ax_pie10)) = plt.subplots(2, 2, figsize=(12,12))
labels=list(df_proteome_taxa.index)
#piecolors=sns.xkcd_palette(['seafoam', 'bubble gum pink', 'pastel yellow', 'sky', 'orchid', 'pinky red', 'pale lilac'])
#piecolors=sns.color_palette("hls", len(set(keylist)))
halfcol=int(len(set(keylist))/2)+1
#piecolors=(sns.hls_palette(halfcol, h=.75, s=1, l=.5))+sns.hls_palette(halfcol, h=.75, s=.6, l=.7)
piecolors=sns.hls_palette(8, l=.5, s=1)+sns.hls_palette(8, l=.45, s=.75)+sns.hls_palette(8, l=.65, s=.9)+sns.hls_palette(8, l=.8, s=.6)+sns.color_palette("Paired", 12)+sns.color_palette("colorblind")+sns.color_palette("pastel")
pie_handle_list=[]
colct=0
for x in labels:
    pie_handle_list.append(patches.Patch(color=piecolors[colct], label=x))
    colct+=1   
ax_pie1.pie(df_proteome_taxa['One_pct'], colors=piecolors, startangle=90)
ax_pie1.axis('equal')
ax_pie1.set_title('Week one', size=fontsize)
ax_pie3.pie(df_proteome_taxa['Three_pct'], colors=piecolors, startangle=90)
ax_pie3.axis('equal')
ax_pie3.set_title('Week three', size=fontsize)
ax_pie5.pie(df_proteome_taxa['Five_pct'], colors=piecolors, startangle=90)
ax_pie5.axis('equal')
ax_pie5.set_title('Week five', size=fontsize)
ax_pie10.pie(df_proteome_taxa['Ten_pct'], colors=piecolors, startangle=90)
ax_pie10.axis('equal')
ax_pie10.set_title('Week ten', size=fontsize)
if taxonomic_level=='phylum':
     plt.legend(handles=pie_handle_list, loc=(-1.1,-.35), title=str(taxonomic_level.title()), ncol=4, fontsize=fontsize)
if taxonomic_level=='class':
     plt.legend(handles=pie_handle_list, loc=(-1.3,-.55), title=str(taxonomic_level.title()), ncol=4, fontsize=fontsize)
if taxonomic_level=='order':
     plt.legend(handles=pie_handle_list, loc=(-1.3,-.55), title=str(taxonomic_level.title()), ncol=4, fontsize=fontsize)    
if taxonomic_level=='genus':
     plt.legend(handles=pie_handle_list, loc=(-1.3,-.55), title=str(taxonomic_level.title()), ncol=4, fontsize=fontsize)  
plt.text(-3.75,3.75,'b)', fontsize=28)
plt.show()
fig_pie.savefig('OUT_pie/proteome_by_taxa_'+taxonomic_level+'_level.png') 
plt.clf()
plt.close(fig_pie)                                                                  #
                                                                #
                                                                #                                                    
##########################################CLASS############################################################################################
##########################################CLASS############################################################################################
print('Creating databases for taxonomic represntation of enzyme classes')
###clean sum_activity_dict####
sum_activity_dict={}
sum_activity_dict=sum_activity_dict_skeleton_copy

poplist=[]
sum_activity_taxa_dict_total_class={}
sum_activity_taxa_dict_w1_class={} 
sum_activity_taxa_dict_w3_class={} 
sum_activity_taxa_dict_w5_class={} 
sum_activity_taxa_dict_w10_class={}
sum_activity_taxa_dict_molpct_w1_class={} 
sum_activity_taxa_dict_molpct_w3_class={} 
sum_activity_taxa_dict_molpct_w5_class={} 
sum_activity_taxa_dict_molpct_w10_class={} 
sum_activity_taxa_dict_molpct_error_w1_class={} 
sum_activity_taxa_dict_molpct_error_w3_class={} 
sum_activity_taxa_dict_molpct_error_w5_class={} 
sum_activity_taxa_dict_molpct_error_w10_class={}
for k in sum_activity_dict:
    sum_activity_taxa_dict_total_class[k]=[]
    sum_activity_taxa_dict_w1_class[k]=[]
    sum_activity_taxa_dict_w3_class[k]=[] 
    sum_activity_taxa_dict_w5_class[k]=[] 
    sum_activity_taxa_dict_w10_class[k]=[]
    sum_activity_taxa_dict_molpct_w1_class[k]={}
    sum_activity_taxa_dict_molpct_w3_class[k]={} 
    sum_activity_taxa_dict_molpct_w5_class[k]={} 
    sum_activity_taxa_dict_molpct_w10_class[k]={}
    sum_activity_taxa_dict_molpct_error_w1_class[k]={}
    sum_activity_taxa_dict_molpct_error_w3_class[k]={} 
    sum_activity_taxa_dict_molpct_error_w5_class[k]={} 
    sum_activity_taxa_dict_molpct_error_w10_class[k]={}

taxonomic_level='class'
total_relevant_taxa_list_class=[]
for x in sum_activity_dict:
    for k in crossref_molarpct:
        x1=str(x)+' '
        x2=str(x)+','
        if bool(re.search(x1, k))==True or bool(re.search(x2, k))==True:
             sum_activity_taxa_dict_total_class[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in total_relevant_taxa_list_class:
                 total_relevant_taxa_list_class.append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             w1temp_1=crossref_molarpct[k][w1list[0]]
             w1temp_2=crossref_molarpct[k][w1list[1]]
             w1temp_3=crossref_molarpct[k][w1list[2]]
             w1temp_avg=(w1temp_1+w1temp_2+w1temp_3)/3
             if not w1temp_1 == 0 or not w1temp_2 == 0 or not w1temp_3 == 0:
                sum_activity_taxa_dict_w1_class[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w1temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] in sum_activity_taxa_dict_molpct_w1_class[x]: 
                     sum_activity_taxa_dict_molpct_w1_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w1temp_avg+sum_activity_taxa_dict_molpct_w1_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in sum_activity_taxa_dict_molpct_w1_class[x]:
                     sum_activity_taxa_dict_molpct_w1_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w1temp_avg                    
             w3temp_1=crossref_molarpct[k][w3list[0]]
             w3temp_2=crossref_molarpct[k][w3list[1]]
             w3temp_3=crossref_molarpct[k][w3list[2]]
             w3temp_avg=(w3temp_1+w3temp_2+w3temp_3)/3
             if not w3temp_1 == 0 or not w3temp_2 == 0 or not w3temp_3 == 0:
                sum_activity_taxa_dict_w3_class[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w3temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] in sum_activity_taxa_dict_molpct_w3_class[x]: 
                     sum_activity_taxa_dict_molpct_w3_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w3temp_avg+sum_activity_taxa_dict_molpct_w3_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in sum_activity_taxa_dict_molpct_w3_class[x]:
                     sum_activity_taxa_dict_molpct_w3_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w3temp_avg
             w5temp_1=crossref_molarpct[k][w5list[0]]
             w5temp_2=crossref_molarpct[k][w5list[1]]
             w5temp_3=crossref_molarpct[k][w5list[2]]
             w5temp_avg=(w5temp_1+w5temp_2+w5temp_3)/3
             if not w5temp_1 == 0 or not w5temp_2 == 0 or not w5temp_3 == 0:
                sum_activity_taxa_dict_w5_class[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w5temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] in sum_activity_taxa_dict_molpct_w5_class[x]: 
                     sum_activity_taxa_dict_molpct_w5_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w5temp_avg+sum_activity_taxa_dict_molpct_w5_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]not in sum_activity_taxa_dict_molpct_w5_class[x]:
                     sum_activity_taxa_dict_molpct_w5_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w5temp_avg                     
             w10temp_1=crossref_molarpct[k][w10list[0]]
             w10temp_2=crossref_molarpct[k][w10list[1]]
             w10temp_3=crossref_molarpct[k][w10list[2]]
             w10temp_avg=(w10temp_1+w10temp_2+w10temp_3)/3
             if not w10temp_1 == 0 or not w10temp_2 == 0 or not w10temp_3 == 0:
                sum_activity_taxa_dict_w10_class[x].append(taxa_master_dict_named[k]['Taxonomy'][taxonomic_level])
             if not w10temp_avg == 0:
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]in sum_activity_taxa_dict_molpct_w10_class[x]: 
                     sum_activity_taxa_dict_molpct_w10_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w10temp_avg+sum_activity_taxa_dict_molpct_w10_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]        
                 if taxa_master_dict_named[k]['Taxonomy'][taxonomic_level] not in sum_activity_taxa_dict_molpct_w10_class[x]:
                     sum_activity_taxa_dict_molpct_w10_class[x][taxa_master_dict_named[k]['Taxonomy'][taxonomic_level]]=w10temp_avg
             for k2 in crossref_molarpct[k]:
                if not k2 == 'Taxonomy':
                     if k2 in sum_activity_dict[x]:
                         sum_activity_dict[x][k2]=(crossref_molarpct[k][k2]+sum_activity_dict[x][k2])
                     elif k2 not in sum_activity_dict[x]:
                         sum_activity_dict[x][k2]=crossref_molarpct[k][k2]
#######################################################################################################                    
for x in sum_activity_dict:
    if sum_activity_dict[x]=={}:
        poplist.append(x)
for key in poplist:
    sum_activity_dict.pop(key)  
    sum_activity_taxa_dict_total.pop(key)
    sum_activity_taxa_dict_w1.pop(key)
    sum_activity_taxa_dict_w3.pop(key)
    sum_activity_taxa_dict_w5.pop(key)
    sum_activity_taxa_dict_w10.pop(key)
    sum_activity_taxa_dict_molpct_w1.pop(key)
    sum_activity_taxa_dict_molpct_w3.pop(key)
    sum_activity_taxa_dict_molpct_w5.pop(key)
    sum_activity_taxa_dict_molpct_w10.pop(key)
    sum_activity_taxa_dict_total_class.pop(key)
    sum_activity_taxa_dict_w1_class.pop(key)
    sum_activity_taxa_dict_w3_class.pop(key)
    sum_activity_taxa_dict_w5_class.pop(key)
    sum_activity_taxa_dict_w10_class.pop(key)
    sum_activity_taxa_dict_molpct_w1_class.pop(key)
    sum_activity_taxa_dict_molpct_w3_class.pop(key) 
    sum_activity_taxa_dict_molpct_w5_class.pop(key)
    sum_activity_taxa_dict_molpct_w10_class.pop(key)
    
taxa_list_map_total={}
taxa_list_map_w1={}
taxa_list_map_w3={}
taxa_list_map_w5={}
taxa_list_map_w10={}
taxa_list_map_total_class={}
taxa_list_map_w1_class={}
taxa_list_map_w3_class={}
taxa_list_map_w5_class={}
taxa_list_map_w10_class={}

count_map=0
for x in sum_activity_taxa_dict_total:
    for x2 in sum_activity_taxa_dict_total[x]:
        if x2 not in taxa_list_map_total:
            count_map+=1
            taxa_list_map_total[x2]=count_map

count_map=0
for x in sum_activity_taxa_dict_total_class:
    for x2 in sum_activity_taxa_dict_total_class[x]:
        if x2 not in taxa_list_map_total:
            count_map+=1
            taxa_list_map_total_class[x2]=count_map
                               
count_map=0
for x in sum_activity_taxa_dict_w1:
    for x2 in sum_activity_taxa_dict_w1[x]:
        if x2 not in taxa_list_map_w1:
            count_map+=1
            taxa_list_map_w1[x2]=count_map     
 
count_map=0
for x in sum_activity_taxa_dict_w1_class:
    for x2 in sum_activity_taxa_dict_w1_class[x]:
        if x2 not in taxa_list_map_w1_class:
            count_map+=1
            taxa_list_map_w1_class[x2]=count_map  
                           
count_map=0
for x in sum_activity_taxa_dict_w3:
    for x2 in sum_activity_taxa_dict_w3[x]:
        if x2 not in taxa_list_map_w3:
            count_map+=1
            taxa_list_map_w3[x2]=count_map  

count_map=0
for x in sum_activity_taxa_dict_w3_class:
    for x2 in sum_activity_taxa_dict_w3_class[x]:
        if x2 not in taxa_list_map_w3_class:
            count_map+=1
            taxa_list_map_w3_class[x2]=count_map
            
count_map=0
for x in sum_activity_taxa_dict_w5:
    for x2 in sum_activity_taxa_dict_w5[x]:
        if x2 not in taxa_list_map_w5:
            count_map+=1
            taxa_list_map_w5[x2]=count_map
 
count_map=0
for x in sum_activity_taxa_dict_w5_class:
    for x2 in sum_activity_taxa_dict_w5_class[x]:
        if x2 not in taxa_list_map_w5_class:
            count_map+=1
            taxa_list_map_w5_class[x2]=count_map
                           
count_map=0
for x in sum_activity_taxa_dict_w10:
    for x2 in sum_activity_taxa_dict_w10[x]:
        if x2 not in taxa_list_map_w10:
            count_map+=1
            taxa_list_map_w10[x2]=count_map                        
 
count_map=0
for x in sum_activity_taxa_dict_w10_class:
    for x2 in sum_activity_taxa_dict_w10_class[x]:
        if x2 not in taxa_list_map_w10_class:
            count_map+=1
            taxa_list_map_w10_class[x2]=count_map
           
print('There are ' + str(len(taxa_list_map_total.keys())) + ' high level taxanomic groups putatively responsible for protein production based on homology with these filtering parameters')    
print('There are ' + str(len(taxa_list_map_w1.keys())) + ' high level taxanomic groups putatively responsible for protein production in week 1 based on homology with these filtering parameters')
print('There are ' + str(len(taxa_list_map_w3.keys())) + ' high level taxanomic groups putatively responsible for protein production in week 3 based on homology with these filtering parameters')
print('There are ' + str(len(taxa_list_map_w5.keys())) + ' high level taxanomic groups putatively responsible for protein production in week 5 based on homology with these filtering parameters')
print('There are ' + str(len(taxa_list_map_w10.keys())) + ' high level taxanomic groups putatively responsible for protein production in week 10 based on homology with these filtering parameters')
#######Re-format to a a dict that can be transferred to dataframe##############
###############################################################################
print('Assigning representative taxonomic values to enzyme classes - total and by week')
sum_activity_taxa_dict_total_format={}
sum_activity_taxa_dict_w1_format={}
sum_activity_taxa_dict_w3_format={}
sum_activity_taxa_dict_w5_format={}
sum_activity_taxa_dict_w10_format={}

sum_activity_taxa_dict_total_format_class={}
sum_activity_taxa_dict_w1_format_class={}
sum_activity_taxa_dict_w3_format_class={}
sum_activity_taxa_dict_w5_format_class={}
sum_activity_taxa_dict_w10_format_class={}

for x, v in sum_activity_taxa_dict_total.items():
    if x not in sum_activity_taxa_dict_total_format:
         sum_activity_taxa_dict_total_format[x]={}
    for x3 in taxa_list_map_total:
             sum_activity_taxa_dict_total_format[x][x3]=v.count(x3)

for x, v in sum_activity_taxa_dict_total_class.items():
    if x not in sum_activity_taxa_dict_total_format_class:
         sum_activity_taxa_dict_total_format_class[x]={}
    for x3 in taxa_list_map_total_class:
             sum_activity_taxa_dict_total_format_class[x][x3]=v.count(x3)
             
for x, v in sum_activity_taxa_dict_w1.items():
    if x not in sum_activity_taxa_dict_w1_format:
         sum_activity_taxa_dict_w1_format[x]={}
    for x3 in taxa_list_map_w1:
             sum_activity_taxa_dict_w1_format[x][x3]=v.count(x3)  
 
for x, v in sum_activity_taxa_dict_w1_class.items():
    if x not in sum_activity_taxa_dict_w1_format_class:
         sum_activity_taxa_dict_w1_format_class[x]={}
    for x3 in taxa_list_map_w1_class:
             sum_activity_taxa_dict_w1_format_class[x][x3]=v.count(x3)
            
for x, v in sum_activity_taxa_dict_w3.items():
    if x not in sum_activity_taxa_dict_w3_format:
         sum_activity_taxa_dict_w3_format[x]={}
    for x3 in taxa_list_map_w3:
             sum_activity_taxa_dict_w3_format[x][x3]=v.count(x3)  

for x, v in sum_activity_taxa_dict_w3_class.items():
    if x not in sum_activity_taxa_dict_w3_format_class:
         sum_activity_taxa_dict_w3_format_class[x]={}
    for x3 in taxa_list_map_w3_class:
             sum_activity_taxa_dict_w3_format_class[x][x3]=v.count(x3) 

for x, v in sum_activity_taxa_dict_w5.items():
    if x not in sum_activity_taxa_dict_w5_format:
         sum_activity_taxa_dict_w5_format[x]={}
    for x3 in taxa_list_map_w5:
             sum_activity_taxa_dict_w5_format[x][x3]=v.count(x3)  

for x, v in sum_activity_taxa_dict_w5_class.items():
    if x not in sum_activity_taxa_dict_w5_format_class:
         sum_activity_taxa_dict_w5_format_class[x]={}
    for x3 in taxa_list_map_w5_class:
             sum_activity_taxa_dict_w5_format_class[x][x3]=v.count(x3) 
             
for x, v in sum_activity_taxa_dict_w10.items():
    if x not in sum_activity_taxa_dict_w10_format:
         sum_activity_taxa_dict_w10_format[x]={}
    for x3 in taxa_list_map_w10:
             sum_activity_taxa_dict_w10_format[x][x3]=v.count(x3) 
              
for x, v in sum_activity_taxa_dict_w10_class.items():
    if x not in sum_activity_taxa_dict_w10_format_class:
         sum_activity_taxa_dict_w10_format_class[x]={}
    for x3 in taxa_list_map_w10_class:
             sum_activity_taxa_dict_w10_format_class[x][x3]=v.count(x3)
###Add missing groups with 0 values for molar percent databases####
###homogenising the data averts issues when converting dict to pandas df###
###barchart script DOES NOT like NANs#######
def replace_missing_taxa (dictionary, total_taxa_list):
     for cazy in dictionary:
         templist=[]
         for tax in total_taxa_list:
            if tax not in dictionary[cazy]:
                templist.append(tax)
         for tax2 in templist:
             dictionary[cazy][tax2]=0

sum_activity_taxa_dict_molpct_w1_raw=sum_activity_taxa_dict_molpct_w1 
sum_activity_taxa_dict_molpct_w3_raw=sum_activity_taxa_dict_molpct_w3
sum_activity_taxa_dict_molpct_w5_raw=sum_activity_taxa_dict_molpct_w5
sum_activity_taxa_dict_molpct_w10_raw=sum_activity_taxa_dict_molpct_w10

sum_activity_taxa_dict_molpct_w1_class_raw=sum_activity_taxa_dict_molpct_w1_class
sum_activity_taxa_dict_molpct_w3_class_raw=sum_activity_taxa_dict_molpct_w3_class
sum_activity_taxa_dict_molpct_w5_class_raw=sum_activity_taxa_dict_molpct_w5_class
sum_activity_taxa_dict_molpct_w10_class_raw=sum_activity_taxa_dict_molpct_w10_class

####_raw files for circos 
           
replace_missing_taxa(sum_activity_taxa_dict_molpct_w1, total_relevant_taxa_list)
replace_missing_taxa(sum_activity_taxa_dict_molpct_w3, total_relevant_taxa_list)
replace_missing_taxa(sum_activity_taxa_dict_molpct_w5, total_relevant_taxa_list)
replace_missing_taxa(sum_activity_taxa_dict_molpct_w10, total_relevant_taxa_list)

replace_missing_taxa(sum_activity_taxa_dict_molpct_w1_class, total_relevant_taxa_list_class)
replace_missing_taxa(sum_activity_taxa_dict_molpct_w3_class, total_relevant_taxa_list_class)
replace_missing_taxa(sum_activity_taxa_dict_molpct_w5_class, total_relevant_taxa_list_class)
replace_missing_taxa(sum_activity_taxa_dict_molpct_w10_class, total_relevant_taxa_list_class)

df_sum_activity_taxa_dict_total_format=pd.DataFrame.from_dict(sum_activity_taxa_dict_total_format, orient='index')
df_sum_activity_taxa_dict_w1_format=pd.DataFrame.from_dict(sum_activity_taxa_dict_w1_format, orient='index')
df_sum_activity_taxa_dict_w3_format=pd.DataFrame.from_dict(sum_activity_taxa_dict_w3_format, orient='index')
df_sum_activity_taxa_dict_w5_format=pd.DataFrame.from_dict(sum_activity_taxa_dict_w5_format, orient='index')
df_sum_activity_taxa_dict_w10_format=pd.DataFrame.from_dict(sum_activity_taxa_dict_w10_format, orient='index')

df_sum_activity_taxa_dict_total_format_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_total_format_class, orient='index')
df_sum_activity_taxa_dict_w1_format_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_w1_format_class, orient='index')
df_sum_activity_taxa_dict_w3_format_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_w3_format_class, orient='index')
df_sum_activity_taxa_dict_w5_format_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_w5_format_class, orient='index')
df_sum_activity_taxa_dict_w10_format_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_w10_format_class, orient='index')

df_sum_activity_taxa_dict_molpct_w1=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w1, orient='index')
df_sum_activity_taxa_dict_molpct_w3=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w3, orient='index')
df_sum_activity_taxa_dict_molpct_w5=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w5, orient='index')
df_sum_activity_taxa_dict_molpct_w10=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w10, orient='index')

df_sum_activity_taxa_dict_molpct_w1_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w1_class, orient='index')
df_sum_activity_taxa_dict_molpct_w3_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w3_class, orient='index')
df_sum_activity_taxa_dict_molpct_w5_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w5_class, orient='index')
df_sum_activity_taxa_dict_molpct_w10_class=pd.DataFrame.from_dict(sum_activity_taxa_dict_molpct_w10_class, orient='index')

############### CIRCOS #######################CIRCOS###########################CIRCOS######################################################
def int_generator(mylist):
    out=[]
    for n in mylist:
        try:
            out.append(int(float(n)))
        except ValueError:
            out.append(n)
            continue
    return out

def seaborn_pal_to_rgb(pal):
     outlist=[]
     for color in pal:
        temptup=''
        col=0
        for val in color:
            if col < 3:
                 nuval=val*255
                 temptup=temptup+(str(int(nuval))+',')
                 col+=1
        outlist.append(temptup[:-1])
     return outlist

#AA=pink, CE=yellows, SLH=reds, GH=greens, PL=blues, CBM=oranges
CEc=seaborn_pal_to_rgb(sns.light_palette("yellow", 12)[::-1])
PLc=seaborn_pal_to_rgb(sns.light_palette("blue", 10)[::-1])
GHc=seaborn_pal_to_rgb(sns.light_palette("highlighter green",65, input="xkcd")[::-1])
AAc=seaborn_pal_to_rgb(sns.light_palette("red",10)[::-1])
SLHc=seaborn_pal_to_rgb(sns.light_palette((0, 90, 50),5, input="husl")[::-1])
CBMc=seaborn_pal_to_rgb(sns.light_palette("cyan", 35, input="xkcd")[::-1])

CEcount=0
PLcount=0
GHcount=0
AAcount=0
SLHcount=0
CBMcount=0

enzyme_class_list=list(df_sum_activity_taxa_dict_molpct_w1.index)
enzyme_class_colors=[]
for x in enzyme_class_list:
         if x.startswith('CE'):
             enzyme_class_colors.append(CEc[CEcount])
             CEcount+=1
         if x.startswith('PL'):
             enzyme_class_colors.append(PLc[PLcount])
             PLcount+=1
         if x.startswith('GH'):
             enzyme_class_colors.append(GHc[GHcount])
             GHcount+=1
         if x.startswith('AA'):
             enzyme_class_colors.append(AAc[AAcount])
             AAcount+=1  
         if x.startswith('SLH'):
             enzyme_class_colors.append(SLHc[SLHcount])
             SLHcount+=1
         if x.startswith('CBM'):
             enzyme_class_colors.append(CBMc[CBMcount])
             CBMcount+=1
enzyme_class_colo_dict=dict(zip(enzyme_class_list, enzyme_class_colors))

            
phyla_colors=['36,122,253',
 '254,44,84',
 '255,255,45',
 '117,253,99',
 '255,165,0',
 '178,91,178',
 '150,150,150',
 '229,25,127',
 '229,25,229',
 '25,25,229',
 '229,25,25',
 '25,127,229',
 '25,229,229',
 '127,229,25',
 '25,229,127']

#phyla_list=list(df_sum_activity_taxa_dict_molpct_w1.sum(axis=0).rank(ascending=False, method='first').order().index)
phyla_list=list(df_sum_activity_taxa_dict_molpct_w1.sum(axis=0).rank(ascending=False, method='first').sort_values().index)
phyla_list=[x.replace('Candidatus_', '')+'_' if x.startswith('Candidatus_') else x for x in phyla_list]
phyla_colo_dict=dict(zip(phyla_list, phyla_colors))

tax_class_colors=['36,122,253',
 '254,44,84',
 '255,255,45',
 '117,253,99',
 '255,165,0',
 '178,91,178',
 '150,150,150',
 '255,170,255',
 '100,0,150',
 '25,25,229',
 '200,150,100',
 '225,245,255',
 '40,150,150',
 '127,229,25',
 '25,229,127',
 '229,25,229',
 '25,25,229',
 '255,220,220',
 '25,127,229',
 '15,200,200',
 '127,229,25',
 '25,229,127']

#cyto=229,25,127
#bac=254,44,84
#class_tax_list=list(df_sum_activity_taxa_dict_molpct_w1_class.sum(axis=0).rank(ascending=False, method='first').order().index)
class_tax_list=list(df_sum_activity_taxa_dict_molpct_w1_class.sum(axis=0).rank(ascending=False, method='first').sort_values().index)
class_tax_list=[x.replace('Candidatus_', '')+'_' if x.startswith('Candidatus_') else x for x in class_tax_list]
class_colo_dict=dict(zip(class_tax_list, tax_class_colors[:len(class_tax_list)]))


###This script is designed to be used with circos with the fllowing options enabled:
    #TICK col with row order, TICK row with col order
    #TICK size with row size, TICK row with col size
    #THICK thicknes, TIGHT spacing, SMALL radium, 
    #BOLD text, NO parallel
    #TRANSPARANCY=2
def dataframe_to_circos(inputdf, input_color_dict, enzyme_colo_dict, taxafilter, classfilter):
    multiplier=1000
    blanklist=['9999', '9998', '9997', '9996']
    circosdf=pd.DataFrame(index=inputdf.index)   
    for index in list(inputdf.columns):
         circosdf[index]=(inputdf[index]*multiplier).astype(int)
    circosdf.fillna(axis='columns', value=0).astype(int)
    circosdf.columns=[x.replace('Candidatus ', '')+'_' if x.startswith('Candidatus ') else x for x in circosdf.columns]
    #Sum rows into column, filter data frame with 0 values, write the remaining sums into a list for appending later 
    circosdf=circosdf[(circosdf.sum(axis=1)>0)]
    circosdf=circosdf[circosdf.columns[circosdf.sum(axis=0) >0]]
    ##############
    class_sum=(circosdf.sum(axis=1))
    taxa_sum=(circosdf.sum(axis=0))
    #filter_based_on_input_value
    classes_to_keep=class_sum/class_sum.sum()*100 > classfilter
    taxa_to_keep=taxa_sum/taxa_sum.sum()*100 > taxafilter
    #########RAW file for TROUBLESHOOTING####
    #circosdfraw=deepcopy(circosdf)
    ##########################################
    
    #Create Other_CAZy and Other_taxon (row and column respectively) - with the sums off all data to be filtered
    classes_to_keep_N=len([tru for tru, fals in enumerate(classes_to_keep) if not fals])
    taxa_to_keep_N=len([tru for tru, fals in enumerate(taxa_to_keep) if not fals])
    ##check the sum of filtered values is not 0
    sum_total_value_filtered_groups=((circosdf[~(classes_to_keep)].sum()).sum()).sum()
    if not sum_total_value_filtered_groups == 0:
         circosdf.loc['CAZy-'+str(classes_to_keep_N)]=circosdf[~(classes_to_keep)].sum()
    circosdf['Other-'+str(taxa_to_keep_N)]=circosdf[circosdf.columns[~(taxa_to_keep)]].sum(axis=1)
    #Update filters class_sum and taxa_sum with the new row and column from above two lines
    class_sum2=(circosdf.sum(axis=1))
    taxa_sum2=(circosdf.sum(axis=0))   
    classes_to_keep=class_sum2/class_sum2.sum()*100 >= classfilter
    taxa_to_keep=taxa_sum2/taxa_sum2.sum()*100 >= taxafilter                         
    circosdf=circosdf[classes_to_keep]
    print(str(len([i for i, x in enumerate(class_sum2/class_sum2.sum()*100 > classfilter) if x]))+' enzyme classes remain of ' +str(len(class_sum2)))
    N_class_in_other=len([i for i, x in enumerate(class_sum2/class_sum2.sum()*100 > classfilter) if x])-(len(class_sum2)-1)
    print('Filtered enzyme classes less than '+str(classfilter)+'% (relative molar percent)')
    circosdf=circosdf[circosdf.columns[taxa_to_keep]]
    filteredclasssum=(circosdf.sum(axis=1))
    class_rank_order=filteredclasssum.rank(ascending=False, method='first')
    print(str(len([i for i, x in enumerate(taxa_sum2/taxa_sum2.sum()*100 > taxafilter) if x]))+' taxa remain of ' +str(len(taxa_sum)))
    N_taxa_in_other=len([i for i, x in enumerate(taxa_sum2/taxa_sum2.sum()*100 > taxafilter) if x])-(len(taxa_sum2)-1)
    print('Filtered taxa less than '+str(taxafilter)+'% (relative molar percent)')
    newtaxasum=(circosdf.sum(axis=0))
    taxa_rank_order=blanklist+list(newtaxasum.rank(ascending=False, method='first'))
    taxa_rank_order=int_generator(taxa_rank_order)
    #The[(0,0,0)] tuple is for the other cazy group
    class_color_order_column=[enzyme_colo_dict[x] for x in circosdf.index if x in enzyme_colo_dict]
    if not sum_total_value_filtered_groups == 0:
         class_color_order_column=class_color_order_column+['255,0,255']
    ####
    colnames=list(circosdf.columns)
    #fill NaNs with 0 and turn all values from floats to integers
    circosdf=circosdf.fillna(axis='columns', value=0).astype(int).reset_index()
    circosdf[111]=class_color_order_column
### class_rank_order=[x+int(max(taxa_rank_order[3:])) for x in class_rank_order]
    class_rank_order=[x+int(max(taxa_rank_order[4:])) for x in class_rank_order]
    class_rank_order=int_generator(class_rank_order)
    circosdf[9955]=class_rank_order
    circosdf['size']=list(filteredclasssum)
    #reorder columns
    colorder=[9955,'size',111, 'index']+colnames
    color_order=blanklist+[input_color_dict[x] for x in colorder[4:] if x in input_color_dict]
    if not len(taxa_sum2)-1 == len([i for i, x in enumerate(taxa_sum2/taxa_sum2.sum()*100 > taxafilter) if x]):
        color_order=color_order+['50,50,50']
    circosdf=circosdf[colorder]
    #change column names to taxa sum, so the dataframes can be concatenated
 ### taxa_rank_order=[new-classrankmax for new in taxa_rank_order]
    circosdf.columns=taxa_rank_order
    #prepend blank data to newtaxsum to fit the dataframe
    newtaxasum=blanklist+list(newtaxasum)
    tempdf=pd.DataFrame([newtaxasum,color_order, colorder], columns=taxa_rank_order)
    circosdf=tempdf.append(circosdf, ignore_index=True)
    return circosdf


def write_circos_file(inputdf, filenamex):
    file1path='Circos_chord_ideograms/'+str(filenamex)+'.circos'
    file2path='Circos_chord_ideograms/Formatted_'+str(filenamex)+'.circos'
    inputdf.to_csv(path_or_buf=file1path, sep=' ',  index=False)#, index_label='index')
    f=open(file1path)
    newfilelines=[]
    while True:
        line=f.readline().strip('\n')
        if line=='':
            break
        newfilelines.append(line)
    f.close()
    print('File construction complete')

#total_PMO_by_phylum_circos=dataframe_to_circos(df_sum_activity_taxa_dict_total_format,phyla_colo_dict,enzyme_class_colo_dict,1,1)
#    
#phyl_tax_filt=0
#phyl_class_filt=0.75
#
#w1_phylum_circos=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w1,phyla_colo_dict,enzyme_class_colo_dict,phyl_tax_filt,phyl_class_filt)  
#w3_phylum_circos=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w3,phyla_colo_dict,enzyme_class_colo_dict,phyl_tax_filt,phyl_class_filt) 
#w5_phylum_circos=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w5, phyla_colo_dict,enzyme_class_colo_dict,phyl_tax_filt,phyl_class_filt)  
#w10_phylum_circos=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w10, phyla_colo_dict,enzyme_class_colo_dict,phyl_tax_filt,phyl_class_filt) 
#
#class_tax_filt=1
#class_class_filt=0.85
#
#w1_phylum_circos_class=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w1_class, class_colo_dict,enzyme_class_colo_dict,class_tax_filt,class_class_filt)  
#w3_phylum_circos_class=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w3_class, class_colo_dict,enzyme_class_colo_dict,class_tax_filt,class_class_filt) 
#w5_phylum_circos_class=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w5_class, class_colo_dict,enzyme_class_colo_dict,class_tax_filt,class_class_filt)  
#w10_phylum_circos_class=dataframe_to_circos(df_sum_activity_taxa_dict_molpct_w10_class, class_colo_dict,enzyme_class_colo_dict,class_tax_filt,class_class_filt) 
#
#write_circos_file(total_PMO_by_phylum_circos, 'total_PMO_phylum')
#
#write_circos_file(w1_phylum_circos, 'w1_phylum_circos')
#write_circos_file(w3_phylum_circos, 'w3_phylum_circos')
#write_circos_file(w5_phylum_circos, 'w5_phylum_circos')
#write_circos_file(w10_phylum_circos, 'w10_phylum_circos')
#
#write_circos_file(w1_phylum_circos_class, 'w1_phylum_circos_class')
#write_circos_file(w3_phylum_circos_class, 'w3_phylum_circos_class')
#write_circos_file(w5_phylum_circos_class, 'w5_phylum_circos_class')
#write_circos_file(w10_phylum_circos_class, 'w10_phylum_circos_class')
############## CIRCOS #######################CIRCOS###########################CIRCOS######################################################
####CRICOS LEGEND####
#
#plot_these_legends=set(list(w1_phylum_circos_class.iloc[2][4:])+list(w3_phylum_circos_class.iloc[2][4:])+list(w5_phylum_circos_class.iloc[2][4:])+list(w10_phylum_circos_class.iloc[2][4:])+['Other'])
#remove=[]
#for x in plot_these_legends:
#    if x.startswith('Other-'):
#        remove.append(x)
#for y in remove:
#    plot_these_legends.remove(y)
#    
#class_colo_dict.update({'Other':'50,50,50'})
#
#changedict={'Gammaproteobacteria':'γ-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria'}
#
#patches5=[]    
#for x in class_colo_dict:
#    if x in plot_these_legends:
#        if x in changedict:
#             patches5.append(patches.Patch(facecolor=tuple([int(x)/255 for x in class_colo_dict[x].split(',')]), label=changedict[x], edgecolor='black'))
#        else:
#             patches5.append(patches.Patch(facecolor=tuple([int(x)/255 for x in class_colo_dict[x].split(',')]), label=str(x), edgecolor='black'))
#
#fig, ax = plt.subplots(figsize=(9,1.5))
#ax.set_axis_off()
#legend4=ax.legend(handles=patches5, loc=(0,0),ncol=5,fontsize=11, frameon=False, edgecolor='black',handlelength=1.5, handleheight=1.5)
#plt.show()
#fig.savefig("Circos_legend.png", dpi=1200, bbox_inches='tight')

#########################################################
def stacked_bar_chart(rank, pivoted_df, stack_vals, level_values_field, chart_title, x_label, y_label, filename, total_taxa_list, y_axis_limit):
    #
    # stacked_bar_chart: draws and saves a barchart figure to filename
    #
    # pivoted_df: dataframe which has been pivoted so columns correspond to the values to be plotted
    # stack_vals: the column names in pivoted_df to plot which will be list(dataframe_name.index)
    # level_values_field: column in the dataframe which has the values to be plotted along the x axis (typically time dimension) (.index)
    # chart_title: how to title chart
    # x_label: label for x axis
    # y_label: label for y axis
    # filename: full path filename to save file
    # color1: first color in spectrum for stacked bars
    # color2: last color in spectrum for stacked bars; routine will select colors from color1 to color2 evenly spaced
    #
    # Implementation: based on (http://randyzwitch.com/creating-stacked-bar-chart-seaborn/; https://gist.github.com/randyzwitch/b71d47e0d380a1a6bef9)
    # this routine draws overlapping rectangles, starting with a full bar reaching the highest point (sum of all values), and then the next shorter bar
    # and so on until the last bar is drawn.  These are drawn largest to smallest with overlap so the visual effect is that the last drawn bar is the
    # bottom of the stack and in effect the smallest rectangle drawn.
    #
    # Here "largest" and "smallest" refer to relationship to foreground, with largest in the back (and tallest) and smallest in front (and shortest).
    # This says nothing about which part of the bar appear large or small after overlap.
    #
    figsize_tuple=(15,8)
    #sns.set_style("ticks")
    plt.figure(figsize=figsize_tuple)
    #sns.set_context({"figure.figsize": figsize_tuple})
    sns.set_context("poster",font_scale=1)
    pivoted_df=pivoted_df.reset_index(drop=True)
    half_n_col=int(len(total_taxa_list)/2)+1
    colour_palette=(sns.xkcd_palette(["orange"])+sns.xkcd_palette(["pinkish red"])+sns.xkcd_palette(["clear blue"])+sns.xkcd_palette(["light green"])+sns.hls_palette(4, h=.8, s=1, l=.5)+sns.xkcd_palette(["light grey","yellow"])+sns.hls_palette(half_n_col+1, h=.75, s=.6, l=.7))
    handle_list=[]
    colour_count=0
    colour_reference_dictionary={}
    for x in total_taxa_list:
         colour_reference_dictionary[x]=colour_palette[colour_count]
         handle_list.append(patches.Patch(color=colour_reference_dictionary[x], label=str(x)))   
         colour_count+=1
    plt.clf()
    plt.gca().set_ylim(0,y_axis_limit)
    #
    stack_total_column = 'Stack_subtotal_xyz'  # placeholder name which should not exist in pivoted_df
    if stack_total_column in list(pivoted_df): #delete if there is carry over (from building script)
        del pivoted_df[stack_total_column]
    bar_num = 0
    legend_rectangles = []
    legend_names = []
    for bar_part in stack_vals:
        sub_count = 0
        pivoted_df[stack_total_column] = 0
        stack_value = ""
        for stack_value in stack_vals:
            stack_color = colour_reference_dictionary[stack_value] # for every item in the stack we create a new subset [stack_total_column] of 1 to N of the sub values
            pivoted_df[stack_total_column] += pivoted_df[stack_value]  # sum up total
            sub_count += 1
            if sub_count >= int(len(stack_vals)) - bar_num:  # we skip out after a certain number of stack values
                break
        # now we have set the subtotal and can plot the bar.  reminder: each bar is overalpped by smaller subsequent bars starting from y=0 axis
        bar_plot = sns.barplot(data=pivoted_df, x=pivoted_df.index,
                           y=stack_total_column, color=stack_color)
        legend_rectangles.append(plt.Rectangle((0,0),1,1,fc=stack_color, edgecolor = 'none'))  
        legend_names.append(stack_value)   # the "last" stack_value is the name of that part of the stack
        bar_num += 1
    if len(legend_rectangles)<=12:
        l = plt.legend(legend_rectangles, legend_names, loc=(.035,-.31), ncol = 6, prop={'size':13.5})
    if len(legend_rectangles)>12:
        l = plt.legend(legend_rectangles, legend_names, loc=(0,-.5), ncol = 4, prop={'size':13.5})
    l.draw_frame(False)
    bar_plot.set_xticklabels(level_values_field)
    for x in bar_plot.get_xticklabels():
        x.set_rotation(90)
    bar_plot.set(xlabel=x_label)
    plt.yticks(fontsize=20)
    bar_plot.set_ylabel(y_label, size=20)
    plt.tight_layout()
    plt.text(35, .25, chart_title, horizontalalignment='center',fontsize=25)
    #plt.title(chart_title, size=25, loc=(0.035, 0.5))
    sns.despine(right=True)
    fig=plt.gcf()
    fig.savefig(("OUT_Taxonomy_phylum_graphs_2/"+str(rank)+str(filename)), dpi=500)
    plt.show()
    time.sleep(2)
    plt.figure(figsize=(36,2))
    k=plt.legend(legend_rectangles, legend_names, ncol = 6, prop={'size':25})
    k.draw_frame(False)
    sns.despine(top=True, right=True, bottom=True, left=True)
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    fig=plt.gcf()
    fig.savefig(("OUT_Taxonomy_phylum_graphs_2/"+str(rank)+'_'+str(filename)+'_legend'), dpi=500)
    plt.show()
############################################################################
    
plt.clf() 
    
#print('Creating stacked barplots for CAZy abundance profile by taxonomic homology')
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_total_format, list(df_sum_activity_taxa_dict_total_format), list(df_sum_activity_taxa_dict_total_format.index), 'Total taxanomic CAZy domain distribution', 'CAZy class', 'Number of CAZy domain containing contigs', 'Total_taxonomic_cazy_distribution', total_relevant_taxa_list, 25)
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_w1_format, list(df_sum_activity_taxa_dict_w1_format), list(df_sum_activity_taxa_dict_w1_format.index), 'Week 1 taxanomic CAZy domain distribution', 'CAZy class', 'Number of CAZy domain containing contigs', 'Week_1_taxonomic_cazy_distribution', total_relevant_taxa_list, 15)
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_w3_format, list(df_sum_activity_taxa_dict_w3_format), list(df_sum_activity_taxa_dict_w3_format.index), 'Week 3 taxanomic CAZy domain distribution', 'CAZy class', 'Number of CAZy domain containing contigs', 'Week_3_taxonomic_cazy_distribution', total_relevant_taxa_list, 15)
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_w5_format, list(df_sum_activity_taxa_dict_w5_format), list(df_sum_activity_taxa_dict_w5_format.index), 'Week 5 taxanomic CAZy domain distribution', 'CAZy class', 'Number of CAZy domain containing contigs', 'Week_5_taxonomic_cazy_distribution', total_relevant_taxa_list, 15)
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_w10_format, list(df_sum_activity_taxa_dict_w10_format), list(df_sum_activity_taxa_dict_w10_format.index), 'Week 10 taxanomic CAZy domain distribution', 'CAZy class', 'Number of CAZy domain containing contigs', 'Week_10_taxonomic_cazy_distribution', total_relevant_taxa_list, 15)
#print('Creating stacked barplots for CAZy molar percentage profile by taxonomic homology')
###ascertain ylim
##Week 1 taxonomic CAZy domain distribution by ' + u"\u2211"r' $\bar{x}$'" molar percentage
#ylimit=max([df_sum_activity_taxa_dict_molpct_w1.max(axis=0).max(),df_sum_activity_taxa_dict_molpct_w3.max(axis=0).max(),df_sum_activity_taxa_dict_molpct_w5.max(axis=0).max(),df_sum_activity_taxa_dict_molpct_w10.max(axis=0).max()])
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_molpct_w1, list(df_sum_activity_taxa_dict_molpct_w1), list(df_sum_activity_taxa_dict_molpct_w1.index), "Week 1", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_1_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list, ylimit)
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_molpct_w3, list(df_sum_activity_taxa_dict_molpct_w3), list(df_sum_activity_taxa_dict_molpct_w3.index), "Week 3", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_3_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list, ylimit)
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_molpct_w5, list(df_sum_activity_taxa_dict_molpct_w5), list(df_sum_activity_taxa_dict_molpct_w5.index), "Week 5", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_5_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list, ylimit)
#stacked_bar_chart('phylum',df_sum_activity_taxa_dict_molpct_w10, list(df_sum_activity_taxa_dict_molpct_w10), list(df_sum_activity_taxa_dict_molpct_w10.index), "Week 10", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_10_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list, ylimit)
##_class
#ylimit=max([df_sum_activity_taxa_dict_molpct_w1_class.max(axis=0).max(),df_sum_activity_taxa_dict_molpct_w3_class.max(axis=0).max(),df_sum_activity_taxa_dict_molpct_w5_class.max(axis=0).max(),df_sum_activity_taxa_dict_molpct_w10_class.max(axis=0).max()])
#stacked_bar_chart('class',df_sum_activity_taxa_dict_total_format_class, list(df_sum_activity_taxa_dict_total_format_class), list(df_sum_activity_taxa_dict_total_format_class.index), 'Total taxanomic (class) CAZy distribution', 'CAZy class', 'Number of CAZy domain containing contigs', 'Total_taxonomic_cazy_distribution', total_relevant_taxa_list_class, 25)
#stacked_bar_chart('class',df_sum_activity_taxa_dict_molpct_w1_class, list(df_sum_activity_taxa_dict_molpct_w1_class), list(df_sum_activity_taxa_dict_molpct_w1_class.index), 'Week 1 taxonomic (class) CAZy domain distribution by ' + u"\u2211"r' $\bar{x}$'" molar percentage", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_1_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list_class, ylimit)
#stacked_bar_chart('class',df_sum_activity_taxa_dict_molpct_w3_class, list(df_sum_activity_taxa_dict_molpct_w3_class), list(df_sum_activity_taxa_dict_molpct_w3_class.index), 'Week 3 taxonomic (class) CAZy domain distribution by ' + u"\u2211"r' $\bar{x}$'" molar percentage", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_3_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list_class, ylimit)
#stacked_bar_chart('class',df_sum_activity_taxa_dict_molpct_w5_class, list(df_sum_activity_taxa_dict_molpct_w5_class), list(df_sum_activity_taxa_dict_molpct_w5_class.index), 'Week 5 taxonomic (class) CAZy domain distribution by ' + u"\u2211"r' $\bar{x}$'" molar percentage", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_5_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list_class, ylimit)
#stacked_bar_chart('class',df_sum_activity_taxa_dict_molpct_w10_class, list(df_sum_activity_taxa_dict_molpct_w10_class), list(df_sum_activity_taxa_dict_molpct_w10_class.index), 'Week 10 taxonomic CAZy (class) domain distribution by ' + u"\u2211"r' $\bar{x}$'" molar percentage", 'CAZy class', u"\u2211"r' $\bar{x}$'" molar percentage", 'Week_10_molpct_taxonomic_cazy_distribution', total_relevant_taxa_list_class, ylimit)



##############################################################
#########################################################################################
df_sum_activity=pd.DataFrame.from_dict(sum_activity_dict, orient='index') 
df_sum_CBM_activity=pd.DataFrame.from_dict(sum_activity_dict, orient='index')
#filter_list=['GT','dockerin','SLH', 'cohesin', 'CBM']
filter_list=['GT', 'CBM']
index_before_filter=len(df_sum_activity.index)
for x in filter_list:
     df_sum_activity=df_sum_activity[~(df_sum_activity.index.str.contains(x))]
print('Filtered ' + str(index_before_filter-(len(df_sum_activity))) +' keys contaning ' + str(filter_list) + ' from sum activity dataframe: ' + str((len(df_sum_activity))) +' keys remain')

df_sum_mean=pd.DataFrame(index=df_sum_activity.index)
df_sum_mean.index.names=['Molar percentage']

df_sum_mean['One']=(df_sum_activity['W1_C201'] + df_sum_activity['W1_3_C429'] + df_sum_activity['W1_4_C429'])/3
df_sum_mean['Three']=(df_sum_activity['W3_C201'] + df_sum_activity['W3_2_C429'] + df_sum_activity['W3_1_C429'])/3
df_sum_mean['Five']=(df_sum_activity['W5_C201'] + df_sum_activity['W5_25_C429'] + df_sum_activity['W5_100_C429'])/3
df_sum_mean['Ten']=(df_sum_activity['W10_C201'] + df_sum_activity['W10_100_i_C429'] + df_sum_activity['W10_100_ii_C429'])/3


filter_all_but_CBM=['GH', 'GT', 'dockerin', 'SLH', 'PL', 'cohesin', 'CE', 'AA']
index_before_filter=len(df_sum_CBM_activity.index)
for x in filter_all_but_CBM:
    df_sum_CBM_activity=df_sum_CBM_activity[~(df_sum_CBM_activity.index.str.contains(x))]
print('Filtered ' + str(index_before_filter-(len(df_sum_CBM_activity))) +' keys contaning ' + str(filter_list) + ' from sum activity dataframe (CBM only): ' + str((len(df_sum_CBM_activity))) +' keys remain')

df_sum_mean_CBM=pd.DataFrame(index=df_sum_CBM_activity.index)
df_sum_mean_CBM.index.names=['Molar percentage']

df_sum_mean_CBM['One']=(df_sum_CBM_activity['W1_C201'] + df_sum_CBM_activity['W1_3_C429'] + df_sum_CBM_activity['W1_4_C429'])/3
df_sum_mean_CBM['Three']=(df_sum_CBM_activity['W3_C201'] + df_sum_CBM_activity['W3_2_C429'] + df_sum_CBM_activity['W3_1_C429'])/3
df_sum_mean_CBM['Five']=(df_sum_CBM_activity['W5_C201'] + df_sum_CBM_activity['W5_25_C429'] + df_sum_CBM_activity['W5_100_C429'])/3
df_sum_mean_CBM['Ten']=(df_sum_CBM_activity['W10_C201'] + df_sum_CBM_activity['W10_100_i_C429'] + df_sum_CBM_activity['W10_100_ii_C429'])/3

##################DATA_FILTERING##############################################################################
#df_mean.index.str.contains('GT')] ===== this searches the index for regex GT, if found, returns TRUE        #
#the ~  symbol invert the boolean, so if found, returns FALSE                                                #
#then use that boolean index to filter the pandas dataframe with df_mean[bool]
#filter_list=['GT','SLH','cohesin','dockerin']   
#for x in filter_list:
#     df_mean=df_mean[~df_mean.index.str.contains(x)]
#num_filter=len(filter_list)-2                                                               
###############################################################################################################
print('\n\nUnique ORF matching peptides in week 1: ' + str(sum(df_mean['One']>0)) +'\nUnique ORF matching peptides in week 3: ' + str(sum(df_mean['Three']>0)) +'\nUnique ORF matching peptides in week 5: ' + str(sum(df_mean['Five']>0)) +'\nUnique ORF matching peptides in week 10: ' + str(sum(df_mean['Ten']>0)))
#############################HEATMAP#################################
ticklist=list(x/20 for x in list(range(0,50)))
ticklabel=list(str(x)+'%' for x in ticklist)
figsize=(13.5,24)
fig, ax = plt.subplots(figsize=figsize)
#sns.heatmap(data=df_mean, cmap='Blues')
ax=sns.heatmap(data=df_mean, cmap='Blues')
ax.set(xlabel = 'Time point', ylabel='Activity(ORF id)')
cbar = ax.collections[0].colorbar
cbar.set_ticks(ticklist)
cbar.set_ticklabels(ticklabel)
cbar.set_label(r'$\bar{x}$''(Molarity)')
#Fig = ax.get_figure()
#fig.savefig('heatmap_composite_unfiltered_blocks.png')
####################################################################
##The GH6 ORF is way too high in comparison, so filter it out, it is the max value in week one, so to return the index (ORF name)
##(df_mean['W1'].argmax())return a_c570829...etc
##Use the index name to return the index integer using df.index.get_loc()
##df.index.get_loc((df_mean['W1'].argmax()))
##Use this to filter from the df into a new df using df_mean.drop(df.index[the index generated from the above code here])
#df_filtered=df_mean.drop(df.index[df.index.get_loc((df_mean['W1'].argmax()))])
#
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')   
#######Create a separate dataframe to harbour values from different classes of enzymes, to assing row colour too############
df_col_row=pd.DataFrame(index=df_mean.index)
#create empty column
df_col_row['color']=np.nan
df_list=[]
num_colors=11
#col_row_pallette = sns.cubehelix_palette(n_colors=num_colors,light=.9, dark=.1, start=1, rot=-2)
#col_row_pallette = sns.color_palette("RdBu_r", num_colors)
#col_row_pallette=sns.diverging_palette(150, 275, s=80, l=55, n=num_colors)
col_row_pallette=sns.diverging_palette(n=num_colors, h_neg=210, h_pos=350, s=90, l=30)
#col_row_palette is a list of lists, the sub lists contain rgb numbers...
#Create a list of rgb codes for each number assigned to CAZy blocks above
row_indices=pd.Series(index=range(0,num_colors))
#Create an empty series with indexes that match the amount of colours needed
index_color_dict = dict(zip(map(str, row_indices.index), col_row_pallette))
legend_dict={}
#to insert value into dataframe syntax = df_col_row[index_integer, column]=value
df_list=[]
for x in df_col_row.index:
    if re.search(r'CBM', x) and not re.search(r'GH', x) and not re.search(r'CE', x) and not re.search(r'PL', x) and not re.search('AA', x):
        #print('CBM/unknown ' + x)
        df_col_row.set_value(x,'color',0)
        df_list.append(col_row_pallette[0])
        legend_dict['CBM/unknown']=col_row_pallette[0]
    elif re.search(r'CBM', x) and re.search(r'GH', x) and not re.search(r'CE', x):
        #print('CBM/GH ' + x)
        df_col_row.set_value(x,'color',1)
        df_list.append(col_row_pallette[1])
        legend_dict['CBM/GH']=col_row_pallette[1]
    elif re.search(r'CBM', x) and re.search(r'GH', x) and re.search(r'CE', x):    
        #print('CBM/GH/CE ' + x)
        df_col_row.set_value(x,'color',2)
        df_list.append(col_row_pallette[2])
        legend_dict['CBM/GH/CE']=col_row_pallette[2]
    elif re.search(r'GH', x):
        #print('GH ' + x)
        df_col_row.set_value(x,'color',3)
        df_list.append(col_row_pallette[3])
        legend_dict['GH']=col_row_pallette[3]
    elif re.search(r'AA', x):
        #print('AA ' + x)
        df_col_row.set_value(x,'color',4)
        df_list.append(col_row_pallette[4])
        legend_dict['AA']=col_row_pallette[4]
    elif re.search(r'CBM', x) and re.search(r'CE', x):
       # print('CBM/CE ' + x)
        df_col_row.set_value(x,'color',5)
        df_list.append(col_row_pallette[5])
        legend_dict['CBM/CE']=col_row_pallette[5]
    elif re.search(r'CE', x):
        #print('CE ' + x)
        df_col_row.set_value(x,'color',6)
        df_list.append(col_row_pallette[6])
        legend_dict['CE']=col_row_pallette[6]
    elif re.search(r'PL', x):
        #print('PL ' + x)
        df_col_row.set_value(x,'color',7)
        df_list.append(col_row_pallette[7])
        legend_dict['PL']=col_row_pallette[7]
    elif re.search(r'GT', x):
        #print('GT ' + x)
        df_col_row.set_value(x,'color',8)
        df_list.append(col_row_pallette[8])
        legend_dict['GT']=col_row_pallette[8]
    elif re.search(r'dockerin', x) or re.search(r'cohesin', x) or re.search(r'SLH', x):    
        #print('Cellulosome ' + x)
        df_col_row.set_value(x,'color',9)
        df_list.append(col_row_pallette[9])
        legend_dict['Cellulosome associated']=col_row_pallette[9]
    else: 
        #print('Other ' + x)
        df_col_row.set_value(x,'color',10)
        df_list.append(col_row_pallette[10])
        legend_dict['Other']=col_row_pallette[10]

#########Domain_level_colours##########
#filter taxa dict to contain only the ORFs in the current dataframe (df_mean etc)
taxa_master_dict_named_filtered={}
for key in taxa_master_dict_named:
    if key in list(df_mean.index):
        taxa_master_dict_named_filtered[key]=taxa_master_dict_named[key]
        
number_domains=4
domain_col_pal=sns.color_palette("husl", number_domains)
#domain_col_pal=sns.cubehelix_palette(number_domains,light=.9, dark=.5, reverse=True, start=.5, rot=-2)
domain_dict={}
domain_list=[]
for x in taxa_master_dict_named_filtered:
    if taxa_master_dict_named_filtered[x]['Taxonomy']['superkingdom']=='Bacteria':
        domain_dict['Bacteria']=domain_col_pal[0]
        domain_list.append(domain_col_pal[0])
    elif taxa_master_dict_named_filtered[x]['Taxonomy']['superkingdom']=='Eukaryota':
        domain_dict['Eukaryota']=domain_col_pal[1]
        domain_list.append(domain_col_pal[1])
    elif taxa_master_dict_named_filtered[x]['Taxonomy']['superkingdom']=='Archaea':
        domain_dict['Archaea']=domain_col_pal[2]
        domain_list.append(domain_col_pal[2])
    elif taxa_master_dict_named_filtered[x]['Taxonomy']['superkingdom']=='NA':
        domain_dict['No homolog']=domain_col_pal[3]
        domain_list.append(domain_col_pal[3])
        
#############Class level colours################
#Get number of classes
num_king_phyl_list=[]
for x in taxa_master_dict_named_filtered:
    if taxa_master_dict_named_filtered[x]['Taxonomy']['phylum'] not in num_king_phyl_list:
        num_king_phyl_list.append(taxa_master_dict_named_filtered[x]['Taxonomy']['phylum'])
        
num_king_phyl=len(num_king_phyl_list)

king_phyl_dict={}
king_phyl_list=[]
king_phyl_col_pal=sns.color_palette("husl", num_king_phyl)
for x in taxa_master_dict_named_filtered:
    for classx in num_king_phyl_list:
        if taxa_master_dict_named_filtered[x]['Taxonomy']['phylum']==classx:
              king_phyl_dict[classx]=king_phyl_col_pal[(num_king_phyl_list.index(classx))]
              king_phyl_list.append(king_phyl_col_pal[(num_king_phyl_list.index(classx))])
             
   
    #####Cluster_map####
#normalise across rows - i.e. per peptide matching ORFs use z_score=0 or 1
#normalise across columns - i.e per timepoint use sandard_scale=1

#df_mean.values returns ROW-wise array
#df_mean.values.T returns COL-wise array
    
#The colour dictionaries are uses to create the legends, the colour lists are used to create the row colours
df_colors=pd.Series(df_list, index=df_mean.index, name='Group')
df_domain_colours=pd.Series(domain_list, index=df_mean.index, name='Taxonomy')

df_phylum_colors=pd.Series(king_phyl_list, index=df_mean.index, name='Taxonomy')

group_domain_colors_dict = pd.DataFrame(dict(Group=df_colors, Domain=df_domain_colours))    
group_phylum_colors_dict = pd.DataFrame(dict(Group=df_colors, Phylum=df_phylum_colors))

legend_handles_group=[]
for x in legend_dict:
     legend_handles_group.append(patches.Patch(color=legend_dict[x], label=x))
     
legend_handles_domain=[]
for x in domain_dict:
     legend_handles_domain.append(patches.Patch(color=domain_dict[x], label=x))
     
legend_handles_phylum=[]
for x in king_phyl_dict:
     legend_handles_phylum.append(patches.Patch(color=king_phyl_dict[x], label=x))    
####Overall clustmap with phylum #####
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_mean.values, df_mean.values.T))
k = sns.clustermap(df_mean, row_linkage=row_linkage, col_linkage=col_linkage,col_cluster=False, figsize=(9,28), standard_scale=0, cmap='Blues', row_colors=group_phylum_colors_dict)
plt.setp(k.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
k.cax.set_position([1.3, .2, .03, .45]) #([x, y, width, height])
legend2=k.ax_heatmap.legend(loc=(0.05,1), handles=legend_handles_phylum, title='Taxonomic homology', ncol=5)
k.ax_heatmap.add_artist(legend2)
k.ax_heatmap.legend(loc=(0.35,1.05),handles=legend_handles_group, title='Group', ncol=4)    
plt.show()
plt.clf()
 
####overall clustermap with domain#####
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_mean.values, df_mean.values.T))
kd = sns.clustermap(df_mean, row_linkage=row_linkage, col_linkage=col_linkage,col_cluster=False, figsize=(9,28), standard_scale=0, cmap='Blues', row_colors=group_domain_colors_dict)
plt.setp(kd.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
kd.cax.set_position([1.3, .2, .03, .45]) #([x, y, width, height])
legend2=kd.ax_heatmap.legend(loc=(0.05,1), handles=legend_handles_domain, title='Taxonomic homology', ncol=5)
kd.ax_heatmap.add_artist(legend2)
kd.ax_heatmap.legend(loc=(0.35,1.05),handles=legend_handles_group, title='Group', ncol=4)    
plt.show()
plt.clf()

#k.savefig('clustermap_GT_cell_filtered_block_derep_row_normalised.png')
######sum clustermap#########
df_sum_mean_col_dict={'AA':0, 'CE':1, 'GH':2, 'PL':3, 'SLH':4, 'cohesin':5, 'dockerin':6}

#df_sum_mean_class_colors=sns.hls_palette(7)#, l=.4, s=.7)
    
df_sum_mean_class_colors=sns.xkcd_palette(['pale orange', 'shocking pink', "robin's egg", 'pastel yellow', 'vivid green', 'pale lilac', 'light green'])    
    
df_sum_mean_col_list=[]
for x in df_sum_mean.index:
    for y in df_sum_mean_col_dict:
         if bool(re.search(y, x))==True:
            df_sum_mean_col_list.append(df_sum_mean_class_colors[df_sum_mean_col_dict[y]])
df_sum_mean_col_df=pd.Series(df_sum_mean_col_list, index=df_sum_mean.index, name='')
###
enzyme_class_legend_handles=[]
for x in df_sum_mean_col_dict:
    enzyme_class_legend_handles.append(patches.Patch(color=df_sum_mean_class_colors[df_sum_mean_col_dict[x]], label=x))

sns.set(font_scale=1.5)
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_sum_mean.iloc[:,:4].values, df_sum_mean.iloc[:,:4].values.T))
tc = sns.clustermap(df_sum_mean[['One', 'Three', 'Five', 'Ten']], row_linkage=row_linkage, col_linkage=col_linkage,row_colors=df_sum_mean_col_df, figsize=(7,12), standard_scale=1, cmap='Blues')
plt.setp(tc.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=12)
leg=tc.ax_heatmap.legend(loc=(-0.155,-0.2), handles=enzyme_class_legend_handles, title='Enzyme class', ncol=4, frameon=False)
plt.setp(leg.get_title(),fontsize='16')
tc.cax.set_position([1.05, .2, .03, .45])
plt.tick_params(which='y', size=10)
col = tc.ax_col_dendrogram.get_position()
row = tc.ax_row_dendrogram.get_position()
tc.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*0.25])
tc.ax_row_dendrogram.set_position([row.x0*-.005, row.y0, row.width*2, row.height])
tc.ax_heatmap.set(xlabel = 'Week', ylabel='Normalised 'u"\u2211"r' $\bar{x}$'" molar percentage")
clustermapstr='OUT_Clustermaps/CAzy_clustermap_by_class_PUB_QUAL_thin'+time.strftime('%H%M%S') + '.png'
tc.savefig(clustermapstr, dpi=1200)

plt.clf()
  
################PAPER VERSION#######################################
################PAPER VERSION#######################################

df_sum_mean_class_colors=sns.xkcd_palette(['pale orange', 'shocking pink', "robin's egg", 'pastel yellow', 'vivid green', 'pale lilac', 'green'])    
enzyme_class_legend_handles2=[]
for x in df_sum_mean_col_dict:
    if not x == 'cohesin' and not x == 'dockerin':
         enzyme_class_legend_handles2.append(patches.Patch(color=df_sum_mean_class_colors[df_sum_mean_col_dict[x]], label=x))

sns.set(font_scale=0.5)
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_sum_mean.iloc[:,:4].values, df_sum_mean.iloc[:,:4].values.T))
tc = sns.clustermap(df_sum_mean[['One', 'Three', 'Five', 'Ten']], row_linkage=row_linkage, col_linkage=col_linkage,row_colors=df_sum_mean_col_df, figsize=(3,10), cmap='Blues', col_cluster=False,linecolor='white', linewidths=1, annot=False, annot_kws={"size": 12}, fmt='.1g')
plt.setp(tc.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=10)
plt.setp(tc.ax_heatmap.xaxis.get_majorticklabels(), rotation=0, size=14)
tc.ax_heatmap.tick_params(axis='y', which='both', labelsize=11)
leg=tc.ax_heatmap.legend(loc=(1.38,-0.025), handles=enzyme_class_legend_handles2, title='CAZy\nclass', ncol=1, frameon=False, fontsize=10, handlelength=0.75, handleheight=0.5)
plt.setp(leg.get_title(),fontsize='13')
tc.ax_heatmap.tick_params(axis='y', which='both', length=0)
tc.ax_heatmap.set_xticklabels(['1','3','5','10'], size=14)
tc.ax_heatmap.set_ylabel('', fontsize=0)
plt.text(-3.25, 0.55, ''u"\u2211"r' $\bar{x}$'" mol%", fontsize=14, rotation=90)
tc.ax_heatmap.set_xlabel('Week', fontsize=14)
tc.cax.tick_params(labelsize=12)
tc.cax.set_position([1.18, .275, .03, .45])
col = tc.ax_col_dendrogram.get_position()
row = tc.ax_row_dendrogram.get_position()
tc.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*.25])
tc.ax_row_dendrogram.set_position([row.x0, row.y0, row.width, row.height])
#tc.savefig('CAZy_clustermap_figsie_5_10.png', dpi=1500)
plt.clf()

####################### CBM HEATMAP #####################################
cell_color=(0.98, 0.37, 0.96)
hemicell_color=(0, 0.8, 0.996)
chitin_color=(0.5, 0.97, 0.67)
    
cellulose_associated_colours=[(1,1,1),cell_color]
hemicellulose_associated_colours=[(1,1,1),hemicell_color]
chitin_associated_colours=[(1,1,1),chitin_color]

cbm_cell_list=[]
cbm_hemi_list=[]
cbm_chitin_list=[]              
              
for x in df_sum_mean_CBM.index:
     if cbm_pickle_cellulose[x]==1:
         cbm_cell_list.append(cellulose_associated_colours[1])
     if cbm_pickle_cellulose[x]==0:
         cbm_cell_list.append(cellulose_associated_colours[0])
     if cbm_pickle_hemicellulose[x]==1:
         cbm_hemi_list.append(hemicellulose_associated_colours[1])
     if cbm_pickle_hemicellulose[x]==0:
         cbm_hemi_list.append(hemicellulose_associated_colours[0])
     if cbm_pickle_chitin[x]==1:
         cbm_chitin_list.append(chitin_associated_colours[1])
     if cbm_pickle_chitin[x]==0:
         cbm_chitin_list.append(chitin_associated_colours[0])
 
cell_row_series=pd.Series(cbm_cell_list, index=df_sum_mean_CBM.index, name='') #Cellulose
hemicell_row_series=pd.Series(cbm_hemi_list, index=df_sum_mean_CBM.index, name='') #Hemicellulose
chitin_row_series=pd.Series(cbm_chitin_list, index=df_sum_mean_CBM.index, name='') #Chitin
  
#CBM_master_row_col_dict=pd.DataFrame(dict(Cellulose=cell_row_series, Hemicellulose=hemicell_row_series, Chitin=chitin_row_series))

CBM_master_row_col_handles=[]
CBM_master_row_col_handles.append(patches.Patch(color=cell_color, label='Cellulose'))
CBM_master_row_col_handles.append(patches.Patch(color=hemicell_color, label='Hemicellulose'))
CBM_master_row_col_handles.append(patches.Patch(color=chitin_color, label='Chitin'))

#sns.set(font_scale=1.5)
#row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_sum_mean_CBM.values[:,:4], df_sum_mean_CBM.values.T[:,:4]))
#cbm = sns.clustermap(df_sum_mean_CBM.ix[:,:4], row_linkage=row_linkage, col_linkage=col_linkage,col_cluster=False, figsize=(7,12), standard_scale=1, cmap='Blues', row_colors=[cell_row_series, hemicell_row_series, chitin_row_series])
#plt.setp(cbm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=12)
#cbm.cax.set_position([1.05, .2, .03, .45])
#cbm.ax_heatmap.set(xlabel = 'Week', ylabel='Normalised 'u"\u2211"r' $\bar{x}$'" molar percentage")
#leg=cbm.ax_heatmap.legend(loc=(0,1), handles=CBM_master_row_col_handles, title='Known binding associations', ncol=5, fontsize=12)
#plt.setp(leg.get_title(),fontsize='16')
#plt.tick_params(which='y', size=10)

sns.set(font_scale=1.3)
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_sum_mean_CBM.iloc[:,:4].values, df_sum_mean_CBM.iloc[:,:4].values.T))
cbm = sns.clustermap(df_sum_mean_CBM.iloc[:,:4], standard_scale=1, row_linkage=row_linkage, col_linkage=col_linkage,col_cluster=False, figsize=(6,10), cmap='Blues', row_colors=[cell_row_series, hemicell_row_series, chitin_row_series],  linecolor='white', linewidths=1, annot=False, annot_kws={"size": 16}, fmt='.1g')
plt.setp(cbm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
cbm.cax.set_position([1.1, .2, .03, .45])
cbm.ax_heatmap.set(xlabel = 'Week', ylabel='Normalised 'u"\u2211"r' $\bar{x}$'" mol%")
leg=cbm.ax_heatmap.legend(loc=(-0.1,1), handles=CBM_master_row_col_handles, title='Known binding associations', ncol=5, frameon=False)
row = cbm.ax_row_dendrogram.get_position()
cbm.ax_row_dendrogram.set_position([row.x0, row.y0, row.width, row.height])
hmp=cbm.ax_heatmap.get_position()
cbm.ax_heatmap.set_position([hmp.x0, hmp.y0, hmp.width, hmp.height])
plt.setp(leg.get_title(),fontsize='16')
plt.tick_params(which='y', size=10)
plt.clf()


#NON NORMALISED PAPER VERSION###############
sns.set(font_scale=1.2)
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_sum_mean_CBM.iloc[:,:4].values, df_sum_mean_CBM.iloc[:,:4].values.T))
cbm = sns.clustermap(df_sum_mean_CBM.iloc[:,:4], row_linkage=row_linkage, col_linkage=col_linkage,col_cluster=False, figsize=(3,13), cmap='Blues', linecolor='white', linewidths=1, annot=False, annot_kws={"size": 12}, fmt='.1g')
plt.setp(cbm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=14)
plt.setp(cbm.ax_heatmap.xaxis.get_majorticklabels(), rotation=0, size=14)
leg=cbm.ax_heatmap.legend(loc=(1.38,-0.025), handles=CBM_master_row_col_handles, title='   Known\n   binding\nassociations', ncol=1, frameon=False, fontsize=14, handlelength=0.35, handleheight=2)
cbm.ax_heatmap.tick_params(axis='y', which='both', length=0)
cbm.ax_heatmap.tick_params(axis='y', which='both', labelsize=14)
cbm.ax_heatmap.set_ylabel('', fontsize=18)
plt.text(-4.5, 0.55, ''u"\u2211"r' $\bar{x}$'" mol%", fontsize=18, rotation=90)
cbm.ax_heatmap.set_xlabel('Week', fontsize=18)
cbm.cax.tick_params(labelsize=16)
cbm.ax_heatmap.set_xticklabels(['1','3','5','10'], size=17)
plt.setp(leg.get_title(),fontsize='16')
cbm.cax.set_position([1.275, .275, .03, .45])

#draw_col_colors with decent thickness
#This returns the reordered row indexes in order they were plotted
reordered_index=cbm.dendrogram_row.reordered_ind
yval=0.05
offset=0.05
for index in reordered_index:
     key=df_sum_mean_CBM.iloc[:,:4].iloc[index].name
     cbm.ax_heatmap.add_patch(plt.Rectangle((-0.1,yval),0.1,0.9,facecolor=cell_row_series[key],clip_on=False,linewidth = 0))
     cbm.ax_heatmap.add_patch(plt.Rectangle((-0.2-offset,yval),0.1,0.9,facecolor=hemicell_row_series[key],clip_on=False,linewidth = 0))
     cbm.ax_heatmap.add_patch(plt.Rectangle((-0.35-offset,yval),0.1,0.9,facecolor=chitin_row_series[key],clip_on=False,linewidth = 0))
     yval+=1
#set_dendro_width
row = cbm.ax_row_dendrogram.get_position()
col = cbm.ax_col_dendrogram.get_position()
cbm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*0.33])
cbm.ax_row_dendrogram.set_position([row.x0-.08, row.y0, row.width, row.height])
hmp=cbm.ax_heatmap.get_position()
cbm.ax_heatmap.set_position([hmp.x0, hmp.y0, hmp.width, hmp.height])
plt.tick_params(which='y', size=10)
#cbm.savefig('CBM_heatmap_non_normalised_paper_version.png', dpi=1500)

cbmclustermapstr='OUT_Clustermaps/CBM_clustermap_PUB_QUAL'+time.strftime('%H%M%S') + '.png' 
cbm.savefig("CBM_clustermap_manual_row_col_3_13.png", dpi=1200, bbox_inches='tight')
#cbm.savefig(cbmclustermapstr, dpi=1200)
plt.clf()
 
############################CBM_BARGRAPHS###############################################

df_sum_mean_CBM['One_err']=df_sum_CBM_activity[['W1_C201','W1_3_C429', 'W1_4_C429']].std(axis=1)
df_sum_mean_CBM['Three_err']=df_sum_CBM_activity[['W3_C201','W3_2_C429','W3_1_C429']].std(axis=1)
df_sum_mean_CBM['Five_err']=df_sum_CBM_activity[['W5_C201','W5_25_C429','W5_100_C429']].std(axis=1)
df_sum_mean_CBM['Ten_err']=df_sum_CBM_activity[['W10_C201','W10_100_i_C429','W10_100_ii_C429']].std(axis=1)

df_sum_mean_CBM['sum']=df_sum_mean_CBM['One']+df_sum_mean_CBM['Three']+df_sum_mean_CBM['Five']+df_sum_mean_CBM['Ten']
df_sum_mean_CBM['class_rank']=df_sum_mean_CBM['sum'].rank(ascending=False)

top_x_CBM_for_barplots=15
print('Plotting mean molar percentage for top ' + str(top_x_CBM_for_barplots) + ' CBMs')
CBM_for_barplots_list=[]
CBM_for_barplots_list=df_sum_mean_CBM[df_sum_mean_CBM['class_rank']<=top_x_CBM_for_barplots].sort_values(by=['class_rank']).index.tolist()
f, axarr = plt.subplots(3, 5, figsize=(15, 9))
sns.set(style="white")
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
colour_barplot=sns.xkcd_palette(['clear blue'])
#axarr[row, col]
for CBM in list(range(top_x_CBM_for_barplots)):
     xvals=list(df_sum_mean_CBM.loc[CBM_for_barplots_list[CBM]][:4])
     width=0.65
     ind=[0, 1, 2, 3]
     yerrx=list(df_sum_mean_CBM.loc[CBM_for_barplots_list[CBM]][4:8])
     if CBM<=4:
          axarr[0, CBM].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
          axarr[0, CBM].set_title(CBM_for_barplots_list[CBM], size=16)
          axarr[0, CBM].set_ylim([0, 0.4])
          axarr[0, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
     elif CBM<=9:
          axarr[1, CBM-5].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
          axarr[1, CBM-5].set_title(CBM_for_barplots_list[CBM], size=16)
          axarr[1, CBM-5].set_ylim([0, 0.4])
          axarr[1, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
     elif CBM<=14:
          axarr[2, CBM-10].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
          axarr[2, CBM-10].set_title(CBM_for_barplots_list[CBM], size=16)
          axarr[2, CBM-10].set_ylim([0, 0.4])
          axarr[2, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
          axarr[2, CBM-10].set_xlabel('Week', fontsize=17)
          axarr[2, CBM-10].set_xticks([0, 1, 2, 3])          
          axarr[2, CBM-10].set_xticklabels(['One', 'Three', 'Five', 'Ten'], fontsize=16)
plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 3]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 4]], visible=False)
sns.despine(top=True, right=True)
f.savefig('OUT_grids/top_15_CBMs.png', dpi=500)
plt.show()
plt.clf()
plt.close(f)  
##########################################################################################
taxonomy_sum_dict={}
taxa_list=[]
for k in df_mean.index:
     dictname_one_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['phylum']), 'One')
     dictname_three_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['phylum']), 'Three')
     dictname_five_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['phylum']),'Five')
     dictname_ten_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['phylum']),'Ten')
     dict_name_list=[dictname_one_tuple,dictname_three_tuple,dictname_five_tuple,dictname_ten_tuple]
     for key in dict_name_list:
          if key not in taxonomy_sum_dict:
             if taxa_master_dict_named[k]['Taxonomy']['phylum'] not in taxa_list:
                 taxa_list.append(taxa_master_dict_named[k]['Taxonomy']['phylum'])
             taxonomy_sum_dict[key]=(df_mean.loc[k][key[1]])
          if key in taxonomy_sum_dict:
              taxonomy_sum_dict[key]=(df_mean.loc[k][key[1]])+taxonomy_sum_dict[key]


taxonomy_sum_dictclass={}
taxa_list_class=[]
for k in df_mean.index:
     dictname_one_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']), 'One')
     dictname_three_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']), 'Three')
     dictname_five_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']),'Five')
     dictname_ten_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']),'Ten')
     dict_name_list=[dictname_one_tuple,dictname_three_tuple,dictname_five_tuple,dictname_ten_tuple]
     for key in dict_name_list:
          if key not in taxonomy_sum_dictclass:
             if taxa_master_dict_named[k]['Taxonomy']['class'] not in taxa_list_class:
                 taxa_list_class.append(taxa_master_dict_named[k]['Taxonomy']['class'])
             taxonomy_sum_dictclass[key]=(df_mean.loc[k][key[1]])
          if key in taxonomy_sum_dictclass:
              taxonomy_sum_dictclass[key]=(df_mean.loc[k][key[1]])+taxonomy_sum_dictclass[key]


##################ADD STANDARD DEVIATION TO COLUMNS IN DF MEAN FOR GRAPHS ETC#################
df_mean['One_err']=df[['W1_C201','W1_3_C429', 'W1_4_C429']].std(axis=1)
df_mean['Three_err']=df[['W3_C201','W3_2_C429','W3_1_C429']].std(axis=1)
df_mean['Five_err']=df[['W5_C201','W5_25_C429','W5_100_C429']].std(axis=1)
df_mean['Ten_err']=df[['W10_C201','W10_100_i_C429','W10_100_ii_C429']].std(axis=1)
##############################################################################################

#NB/ samplenamelist can be used if needed to split into replicates

week_list=['One','Three','Five','Ten']
taxa_index=pd.MultiIndex.from_product([taxa_list, week_list])
df2=pd.DataFrame(index=taxa_index)

df_taxa_sum=pd.DataFrame.from_dict(taxonomy_sum_dict, orient='index')
df2=pd.concat([df2, df_taxa_sum], axis=1, join='inner')
#df2.rename(columns={0:'mean_value'}) ----this is superfluous
df2.reset_index(drop=False, inplace=True)
df2=df2.rename(columns={'level_0': 'Taxonomy', 'level_1':'Week', 0:'Molar_percentage'})
#create an error columns of zeros, as we cannot generate error bars for this data...
df2['err']=pd.Series(np.zeros(len(df2.index)), index=df2.index)
#Sort the week column so that they display in the correct order in the graphs...
custom_sort_dict={'One':1,'Three':3,'Five':5,'Ten':10}
df2['temp_rank']=df2['Week'].map(custom_sort_dict)
df2.sort_values(by=['temp_rank'],inplace=True)
df2.drop(labels=['temp_rank'],axis=1)
sns.set(style="white")





week_list=['One','Three','Five','Ten']
taxa_index_class=pd.MultiIndex.from_product([taxa_list_class, week_list])
dfclass=pd.DataFrame(index=taxa_index_class)

df_taxa_sum_class=pd.DataFrame.from_dict(taxonomy_sum_dictclass, orient='index')
dfclass=pd.concat([dfclass, df_taxa_sum_class], axis=1, join='inner')
#df2.rename(columns={0:'mean_value'}) ----this is superfluous
dfclass.reset_index(drop=False, inplace=True)
dfclass=dfclass.rename(columns={'level_0': 'Taxonomy', 'level_1':'Week', 0:'Molar_percentage'})
#create an error columns of zeros, as we cannot generate error bars for this data...
dfclass['err']=pd.Series(np.zeros(len(dfclass.index)), index=dfclass.index)
#Sort the week column so that they display in the correct order in the graphs...
custom_sort_dict={'One':1,'Three':3,'Five':5,'Ten':10}
dfclass['temp_rank']=dfclass['Week'].map(custom_sort_dict)
dfclass.sort_values(by=['temp_rank'],inplace=True)
dfclass.drop(labels=['temp_rank'],axis=1)
sns.set(style="white")




print('Responsible taxa based on sequence homology at the phylum level, less than ' + str(phylum_percent_filter) + ' filtered')
 
df2_filter=df2.groupby(['Taxonomy'])['Molar_percentage'].sum().reset_index()
df2_filter['contribution']=df2_filter['Molar_percentage']/(df2_filter['Molar_percentage'].sum())*100
df2_filter.set_index('Taxonomy')
df2_filter[df2_filter['contribution']>phylum_percent_filter]['Taxonomy']
df2_filter=df2[df2['Taxonomy'].isin(df2_filter[df2_filter['contribution']>phylum_percent_filter]['Taxonomy'])]
####Just the dominant phyla#########
sns.set(style='white', font_scale=1.6)
t=sns.factorplot(data=df2, x="Week", y="Molar_percentage", col="Taxonomy", kind="bar", col_wrap=5, palette=sns.husl_palette(6, h=.5))
t.set_yticklabels(fontsize=20)
t.set_xticklabels(fontsize=20)
t.set_ylabels(u"\u2211"' ' r'$\bar{x}$'" Molar percentage", fontsize=20)
t.set_xlabels(fontsize=20)
t.savefig("OUT_taxonomy_phylum_graphs/Top_phylum_ALL_unfiltered_dpi_500.png", dpi=500)
plt.show(t)
plt.clf()

print('Top three responsible phyla only') 

tax_filt=set(list(df2['Taxonomy']))
tax_retain=['Bacteroidetes','Verrucomicrobia','Proteobacteria']
for x in tax_retain:
    tax_filt.remove(x)
for x in tax_filt:
    df2=df2[~df2['Taxonomy'].str.contains(x)]
t_f=sns.factorplot(data=df2, x="Week", y="Molar_percentage", col="Taxonomy", kind="bar", col_wrap=3, palette=sns.husl_palette(6, h=.5))
t_f.set_ylabels(u"\u2211""("r'$\bar{x}$'"(""Molar percentage""))")
plt.show(t_f)
plt.clf()

####################class level taxa################

taxonomy_sum_dict_class={}
taxa_list_class=[]
for k in df_mean.index:
      dictname_one_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']), 'One')
      dictname_three_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']), 'Three')
      dictname_five_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']),'Five')
      dictname_ten_tuple=(str(taxa_master_dict_named[k]['Taxonomy']['class']),'Ten')
      dict_name_list=[dictname_one_tuple,dictname_three_tuple,dictname_five_tuple,dictname_ten_tuple]
      for key in dict_name_list:
          if key not in taxonomy_sum_dict_class:
              if taxa_master_dict_named[k]['Taxonomy']['class'] not in taxa_list_class:
                  taxa_list_class.append(taxa_master_dict_named[k]['Taxonomy']['class'])
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])
          if key in taxonomy_sum_dict_class:
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])+taxonomy_sum_dict_class[key]

week_list=['One','Three','Five','Ten']
taxa_index=pd.MultiIndex.from_product([taxa_list_class, week_list])
df3=pd.DataFrame(index=taxa_index)

df_taxa_sum=pd.DataFrame.from_dict(taxonomy_sum_dict_class, orient='index')
df3=pd.concat([df3, df_taxa_sum], axis=1, join='inner')
df3.reset_index(drop=False, inplace=True)
df3=df3.rename(columns={'level_0': 'Taxonomy', 'level_1':'Week', 0:'Molar_percentage'})
#create an error columns of zeros, as we cannot generate error bars for this data...
df3['err']=pd.Series(np.zeros(len(df3.index)), index=df3.index)
#Sort the week column so that they display in the correct order in the graphs...
custom_sort_dict={'One':1,'Three':3,'Five':5,'Ten':10}
df3['temp_rank']=df3['Week'].map(custom_sort_dict)
df3.sort_values(by=['temp_rank'],inplace=True)
df3.drop(labels=['temp_rank'],axis=1)
sns.set(style="white")

###########################TAXA clustermap
df_class_taxa_heatmap=df3.pivot(index='Taxonomy', columns='Week', values='Molar_percentage').reset_index().set_index('Taxonomy')
df_class_taxa_heatmap=df_class_taxa_heatmap.reindex(df_class_taxa_heatmap.sum(axis=1).rank(ascending=False).sort_values().index)
df_class_taxa_heatmap=df_class_taxa_heatmap[['One', 'Three', 'Five', 'Ten']]

sns.set(font_scale=1.6)
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_class_taxa_heatmap.values, df_class_taxa_heatmap.values.T))
tchm = sns.clustermap(df_class_taxa_heatmap, row_linkage=row_linkage, col_linkage=col_linkage, figsize=(7,12), cmap='viridis_r')
plt.setp(tchm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=15)
tchm.cax.set_position([1.275, .2, .03, .45])
plt.tick_params(which='y', size=14)
col = tchm.ax_col_dendrogram.get_position()
row = tchm.ax_row_dendrogram.get_position()
tchm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*0.25])
tchm.ax_row_dendrogram.set_position([row.x0*-.005, row.y0, row.width*2, row.height])
tchm.ax_heatmap.set(xlabel = 'Week', ylabel=''u"\u2211"r' $\bar{x}$'" molar percentage")
clustermapstr='OUT_Taxonomy_phylum_graphs/TAXA_Class_clustermap_molar_percentage'+time.strftime('%H%M%S') + '.png'
#tchm.savefig(clustermapstr, dpi=1200)                               
plt.clf()
                                
####################order level taxa################

taxonomy_sum_dict_class={}
taxa_list_class=[]
taxlevel='order'
for k in df_mean.index:
      dictname_one_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'One')
      dictname_three_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'Three')
      dictname_five_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Five')
      dictname_ten_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Ten')
      dict_name_list=[dictname_one_tuple,dictname_three_tuple,dictname_five_tuple,dictname_ten_tuple]
      for key in dict_name_list:
          if key not in taxonomy_sum_dict_class:
              if taxa_master_dict_named[k]['Taxonomy'][taxlevel] not in taxa_list_class:
                  taxa_list_class.append(taxa_master_dict_named[k]['Taxonomy'][taxlevel])
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])
          if key in taxonomy_sum_dict_class:
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])+taxonomy_sum_dict_class[key]

week_list=['One','Three','Five','Ten']
taxa_index=pd.MultiIndex.from_product([taxa_list_class, week_list])
df3=pd.DataFrame(index=taxa_index)

df_taxa_sum=pd.DataFrame.from_dict(taxonomy_sum_dict_class, orient='index')
df3=pd.concat([df3, df_taxa_sum], axis=1, join='inner')
df3.reset_index(drop=False, inplace=True)
df3=df3.rename(columns={'level_0': 'Taxonomy', 'level_1':'Week', 0:'Molar_percentage'})
#create an error columns of zeros, as we cannot generate error bars for this data...
df3['err']=pd.Series(np.zeros(len(df3.index)), index=df3.index)
#Sort the week column so that they display in the correct order in the graphs...
custom_sort_dict={'One':1,'Three':3,'Five':5,'Ten':10}
df3['temp_rank']=df3['Week'].map(custom_sort_dict)
df3.sort_values(by=['temp_rank'],inplace=True)
df3.drop(labels=['temp_rank'],axis=1)
sns.set(style="white")


###########################TAXA clustermap
df_class_taxa_heatmap=df3.pivot(index='Taxonomy', columns='Week', values='Molar_percentage').reset_index().set_index('Taxonomy')
df_class_taxa_heatmap=df_class_taxa_heatmap.reindex(df_class_taxa_heatmap.sum(axis=1).rank(ascending=False).sort_values().index)
df_class_taxa_heatmap=df_class_taxa_heatmap[['One', 'Three', 'Five', 'Ten']]
###filter###
#filter_pct=1.5
#df_class_taxa_heatmap=df_class_taxa_heatmap[df_class_taxa_heatmap.sum(axis=1)*max(df_class_taxa_heatmap.sum(axis=1))*100>filter_pct]

####
sns.set(font_scale=1.6)
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_class_taxa_heatmap.values, df_class_taxa_heatmap.values.T))
tchm = sns.clustermap(df_class_taxa_heatmap, row_linkage=row_linkage, col_linkage=col_linkage, col_cluster=False, row_cluster=False, figsize=(7,12), cmap='Blues', linecolor='white', linewidths=2, annot=True, annot_kws={"size": 16}, fmt='.2g')
plt.setp(tchm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=15)
tchm.cax.set_position([1.275, .2, .03, .45])
plt.tick_params(which='y', size=14)
col = tchm.ax_col_dendrogram.get_position()
row = tchm.ax_row_dendrogram.get_position()
tchm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*0.25])
tchm.ax_row_dendrogram.set_position([row.x0*-.005, row.y0, row.width*2, row.height])
tchm.ax_heatmap.set(xlabel = 'Week', ylabel=''u"\u2211"r' $\bar{x}$'" molar percentage")
clustermapstr='OUT_Taxonomy_phylum_graphs/TAXA_'+str(taxlevel)+'_clustermap_molar_percentage'+time.strftime('%H%M%S') + '.png'

tchm.savefig(clustermapstr, dpi=1200)  
plt.clf()
                           

############family level taxa##############################
taxonomy_sum_dict_class={}
taxa_list_class=[]
taxlevel='family'
for k in df_mean.index:
      dictname_one_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'One')
      dictname_three_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'Three')
      dictname_five_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Five')
      dictname_ten_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Ten')
      dict_name_list=[dictname_one_tuple,dictname_three_tuple,dictname_five_tuple,dictname_ten_tuple]
      for key in dict_name_list:
          if key not in taxonomy_sum_dict_class:
              if taxa_master_dict_named[k]['Taxonomy'][taxlevel] not in taxa_list_class:
                  taxa_list_class.append(taxa_master_dict_named[k]['Taxonomy'][taxlevel])
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])
          if key in taxonomy_sum_dict_class:
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])+taxonomy_sum_dict_class[key]

week_list=['One','Three','Five','Ten']
taxa_index=pd.MultiIndex.from_product([taxa_list_class, week_list])
df3=pd.DataFrame(index=taxa_index)

df_taxa_sum=pd.DataFrame.from_dict(taxonomy_sum_dict_class, orient='index')
df3=pd.concat([df3, df_taxa_sum], axis=1, join='inner')
df3.reset_index(drop=False, inplace=True)
df3=df3.rename(columns={'level_0': 'Taxonomy', 'level_1':'Week', 0:'Molar_percentage'})
#create an error columns of zeros, as we cannot generate error bars for this data...
df3['err']=pd.Series(np.zeros(len(df3.index)), index=df3.index)
#Sort the week column so that they display in the correct order in the graphs...
custom_sort_dict={'One':1,'Three':3,'Five':5,'Ten':10}
df3['temp_rank']=df3['Week'].map(custom_sort_dict)
df3.sort_values(by=['temp_rank'],inplace=True)
df3.drop(labels=['temp_rank'],axis=1)
sns.set(style="white")

df_class_taxa_heatmap=df3.pivot(index='Taxonomy', columns='Week', values='Molar_percentage').reset_index().set_index('Taxonomy')
df_class_taxa_heatmap=df_class_taxa_heatmap.reindex(df_class_taxa_heatmap.sum(axis=1).rank(ascending=False).sort_values().index)
df_class_taxa_heatmap=df_class_taxa_heatmap[['One', 'Three', 'Five', 'Ten']]
###filter###
#filter_pct=1.5
#df_class_taxa_heatmap=df_class_taxa_heatmap[df_class_taxa_heatmap.sum(axis=1)*max(df_class_taxa_heatmap.sum(axis=1))*100>filter_pct]

############GENUS level taxa##############################
taxonomy_sum_dict_class={}
taxa_list_class=[]
taxlevel='genus'
for k in df_mean.index:
      dictname_one_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'One')
      dictname_three_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'Three')
      dictname_five_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Five')
      dictname_ten_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Ten')
      dict_name_list=[dictname_one_tuple,dictname_three_tuple,dictname_five_tuple,dictname_ten_tuple]
      for key in dict_name_list:
          if key not in taxonomy_sum_dict_class:
              if taxa_master_dict_named[k]['Taxonomy'][taxlevel] not in taxa_list_class:
                  taxa_list_class.append(taxa_master_dict_named[k]['Taxonomy'][taxlevel])
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])
          if key in taxonomy_sum_dict_class:
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])+taxonomy_sum_dict_class[key]

week_list=['One','Three','Five','Ten']
taxa_index=pd.MultiIndex.from_product([taxa_list_class, week_list])
df3=pd.DataFrame(index=taxa_index)

df_taxa_sum=pd.DataFrame.from_dict(taxonomy_sum_dict_class, orient='index')
df3=pd.concat([df3, df_taxa_sum], axis=1, join='inner')
df3.reset_index(drop=False, inplace=True)
df3=df3.rename(columns={'level_0': 'Taxonomy', 'level_1':'Week', 0:'Molar_percentage'})
#create an error columns of zeros, as we cannot generate error bars for this data...
df3['err']=pd.Series(np.zeros(len(df3.index)), index=df3.index)
#Sort the week column so that they display in the correct order in the graphs...
custom_sort_dict={'One':1,'Three':3,'Five':5,'Ten':10}
df3['temp_rank']=df3['Week'].map(custom_sort_dict)
df3.sort_values(by=['temp_rank'],inplace=True)
df3.drop(labels=['temp_rank'],axis=1)
sns.set(style="white")

df_genus_taxa_heatmap=df3.pivot(index='Taxonomy', columns='Week', values='Molar_percentage').reset_index().set_index('Taxonomy')
df_genus_taxa_heatmap=df_genus_taxa_heatmap.reindex(df_genus_taxa_heatmap.sum(axis=1).rank(ascending=False).sort_values().index)
df_genus_taxa_heatmap=df_genus_taxa_heatmap[['One', 'Three', 'Five', 'Ten']]
sns.heatmap(df_genus_taxa_heatmap)
############REAL CLASS LEVEL HEATMAP#######################
taxonomy_sum_dict_class={}
taxa_list_class=[]
taxlevel='class'
for k in df_mean.index:
      dictname_one_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'One')
      dictname_three_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]), 'Three')
      dictname_five_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Five')
      dictname_ten_tuple=(str(taxa_master_dict_named[k]['Taxonomy'][taxlevel]),'Ten')
      dict_name_list=[dictname_one_tuple,dictname_three_tuple,dictname_five_tuple,dictname_ten_tuple]
      for key in dict_name_list:
          if key not in taxonomy_sum_dict_class:
              if taxa_master_dict_named[k]['Taxonomy'][taxlevel] not in taxa_list_class:
                  taxa_list_class.append(taxa_master_dict_named[k]['Taxonomy'][taxlevel])
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])
          if key in taxonomy_sum_dict_class:
              taxonomy_sum_dict_class[key]=(df_mean.loc[k][key[1]])+taxonomy_sum_dict_class[key]

week_list=['One','Three','Five','Ten']
taxa_index=pd.MultiIndex.from_product([taxa_list_class, week_list])
df3=pd.DataFrame(index=taxa_index)

df_taxa_sum=pd.DataFrame.from_dict(taxonomy_sum_dict_class, orient='index')
df3=pd.concat([df3, df_taxa_sum], axis=1, join='inner')
df3.reset_index(drop=False, inplace=True)
df3=df3.rename(columns={'level_0': 'Taxonomy', 'level_1':'Week', 0:'Molar_percentage'})
#create an error columns of zeros, as we cannot generate error bars for this data...
df3['err']=pd.Series(np.zeros(len(df3.index)), index=df3.index)
#Sort the week column so that they display in the correct order in the graphs...
custom_sort_dict={'One':1,'Three':3,'Five':5,'Ten':10}
df3['temp_rank']=df3['Week'].map(custom_sort_dict)
df3.sort_values(by=['temp_rank'],inplace=True)
df3.drop(labels=['temp_rank'],axis=1)
sns.set(style="white")

df_classx_taxa_heatmap=df3.pivot(index='Taxonomy', columns='Week', values='Molar_percentage').reset_index().set_index('Taxonomy')
df_classx_taxa_heatmap=df_classx_taxa_heatmap.reindex(df_classx_taxa_heatmap.sum(axis=1).rank(ascending=False).sort_values().index)
df_classx_taxa_heatmap=df_classx_taxa_heatmap[['One', 'Three', 'Five', 'Ten']]
sns.heatmap(df_classx_taxa_heatmap)




####
sns.set(font_scale=1.6)
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_class_taxa_heatmap.values, df_class_taxa_heatmap.values.T))
tchm = sns.clustermap(df_class_taxa_heatmap, row_linkage=row_linkage, col_linkage=col_linkage, col_cluster=False, row_cluster=False, figsize=(7,12), cmap='Blues', linecolor='white', linewidths=2, annot=True, annot_kws={"size": 16}, fmt='.2g')
plt.setp(tchm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=15)
tchm.cax.set_position([1.275, .2, .03, .45])
plt.tick_params(which='y', size=14)
col = tchm.ax_col_dendrogram.get_position()
row = tchm.ax_row_dendrogram.get_position()
tchm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height*0.25])
tchm.ax_row_dendrogram.set_position([row.x0*-.005, row.y0, row.width*2, row.height])
tchm.ax_heatmap.set(xlabel = 'Week', ylabel=''u"\u2211"r' $\bar{x}$'" molar percentage")
clustermapstr='OUT_Taxonomy_phylum_graphs/TAXA_'+str(taxlevel)+'_clustermap_molar_percentage'+time.strftime('%H%M%S') + '.png'

tchm.savefig(clustermapstr, dpi=1200)  
plt.clf()

       #################################################################################################################                    


###############################COMMUNITY PROFILE DATA REQUIRED#########################################

#gamprod=df_class_taxa_heatmap.loc['Gammaproteobacteria']
#delprod=df_class_taxa_heatmap.loc['Deltaproteobacteria']
#bacteroidiaprod=df_class_taxa_heatmap.loc['Bacteroidia']
#flavoprod=df_class_taxa_heatmap.loc['Flavobacteriia']
#cytoprod=df_class_taxa_heatmap.loc['Cytophagia']
#sphingoprod=df_class_taxa_heatmap.loc['Sphingobacteriia']
#verrucoprod=df_class_taxa_heatmap.loc['Verrucomicrobiae']
#
#gamabd=np.array(df_filtered_normalised.loc['Gammaproteobacteria'])[[1,3,5,8]]
#delabd=array(df_filtered_normalised.loc['Deltaproteobacteria'])[[1,3,5,8]]
#bacteroidiaabd=array(df_filtered_normalised.loc['Bacteroidia'])[[1,3,5,8]]
#flavoabd=array(df_filtered_normalised.loc['Flavobacteriia'])[[1,3,5,8]]
#cytoabd=array(df_filtered_normalised.loc['Cytophagia'])[[1,3,5,8]]
#sphingoabd=array(df_filtered_normalised.loc['Sphingobacteriia'])[[1,3,5,8]]
#verrucoabd=array(df_filtered_normalised.loc['Verrucomicrobiae'])[[1,3,5,8]]
#
#tax=['Gammaproteobacteria','Deltaproteobacteria', 'Bacteroidia', 'Flavobacteriia', 'Cytophagia','Sphingobacteriia','Verrucomicrobiae']
#taxdict={}
#for x in tax:
#    taxdict[x]={}
#    taxdict[x]['abd']=np.array(df_filtered_normalised.loc[x])[[1,3,5,8]]
#    taxdict[x]['mol']=df_class_taxa_heatmap.loc[x]
#
#
#def abundance_vs_productivity_plot(abd, mol, name):
#     fig, ax = plt.subplots(figsize=(8,6),dpi=700) 
#     ax.plot([1,2,3,4], mol, marker='^', color='black',markersize=8)
#     ax2=ax.twinx()
#     ax2.plot([1,2,3,4], abd, marker='s', color='blue',markersize=8)
#     ax2.set_ylabel('Relative abundance (%)', color='blue',  fontsize=20)
#     ax.set_ylabel(''u"\u2211"r' $\bar{x}$'" molar percentage", color='black',  fontsize=20)
#     plt.tick_params(axis='both', labelsize=20)
#     ax2.tick_params(axis='both', direction='out', length=3, width=2,  labelsize=20)
#     ax.tick_params(axis='both', direction='out', length=3, width=2, labelsize=20)
#     ax.spines['top'].set_visible(False)
#     ax2.spines['top'].set_visible(False)
#     ax.set_xticks([1,2,3,4])
#     ax2.set_xticks([1,2,3,4])
#     ax.set_xticklabels(['One','Three', 'Five', 'Ten'])
#     ax2.set_xticklabels(['One','Three', 'Five', 'Ten'])
#     ax.set_xlabel('Week', fontsize=20)
#     plt.title(name, size=20)
#     return fig
#    
#    
#    
#fig, ax = plt.subplots(figsize=(8,6),dpi=700) 
#ax.plot([1,2,3,4], gamprod, marker='^', color='blue',markersize=8)
#ax2=ax.twinx()
#ax2.plot([1,2,3,4], gamabd, marker='s', color='black',markersize=8)
#ax2.set_ylabel('Relative abundance (%)', color='blue',  fontsize=20)
#ax.set_ylabel(''u"\u2211"r' $\bar{x}$'" molar percentage", color='black',  fontsize=20)
#plt.tick_params(axis='both', labelsize=20)
#ax2.tick_params(axis='both', direction='out', length=3, width=2,  labelsize=20)
#ax.tick_params(axis='both', direction='out', length=3, width=2, labelsize=20)
#ax.spines['top'].set_visible(False)
#ax2.spines['top'].set_visible(False)
#ax.set_xticks([1,2,3,4])
#ax2.set_xticks([1,2,3,4])
#ax.set_xticklabels(['One','Three', 'Five', 'Ten'])
#ax2.set_xticklabels(['One','Three', 'Five', 'Ten'])
#ax.set_xlabel('Week', fontsize=20)

######################################################################################################

############################
#class_percent_filter=1
#print('Responsible taxa based on sequence homology at the class level, less than ' + str(class_percent_filter) + ' filtered')
#
#df3_filter=df3.groupby(['Taxonomy'])['Molar_percentage'].sum().reset_index()
#df3_filter['contribution']=df3_filter['Molar_percentage']/(df3_filter['Molar_percentage'].sum())*100
#df3_filter.set_index('Taxonomy')
#df3_filter[df3_filter['contribution']>class_percent_filter]['Taxonomy']
#df3_filter=df3[df3['Taxonomy'].isin(df3_filter[df3_filter['contribution']>class_percent_filter]['Taxonomy'])]
#sns.set(style="white", font_scale=1.6)
#t3_f = sns.factorplot(data=df3_filter, x="Week", y="Molar_percentage", col="Taxonomy", kind="bar", col_wrap=3, palette=sns.husl_palette(6, h=.5))
#t3_f.set_yticklabels(fontsize=20)
#t3_f.set_xticklabels(fontsize=20)
#t3_f.set_ylabels(u"\u2211" r'$\bar{x}$'" Molar percentage", fontsize=20)
#t.set_xlabels(fontsize=20)
#t3_f.savefig("OUT_taxonomy_phylum_graphs/Top_class_less_than_"+str(class_percent_filter)+"_dpi_500.png")
#plt.show(t3_f)
#plt.clf()


#print('Top three responsible class only') 
##tax_retain_class=['Bacteroidia','Cytophagia','Flavobacteriia', 'Deltaproteobacteria', 'Gammaproteobacteria']
#tax_retain_class=['Vibrionales','Bacteroidales','Chitinophagales', 'Cellvibrionales', 'Flavobacteriales']
#df3_top=df3_filter[df3_filter['Taxonomy'].isin(tax_retain_class)]
#t_f_x=sns.factorplot(data=df3_top, x="Week", y="Molar_percentage", col="Taxonomy", kind="bar", col_wrap=3, palette=sns.husl_palette(6, h=.5))
#t_f_x.set_ylabels(u"\u2211" r'$\bar{x}$'" Molar percentage")
#plt.show(t_f_x)



#############Create a list of targets by bringing all this info together##############
#Create a rank based on overall abundance of enzyme groups, e.g. GH6 is #1, CE1 is #2 etc...
#new dataframe
#make new folder to store graphs in...
df_rank=pd.DataFrame(index=df_sum_mean.index)    
df_rank['total_sum']=df_sum_mean['One']+df_sum_mean['Three']+df_sum_mean['Five']+df_sum_mean['Ten']
df_rank['rank']=df_rank['total_sum'].rank(ascending=False)
df_rank[df_rank['rank']<extract_top_n_class].sort_values(by=['rank'])
top_target_list=[]
top_target_list=df_rank[df_rank['rank']<=extract_top_n_class].sort_values(by=['rank']).index.tolist()

for x in manual_target_list:
    for y in df_rank.index:
        if bool(re.search(x, y))==True:
            if y not in top_target_list:
                top_target_list.append(y)

#add a comma (,) to the end of the enzyme class so they are picked up in the names e.g. (GH6, CBM...)
top_target_list_comma=list(x +',' for x in top_target_list)
#add a space to the class so singular domain ORFs are picked up e.g. CE1 (a_...)
top_target_list_space=list(x +' ' for x in top_target_list)
top_target_list_none=list(x +'' for x in top_target_list)
top_target_list_final=top_target_list_comma+top_target_list_space+top_target_list_none

########################################################################################
target_ORF_dict={}
for key in copy_crossref_nested_dict_named:
    for ki in top_target_list_final:
        if re.search(ki, key):
            target_ORF_dict[key]={}
############
#Append CBM only ORFs if required########
if manual_entry_of_CBM_only_ORFs.lower()=='yes':
    for x in CBM_only_ORF_list:
        if x not in target_ORF_dict:
            target_ORF_dict[x]={}
            target_ORF_dict[x]['Class_rank']='None'

#####CLASS RANK########
for x in df_rank.index:
    for y in target_ORF_dict:
         if bool(re.search(x, y))==True:
             if re.search(x, y):
                 target_ORF_dict[y]['Class_rank']=df_rank['rank'][x]
#######ORF_RANK############
df_ORF_rank=pd.DataFrame(index=df_mean.index)
df_ORF_rank['sum']=df_mean['One']+df_mean['Three']+df_mean['Five']+df_mean['Ten']
df_ORF_rank['ORFrank']=df_ORF_rank['sum'].rank(ascending=False)

for x in df_ORF_rank.index:
    for y in target_ORF_dict:
        if x == y:
            target_ORF_dict[y]['ORF_rank']=df_ORF_rank['ORFrank'][x]
###############PDF_COUNT_DICT##############
pdf_count_dict={}
for x in top_target_list:
    count=0
    for y in target_ORF_dict:
        if bool(re.search(x+',', y))==True or bool(re.search(x+' ', y))==True:
            count+=1
    if df_rank.loc[x]['rank'] not in pdf_count_dict:
         pdf_count_dict[df_rank.loc[x]['rank']]={}
    pdf_count_dict[df_rank.loc[x]['rank']]['count']=count
    pdf_count_dict[df_rank.loc[x]['rank']]['class']=str(x)        
############################################
###############################Draw GRAPHS and VECTOR IMAGES##################

#key='PL8 (a_c2293733_g1_i1_1)'
def plot_barplot(key):
   sample_to_week_dict={'W1_C201':'One','W1_3_C429':'One','W1_4_C429':'One','W3_C201':'Three','W3_2_C429':'Three','W3_1_C429':'Three','W5_C201':'Five','W5_25_C429':'Five','W5_100_C429':'Five','W10_C201':'Ten','W10_100_i_C429':'Ten','W10_100_ii_C429':'Ten'}
   key=str(key)
   target_index=df.index.get_loc(key)
   temp_df=df.iloc[target_index][['W1_C201','W1_3_C429', 'W1_4_C429','W3_C201','W3_2_C429','W3_1_C429', 'W5_C201','W5_25_C429','W5_100_C429','W10_C201','W10_100_i_C429','W10_100_ii_C429']]
   temp_df=temp_df.reset_index()
   temp_df['Week']=temp_df['index'].map(sample_to_week_dict)
   keyname=r'$\bar{x}$''Molar percentage\n'+str(key)
   temp_df=temp_df.rename(columns={key:keyname})
   tempgraph=sns.barplot(data=temp_df, x='Week', y=keyname, palette='hls', ci=68)
   tempgraph.set_ylabel(keyname)
   sns.despine(bottom=False)
   savnam=((''.join(reversed((''.join(reversed(key))).split(' ', 1)[0])).strip('()')))
   savnam='OUT_target_images/'+savnam+'_molpct_graph_sterr.png'
   plt.savefig(savnam)
   #return plt.show(tempgraph)
   plt.clf()
   plt.cla()
   plt.close()
   
#plot_barplot(key)

def draw_skeleton(length, ORFlength, ORF_name):
   turtle.hideturtle()
   turtle.speed(0)
   turtle.tracer(0,0)
   turtle.color("blue")
   turtle.back(0.5*length)
   turtle.forward(length)
   turtle.penup()
   turtle.setpos(0,70)
   turtle.pendown()
   turtle.back(0.5*length)
   turtle.forward(length)
   turtle.penup()
   turtle.color("red")
   turtle.setpos(0,-70)
   turtle.pendown()
   turtle.back(0.5*length)
   turtle.forward(length)
   turtle.penup()
#   turtle.color("green")
#   turtle.setpos(0,-140)
#   turtle.pendown()
#   turtle.back(0.5*length)
#   turtle.forward(length)
#   turtle.penup()
   turtle.color("blue")
   bars=18
   x=ORFlength/bars
   xlist=list(range(0, bars))
   xtickdist=length/(10*bars)
   turtle.setpos(-0.5*length,70)
   for x in xlist:
     turtle.pendown()
     turtle.rt(270)
     turtle.forward(15)
     number=round(ORFlength/bars*x)
     txt=str(number)
     turtle.write(txt, align='center')
     turtle.back(15)
     turtle.rt(90)
     turtle.forward(xtickdist)
     for y in list(range(9)):
           turtle.rt(270)
           turtle.forward(8)
           turtle.back(8)
           turtle.rt(90)
           turtle.forward(xtickdist)
   turtle.rt(90)
   turtle.color("blue")
   turtle.begin_fill()
   turtle.forward(10)
   turtle.rt(225)
   turtle.forward(20)
   turtle.rt(270)
   turtle.forward(20)
   turtle.write(str(ORFlength), align='center')
   turtle.rt(225)
   turtle.forward(20)
   turtle.end_fill() 
   turtle.penup()
   turtle.setpos(0, 100)
   name=str(ORF_name)
   turtle.write(name, align='center',font=("Arial", 8, "bold"))
   turtle.setpos(0, 20)
   turtle.write('dbCAN domain architecture', align='center',font=("Arial", 8, "bold"))
   turtle.setpos(0, -50)
   turtle.color("red")
   turtle.write('Pfam domain architecture', align='center',font=("Arial", 8, "bold"))
#   turtle.setpos(0, -120)
#   turtle.color("green")
#   turtle.write('BlastP domain architecture', align='center',font=("Arial", 8, "bold"))
   turtle.update()

def plot_box(start, end, length, ORFlength, domain_name, e_val, colour, level):
    if level == 0:
        y=0
    elif level == 1:
        y=-70
    elif level ==2:
        y=-140
    elif level==3:
        y=-150
    turtle.penup()
    width=5
    turtle.hideturtle()
    turtle.speed(0)
    turtle.tracer(0,0)
    if start >= 0.5*ORFlength and end >= 0.5*ORFlength:
         start_coord=length/ORFlength*(start-0.5*ORFlength)
         end_coord=length/ORFlength*(end-0.5*ORFlength)
         turtle.setpos(start_coord,y)
         turtle.setheading(180)
         turtle.pendown()
         turtle.color(colour)
         turtle.begin_fill()
         turtle.rt(90)
         turtle.fd(width)
         startstr=str(start)
         turtle.write(startstr, align='left')
         turtle.rt(90)
         turtle.fd(end_coord-start_coord)
         endstr=str(end)
         turtle.write(endstr, align='right')
         turtle.rt(90)
         turtle.fd(width*2)
         turtle.rt(90)
         turtle.fd(end_coord-start_coord)
         turtle.rt(90)
         turtle.fd(width)
         turtle.end_fill()
         turtle.penup()
         turtle.setpos(start_coord+(0.5*(end_coord-start_coord)), y-20)
         domain_name=(str(domain_name))
         turtle.write(domain_name, align='center')
         turtle.setpos(start_coord+(0.5*(end_coord-start_coord)), y-35)
         evalu=str(e_val)
         turtle.write(evalu, align='center')
    elif start <= 0.5*ORFlength and end <= 0.5*ORFlength:
         start_coord=-0.5*length+((length/ORFlength*start))#-ORFlength
         end_coord=-0.5*length+((length/ORFlength*end))
         turtle.setpos(start_coord,y)
         turtle.setheading(180)
         turtle.pendown()
         turtle.color(colour)
         turtle.begin_fill()
         turtle.rt(90)
         turtle.fd(width)
         startstr=str(start)
         turtle.write(startstr, align='left')
         turtle.rt(90)
         turtle.fd((start_coord-end_coord)*-1)
         endstr=str(end)
         turtle.write(endstr, align='right')
         turtle.rt(90)
         turtle.fd(width*2)
         turtle.rt(90)
         turtle.fd((start_coord-end_coord)*-1)
         turtle.rt(90)
         turtle.fd(width)
         turtle.end_fill()
         turtle.penup()
         turtle.setpos(start_coord-(-0.5*(end_coord+(start_coord*-1))), y-20)
         domain_name=(str(domain_name))
         turtle.write(domain_name, align='center') 
         turtle.setpos(start_coord-(-0.5*(end_coord+(start_coord*-1))), y-35)
         evalu=str(e_val)
         turtle.write(evalu, align='center')
    elif start <= 0.5*ORFlength and end >= 0.5*ORFlength:
         start_coord=-0.5*length+((length/ORFlength*start))#-ORFlength
         end_coord=length/ORFlength*(end-0.5*ORFlength)
         turtle.setpos(start_coord,y)
         turtle.setheading(180)
         turtle.pendown()
         turtle.color(colour)
         turtle.begin_fill()
         turtle.rt(90)
         turtle.fd(width)
         startstr=str(start)
         turtle.write(startstr, align='left')
         turtle.rt(90)
         turtle.fd(end_coord+start_coord*-1)
         endstr=str(end)
         turtle.write(endstr, align='right')
         turtle.rt(90)
         turtle.fd(width*2)
         turtle.rt(90)
         turtle.fd(end_coord+start_coord*-1)
         turtle.rt(90)
         turtle.fd(width)
         turtle.end_fill()
         turtle.penup()
         turtle.setpos(start_coord-(-0.5*(end_coord+(start_coord*-1))), y-20)
         domain_name=(str(domain_name))
         turtle.write(domain_name, align='center')
         turtle.setpos(start_coord-(-0.5*(end_coord+(start_coord*-1))), y-35)
         evalu=str(e_val)
         turtle.write(evalu, align='center')
    turtle.update()
    

#draw_skeleton(900, 450, 'titleoffile')
#plot_box(10, 130, 900, 450, 'GH6', 1e-21, "blue", 0)
#plot_box(10, 130, 900, 450, 'GH6', 1e-21, "red", 1)
#plot_box(10, 130, 900, 450, 'GH6', 1e-21, "green", 2)
#plot_box(371, 417, 900, 450, 'CBM', 1e-160, "blue", 0)
#plot_box(371, 417, 900, 450, 'CBM', 1e-160, "red", 1)
#plot_box(371, 417, 900, 450, 'CBM', 1e-160, "green", 2)
#plot_box(170, 269, 900, 450, 'AA2', 9e-999, "blue", 0)
#plot_box(170, 269, 900, 450, 'AA2', 9e-999, "red", 1)
#plot_box(170, 269, 900, 450, 'AA2', 9e-999, "green", 2)

#def save_the_turtles(ORF_name):
#    name='OUT_target_images/'+str(ORF_name)+'_file'
#    ts=turtle.getscreen().getcanvas()
#    canvasvg.saveall(name+'.svg', ts)
#    with open(name+'.svg') as svg_input, open(name+'.png', 'wb') as png_out:
#        cairosvg.svg2png(bytestring=svg_input.read(), write_to=png_out)

#draw_skeleton(900, 450, 'titleoffile')
#plot_box(10, 130, 900, 450, 'GH6', 1e-2110, "blue")
#plot_box(371, 417, 900, 450, 'CBM', 1e-160, "blue")
#plot_box(170, 269, 900, 450, 'AA2', 9e-999, "blue")
#save_the_turtles('titleoffile')  
#turtle.reset()

#####################Extract ORF length for diagrammatic representation################
seqs={}
annotation_dict_length={}#copy.deepcopy(annotation_dict)
for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
    for key in target_ORF_dict:
        ki=((''.join(reversed((''.join(reversed(key))).split(' ', 1)[0])).strip('()')))
        if seq.id == ki:
            annotation_dict_length[ki]=annotation_dict[ki]
            annotation_dict_length[ki]['length']=len(seq.seq)
            target_ORF_dict[key]['sequence']=seq.seq
################TURTLES#################################################################



#for key in annotation_dict_length:
#    lengthbp=annotation_dict_length[key]['length']
#    titlebp=str(key)
#    draw_skeleton(vector_length, lengthbp, titlebp)
#    for key2 in annotation_dict_length[key]:
#        if not key2 == 'length':
#            startbp=round(annotation_dict_length[key][key2]['start'])
#            endbp=round(annotation_dict_length[key][key2]['end'])
#            domainx=str(key2.split('_',1)[0])
#            evalu=annotation_dict_length[key][key2]['e_val']
#            plot_box(startbp, endbp, vector_length, lengthbp, domainx, evalu, "blue", 0)
#    if key in pfam_annotation_dict:
#        for key2 in pfam_annotation_dict[key]:
#            startbp=round(pfam_annotation_dict[key][key2]['start'])
#            endbp=round(pfam_annotation_dict[key][key2]['end'])
#            domainx='_'.join(key2.split('_', (key2.count('_')))[:(key2.count('_'))])
#            evalu=pfam_annotation_dict[key][key2]['e_val']
#            plot_box(startbp, endbp, vector_length, lengthbp, domainx, evalu, "red", 1)          
#    filsav=str(key)
#    save_the_turtles(filsav)
#    turtle.reset()
#      
###########Compile all data into target dictionary##################    
#f=xlrd.open_workbook(blastp_tophit)
#sh=f.sheet_by_index(0)
#for rownum in range(0,sh.nrows):
#    for k in target_ORF_dict:
#        ki=((''.join(reversed((''.join(reversed(k))).split(' ', 1)[0])).strip('()')))
#        graphname='OUT_target_images/'+ki+'_molpct_graph_sterr.png'
#        vectorname='OUT_target_images/'+ki+'_file.png'
#        if sh.cell_value(rownum,0) == ki:
#            target_ORF_dict[k]['BlastP_tophit']=sh.cell_value(rownum,3)
#            target_ORF_dict[k]['BlastP_e_value']=sh.cell_value(rownum,4)
#            target_ORF_dict[k]['BlastP_ID']=sh.cell_value(rownum,5)
#            target_ORF_dict[k]['BlastP_coverage']=sh.cell_value(rownum,6)
#            target_ORF_dict[k]['lineage']=taxa_master_dict[ki]['Taxonomy']['lineage']
#            target_ORF_dict[k]['totalrank']=df_master_mean.loc[ki]['rank']
#            plot_barplot(k)
#            target_ORF_dict[k]['graph']=graphname
#            target_ORF_dict[k]['vector']=vectorname   
#non_blastP_annot_list=[]
#for rownum in range(0,sh.nrows):
#       non_blastP_annot_list.append(sh.cell_value(rownum,0))
#for k in target_ORF_dict:
#        ki=((''.join(reversed((''.join(reversed(k))).split(' ', 1)[0])).strip('()')))
#        if ki not in non_blastP_annot_list:
#            target_ORF_dict[k]['BlastP_tophit']='None'
#            target_ORF_dict[k]['BlastP_e_value']='None'
#            target_ORF_dict[k]['BlastP_ID']='None'
#            target_ORF_dict[k]['BlastP_coverage']='None'
#            target_ORF_dict[k]['lineage']='None'
#            target_ORF_dict[k]['totalrank']=df_master_mean.loc[ki]['rank']
#            plot_barplot(k)
#            target_ORF_dict[k]['graph']=graphname
#            target_ORF_dict[k]['vector']=vectorname  
#df_master_ts=pd.DataFrame()
#df_ts=copy.deepcopy(df_sum_activity)
#ts_col_pal=sns.color_palette("husl", len(df_sum_activity.index))
#for x in range(0,len(df_sum_activity.index)):    
#   df_ts=copy.deepcopy(df_sum_activity)
#  #del df_ts['taxonomy']
#   df_temp=df_ts.iloc[x].reset_index()
#   classx=str(df_temp[[1]].columns.tolist()[0])
#   df_temp=df_temp.rename(columns={'index': 'sample', classx:'molar_percentage'})
#   class_series=pd.Series(classx, index=range(0,len(df_temp[[1]])))
#   df_temp['class']=class_series
#   sample_to_week_dict={'W1_C201':'1','W1_3_C429':'1','W1_4_C429':'1','W3_C201':'3','W3_2_C429':'3','W3_1_C429':'3','W5_C201':'5','W5_25_C429':'5','W5_100_C429':'5','W10_C201':'10','W10_100_i_C429':'10','W10_100_ii_C429':'10'}
#   df_temp['week']=df_temp['sample'].map(sample_to_week_dict)       
#   df_master_ts=pd.concat([df_master_ts, df_temp])   
#   sns.set_style()
#   plt.plot(df_temp['week'], df_temp['molar_percentage'], label=classx, color=ts_col_pal[x])    
#df_master_ts=df_master_ts.reset_index() 
#plt.show()
#
#ax=sns.tsplot(time='week', value='molar_percentage', unit='sample', condition='class', data=df_master_ts)
#
#if len(df_master_ts.index)/12 == len(df_sum_activity.index):
#    print('Successfully merged all data into a timeseries dataframe with ' + str(len(df_sum_activity.index)) + ' enzyme classes')
#if not len(df_master_ts.index)/12 == len(df_sum_activity.index):
#    print('***WARNING***\nFailed to merge all data into a timeseries dataframe...')

###############POP TARGETS####################
popkeys=['GH6', 'GH5', 'GH3', 'GH2', 'GH10', 'GH11', 'GH16', 'GH13', 'GH9']
popkeycomma=[x+',' for x in popkeys]
popkeyspace=[x+' ' for x in popkeys]
popkeyfinal=popkeycomma+popkeyspace
keeplist=[]
poplist=[]
for x in target_ORF_dict:
    for y in popkeyfinal:
        if bool(re.search(y, x))==True:
            if x not in keeplist:
                keeplist.append(x)
for k in target_ORF_dict:
    if k not in keeplist:
        poplist.append(k)
for k in poplist:
    target_ORF_dict.pop(k)
    
#################################################  
Proteobacteria_classes=list(df_sum_activity_taxa_dict_molpct_w1[df_sum_activity_taxa_dict_molpct_w1['Proteobacteria']>0].index)
Proteobacteria_classes=Proteobacteria_classes+list(df_sum_activity_taxa_dict_molpct_w3[df_sum_activity_taxa_dict_molpct_w3['Proteobacteria']>0].index)
Proteobacteria_classes=Proteobacteria_classes+list(df_sum_activity_taxa_dict_molpct_w5[df_sum_activity_taxa_dict_molpct_w5['Proteobacteria']>0].index)
Proteobacteria_classes=Proteobacteria_classes+list(df_sum_activity_taxa_dict_molpct_w10[df_sum_activity_taxa_dict_molpct_w10['Proteobacteria']>0].index)

Bacteroidetes_classes=list(df_sum_activity_taxa_dict_molpct_w1[df_sum_activity_taxa_dict_molpct_w1['Bacteroidetes']>0].index)
Bacteroidetes_classes=Bacteroidetes_classes+list(df_sum_activity_taxa_dict_molpct_w3[df_sum_activity_taxa_dict_molpct_w3['Bacteroidetes']>0].index)
Bacteroidetes_classes=Bacteroidetes_classes+list(df_sum_activity_taxa_dict_molpct_w5[df_sum_activity_taxa_dict_molpct_w5['Bacteroidetes']>0].index)
Bacteroidetes_classes=Bacteroidetes_classes+list(df_sum_activity_taxa_dict_molpct_w10[df_sum_activity_taxa_dict_molpct_w10['Bacteroidetes']>0].index)

Verrucomicrobia_classes=list(df_sum_activity_taxa_dict_molpct_w1[df_sum_activity_taxa_dict_molpct_w1['Verrucomicrobia']>0].index)
Verrucomicrobia_classes=Verrucomicrobia_classes+list(df_sum_activity_taxa_dict_molpct_w3[df_sum_activity_taxa_dict_molpct_w3['Verrucomicrobia']>0].index)
Verrucomicrobia_classes=Verrucomicrobia_classes+list(df_sum_activity_taxa_dict_molpct_w5[df_sum_activity_taxa_dict_molpct_w5['Verrucomicrobia']>0].index)
Verrucomicrobia_classes=Verrucomicrobia_classes+list(df_sum_activity_taxa_dict_molpct_w10[df_sum_activity_taxa_dict_molpct_w10['Verrucomicrobia']>0].index)

def construct_string_with_endlines(input_list, n_per_line):
    fullstr=''
    if n_per_line==1:
         for x in input_list:
             fullstr=fullstr+x+'\n'
    if n_per_line>1:
        count=0
        for x in input_list:
             if count%n_per_line==0:
                  fullstr=fullstr+x+', '
                  count+=1
             else:
                  fullstr=fullstr+x+'\n'
                  count+=1
    return fullstr

#fig=plt.figure(1, figsize=(18, 18))
#v=venn3([set(Proteobacteria_classes), set(Bacteroidetes_classes), set(Verrucomicrobia_classes)])
#Prot_only=list(set(Proteobacteria_classes)-set(Bacteroidetes_classes))
#Bac_only=list(set(Bacteroidetes_classes)-set(Proteobacteria_classes))
#Prot_Bac_shared=list(set(Proteobacteria_classes).intersection(Bacteroidetes_classes)) 
#      
#protstr=construct_string_with_endlines(Prot_only, 1)
#bacstr=construct_string_with_endlines(Bac_only, 1)
#protbacsharestr=construct_string_with_endlines(Prot_Bac_shared, 2)
#
#v.get_label_by_id('A').set_text('Proteobacteria')
#v.get_label_by_id('B').set_text('Bacteroidetes')
#v.get_label_by_id('C').set_text('Verrucomicrobia')
#alpha=0.6
#
##red
#v.get_patch_by_id('100').set_alpha(alpha)
#v.get_patch_by_id('100').set_color((1, 0, 0))
##green
#v.get_patch_by_id('010').set_color((0, 1, 0))
#v.get_patch_by_id('010').set_alpha(alpha)
##blue
#v.get_patch_by_id('001').set_color((0, 0, 1))
#v.get_patch_by_id('001').set_alpha(alpha)
##yellow
#v.get_patch_by_id('110').set_color((1, 1, 0))
#v.get_patch_by_id('110').set_alpha(alpha)
##white
#v.get_patch_by_id('111').set_color((1, 1, 1))
#v.get_patch_by_id('111').set_alpha(alpha)
##orange
#v.get_patch_by_id('101').set_color((1, .75, .25))
#v.get_patch_by_id('101').set_alpha(alpha)
#
#v.get_patch_by_id('011').set_color((.75, .5, 1))
#v.get_patch_by_id('011').set_alpha(alpha)
#
#for text in v.set_labels:
#    text.set_fontsize(25)
#for text in v.subset_labels:
#    text.set_fontsize(25)
#
#plt.text(-.575, -0.2125, protstr, fontsize=17.5)
#plt.text(0.475, -0.2375, bacstr, fontsize=17.5)
#plt.text(-0.075,0.1, protbacsharestr, fontsize=17.5)
#plt.savefig('OUT_venn/top_phyla.png')
#plt.show()

##################ADD STANDARD DEVIATION TO COLUMNS IN DF MEAN FOR GRAPHS ETC#################
df_sum_mean['One_err']=df_sum_activity[['W1_C201','W1_3_C429', 'W1_4_C429']].std(axis=1)
df_sum_mean['Three_err']=df_sum_activity[['W3_C201','W3_2_C429','W3_1_C429']].std(axis=1)
df_sum_mean['Five_err']=df_sum_activity[['W5_C201','W5_25_C429','W5_100_C429']].std(axis=1)
df_sum_mean['Ten_err']=df_sum_activity[['W10_C201','W10_100_i_C429','W10_100_ii_C429']].std(axis=1)


df_sum_mean['sum']=df_sum_mean['One']+df_sum_mean['Three']+df_sum_mean['Five']+df_sum_mean['Ten']
df_sum_mean['class_rank']=df_sum_mean['sum'].rank(ascending=False)
top_x_for_barplots=15
print('Plotting mean molar percentage for top ' + str(top_x_for_barplots) + ' enzyme classes')
class_for_barplots_list=[]
class_for_barplots_list=df_sum_mean[df_sum_mean['class_rank']<=top_x_for_barplots].sort_values(by=['class_rank']).index.tolist()
f, axarr = plt.subplots(3, 5, figsize=(15, 9))
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)
colour_barplot=sns.xkcd_palette(['clear blue'])
#axarr[row, col]
sns.set_style="white"
for enz_class in list(range(top_x_for_barplots)):
     xvals=list(df_sum_mean.loc[class_for_barplots_list[enz_class]][:4])
     width=0.65
     ind=[0, 1, 2, 3]
     yerrx=list(df_sum_mean.loc[class_for_barplots_list[enz_class]][4:8])
     if enz_class<=4:
          axarr[0, enz_class].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
          axarr[0, enz_class].set_title(class_for_barplots_list[enz_class], size=15)
          axarr[0, enz_class].set_ylim([0, 0.35])
          axarr[0, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=17)
     elif enz_class<=9:
          axarr[1, enz_class-5].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
          axarr[1, enz_class-5].set_title(class_for_barplots_list[enz_class], size=15)
          axarr[1, enz_class-5].set_ylim([0, 0.35])
          axarr[1, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=17)
     elif enz_class<=14:
          axarr[2, enz_class-10].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
          axarr[2, enz_class-10].set_title(class_for_barplots_list[enz_class], size=15)
          axarr[2, enz_class-10].set_ylim([0, 0.35])
          axarr[2, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=17)
          axarr[2, enz_class-10].set_xlabel('Week', fontsize=17)
          axarr[2, enz_class-10].set_xticks([0, 1, 2, 3])          
          axarr[2, enz_class-10].set_xticklabels(['One', 'Three', 'Five', 'Ten'])
plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 3]], visible=False)
plt.setp([a.get_yticklabels() for a in axarr[:, 4]], visible=False)
sns.despine(top=True, right=True)
f.savefig('OUT_grids/top_'+str(top_x_for_barplots)+'_enzyme_classes.png')
plt.show()
plt.clf()
plt.close(f)  
######################################################################################################################
def make_barplot_all_index(nam, input_dataframe):
     num_plots=len(input_dataframe)
     if num_plots==1:
          f, axarr = plt.subplots(1, sharey=True, figsize=(3, 3))
     elif num_plots<=5:
          f, axarr = plt.subplots(1, num_plots, sharey=True, figsize=(num_plots*3, 3))
     elif num_plots<=6:
          f, axarr = plt.subplots(2, 3, figsize=(9, 6))
     elif num_plots<=10:
          f, axarr = plt.subplots(2, 5, figsize=(9, 6))
     elif num_plots<=15:
          f, axarr = plt.subplots(3, 5, figsize=(15, 9))     
     matplotlib.rc('xtick', labelsize=16)
     matplotlib.rc('ytick', labelsize=18)
     colour_barplot=sns.xkcd_palette(['clear blue'])
     ylimit=input_dataframe[input_dataframe.columns[:4]].values.max()*1.25
     #axarr[row, col]
     sns.set(style="white")
     print(num_plots)
     for enz_class in list(range(num_plots)):
          xvals=list(input_dataframe.loc[input_dataframe.index[enz_class]][:4])
          width=0.65
          ind=[0, 1, 2, 3]
          yerrx=list(input_dataframe.loc[input_dataframe.index[enz_class]][4:8])
          name=str(input_dataframe.iloc[enz_class].name)
          if num_plots==1:
              axarr.bar(ind, xvals, width, color=colour_barplot, yerr=yerrx) 
              axarr.set_title(name, size=15)
              axarr.set_ylim([0, ylimit])
              axarr.set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
              axarr.set_xticks([0, 1, 2, 3])          
              axarr.set_xticklabels(['One', 'Three', 'Five', 'Ten'], fontsize=16)
              axarr.set_xlabel('Week', fontsize=18)
          elif num_plots<=5:
                    axarr[enz_class].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                    axarr[enz_class].set_title(name, size=15)
                    axarr[enz_class].set_ylim([0, ylimit])
                    axarr[0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
                    axarr[enz_class].set_xticks([0, 1, 2, 3])          
                    axarr[enz_class].set_xticklabels(['One', 'Three', 'Five', 'Ten'], fontsize=16)
                    axarr[enz_class].set_xlabel('Week', fontsize=18)
          elif num_plots<=6:
              if enz_class<=2:
                    axarr[0, enz_class].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                    axarr[0, enz_class].set_title(name, size=16)
                    axarr[0, enz_class].set_ylim([0, ylimit])
                    axarr[0, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16) 
              elif enz_class<=5:
                    axarr[1, enz_class-3].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                    axarr[1, enz_class-3].set_title(name, size=16)
                    axarr[1, enz_class-3].set_ylim([0, ylimit])
                    axarr[1, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
                    axarr[1, enz_class-3].set_xticks([0, 1, 2, 3])          
                    axarr[1, enz_class-3].set_xticklabels(['One', 'Three', 'Five', 'Ten'],fontsize=16)
                    axarr[1, enz_class-3].set_xlabel('Week', fontsize=18)
          elif num_plots<=10:
              if enz_class<=4:
                    axarr[0, enz_class].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                    axarr[0, enz_class].set_title(name, size=16)
                    axarr[0, enz_class].set_ylim([0, ylimit])
                    axarr[0, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
              elif enz_class<=9:
                   axarr[1, enz_class-5].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                   axarr[1, enz_class-5].set_title(name, size=16)
                   axarr[1, enz_class-5].set_ylim([0, ylimit])
                   axarr[1, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
                   axarr[1, enz_class-5].set_xlabel('Week', fontsize=18)
                   axarr[1, enz_class-5].set_xticks([0, 1, 2, 3])          
                   axarr[1, enz_class-5].set_xticklabels(['One', 'Three', 'Five', 'Ten'],fontsize=16)
          elif num_plots<=15:
              if enz_class<=4:
                    axarr[0, enz_class].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                    axarr[0, enz_class].set_title(name, size=16)
                    axarr[0, enz_class].set_ylim([0, ylimit])
                    axarr[0, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
              elif enz_class<=9:
                   axarr[1, enz_class-5].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                   axarr[1, enz_class-5].set_title(name, size=16)
                   axarr[1, enz_class-5].set_ylim([0, ylimit])
                   axarr[1, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
              elif enz_class<=14:
                   axarr[2, enz_class-10].bar(ind, xvals, width, color=colour_barplot, yerr=yerrx)
                   axarr[2, enz_class-10].set_title(name, size=16)
                   axarr[2, enz_class-10].set_ylim([0, ylimit], fontsize=16)
                   axarr[2, 0].set_ylabel(r' $\bar{x}$'" molar percentage", fontsize=16)
                   axarr[2, enz_class-10].set_xlabel('Week', fontsize=18)
                   axarr[2, enz_class-10].set_xticks([0, 1, 2, 3])          
                   axarr[2, enz_class-10].set_xticklabels(['One', 'Three', 'Five', 'Ten'],fontsize=16)
     if num_plots<=6 and not num_plots<=5:
          plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
          #plt.setp([a.get_yticklabels() for a in axarr[:, 0]], visible=False)
          plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
          plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)
     elif num_plots>=10:
          plt.setp([a.get_yticklabels() for a in axarr[:, 3]], visible=False)
          plt.setp([a.get_yticklabels() for a in axarr[:, 4]], visible=False)
     matplotlib.rc('ytick', labelsize=18)
     sns.despine(top=True, right=True)
     plt.tight_layout()
     plt.savefig('OUT_grids/'+nam+'_dpi_500.png', dpi=500)
     plt.show()
     return f

df_AA_sum_mean=deepcopy(df_sum_mean)
df_AA_sum_mean=df_AA_sum_mean[(df_AA_sum_mean.index.str.contains('AA'))]
f=make_barplot_all_index('AA',df_AA_sum_mean)
plt.clf()
plt.close(f)  

df_CE_sum_mean=deepcopy(df_sum_mean)
df_CE_sum_mean=df_CE_sum_mean[(df_CE_sum_mean.index.str.contains('CE'))]
f=make_barplot_all_index('CE',df_CE_sum_mean)
plt.clf()
plt.close(f)  

df_PL_sum_mean=deepcopy(df_sum_mean)
df_PL_sum_mean=df_PL_sum_mean[(df_PL_sum_mean.index.str.contains('PL'))]
f=make_barplot_all_index('PL',df_PL_sum_mean)
plt.clf()
plt.close(f)  

#df_cellulosome_associated_sum_mean=deepcopy(df_sum_mean)
#cellulosome_associated_filter_list=['SLH', 'dockerin', 'cohesin']
#cellulosome_associated_filter_list_crossreferenced=[]
#for x in cellulosome_associated_filter_list:
#    if x in df_cellulosome_associated_sum_mean.index:
#        cellulosome_associated_filter_list_crossreferenced.append(x)    
#for x in cellulosome_associated_filter_list_crossreferenced:
#    df_cellulosome_associated_sum_mean=df_cellulosome_associated_sum_mean[(df_cellulosome_associated_sum_mean.index.str.contains(x))]
#make_barplot_all_index('cell_asoc',df_cellulosome_associated_sum_mean)

##################################################################################################################### 

###################PDF out######################
###TITLE PAGE###########
#canvasname="ISDE_auto_targets_e_val_filter_"+str(e_value_filter)+"_rank_top_n_class_"+str(extract_top_n_class)+'_TS_'+time.strftime('%H%M%S')+'.pdf'
#c=canvas.Canvas(canvasname) 
#c.setFont("Times-Roman",10)
#c.drawCentredString(300,750,"Potential targets")
#c.drawCentredString(550,750,"1")
#c.drawCentredString(300,735,"Overview")
#c.line(50, 730, 550, 730)
#c.drawCentredString(300, 715, "Filtering parameters")
#c.drawString(50, 700, "E-value: "+str(e_value_filter))
#c.drawCentredString(300, 700, "Significant sequences: "+str(sig_seq_filter))
#c.drawCentredString(500, 700, "Phylum percent filter: "+str(phylum_percent_filter))
#c.drawString(50, 690, "Class percent filter: " +str(class_percent_filter))
#c.drawCentredString(300, 690, "Manual target entry regex: "+str(', '.join(manual_target_list)))
#c.drawCentredString(500, 690, "Enzyme class filter: " + str(enzyme_class_filter))
#c.drawString(50, 680, "Top N enzyme classes by total "+ u"\u2211""mean molar percentage" +': '+ str(extract_top_n_class))
#c.drawString(400, 680, "CBM only ORFs included: "+str(manual_entry_of_CBM_only_ORFs))
#c.line(50,675,550,675)
#c.drawCentredString(300, 665, "Breakdown")
#c.drawString(50, 652, "Top " + str(extract_top_n_class) + ' most abundant enzyme classes: ')
#classonline=12
#classloop=0
#starty=652
#for x in list(range(1,math.ceil(len(top_target_list)/classonline))):
#      classloop+=1
#      c.drawString(215, (starty-(10*(classloop-1))), str(' '.join(top_target_list_comma[classonline*(classloop-1):classonline*(classloop)])))
#c.drawString(50, 625, "Total number of targets: " + str(len(target_ORF_dict.keys())))
#c.line(50, 620, 550, 620)
#spacer=spacer=(610-50)/(len(pdf_count_dict)+1)
#countx=1
#c.drawCentredString(60,(620-(spacer*countx)),'Rank')
#c.drawCentredString(100,(620-(spacer*countx)),'Class')
#c.drawCentredString(140,(620-(spacer*countx)),'Targets')
#pagen=2
#for x in pdf_count_dict:
#     countx+=1
#     c.drawCentredString(60,(620-(spacer*countx)),str(x))
#     c.drawCentredString(100,(620-(spacer*countx)),str(pdf_count_dict[x]['class']))
#     c.drawCentredString(140,(620-(spacer*countx)),str(pdf_count_dict[x]['count']))
#heatmap=ImageReader(clustermapstr)
#distribution=ImageReader(evaluestr)
##c.drawImage(img, x, y (of lower left corner or img to start drawing), imgwidth, imgheight, preserveAspectRatio=True)
#c.drawImage(heatmap,200,340, width=300, height=275,preserveAspectRatio=False)
#c.drawImage(distribution, 200, 75, width=300, height=200, preserveAspectRatio=False, mask='auto')
#c.showPage()
###########################top of page
#countloop=0
#for x in target_ORF_dict:
#    c.setFont("Times-Roman",10)
#    countloop+=1
#    if not countloop%2==0:
#         c.drawCentredString(300,770,str(x))
#         c.drawCentredString(550,770,str(pagen))
#         c.drawCentredString(50,770,str(pagen)+"A")
#         c.drawString(50, 750,"Post filter CAZy rank: "+str(target_ORF_dict[x]['ORF_rank']))
#         c.drawCentredString(300, 750, "Class rank: "+str(target_ORF_dict[x]['Class_rank']))
#         c.drawCentredString(480, 750, "Rank within all hits("+str(len(master_nested_dict)) +'): ' + str(target_ORF_dict[x]['totalrank']))
#         c.setFont('Times-Roman',8)    
#         c.drawString(50, 550,"BlastP tophit: " + str(target_ORF_dict[x]['BlastP_tophit']))
#         c.setFont('Times-Roman', 10)
#         c.drawString(50, 735, "Lineage(BlastP): ") 
#         p=(str(target_ORF_dict[x]['lineage']))
#         countk=0
#         for y in list(p.split(';',p.count(';'))):
#            countk+=1
#            c.drawString(140,750-(15*countk), str(y))
#         vector=ImageReader(target_ORF_dict[x]['vector'])
#         graph=ImageReader(target_ORF_dict[x]['graph'])
##c.drawImage(img, x, y (of lower left corner or img to start drawing), imgwidth, imgheight, preserveAspectRatio=True)
#         c.drawImage(vector,50,420, width=480, height=130,preserveAspectRatio=False, mask='auto')
#         c.drawImage(graph, 350, 530, width=250, height=250,preserveAspectRatio=True)
#         c.line(50,422, 550,422)
#    else:
#################################bottom of page
#         c.drawCentredString(300,410,str(x))
#         c.drawCentredString(50,410,str(pagen)+"B")
#         c.drawString(50, 390,"Post filter CAZy rank: "+str(target_ORF_dict[x]['ORF_rank']))
#         c.drawCentredString(300, 390, "Class rank: "+str(target_ORF_dict[x]['Class_rank']))
#         c.drawCentredString(480, 390, "Rank within all hits("+str(len(master_nested_dict)) +'): ' + str(target_ORF_dict[x]['totalrank']))
#         c.setFont('Times-Roman',8)    
#         c.drawString(50, 205,"BlastP tophit: " + str(target_ORF_dict[x]['BlastP_tophit']))
#         c.setFont("Times-Roman",10)
#         c.drawString(50, 375, "Lineage(BlastP): ") 
#         p=(str(target_ORF_dict[x]['lineage']))
#         countk=0
#         for y in list(p.split(';',p.count(';'))):
#              countk+=1
#              c.drawString(140,390-(15*countk), str(y))
#         vector=ImageReader(target_ORF_dict[x]['vector'])
#         graph=ImageReader(target_ORF_dict[x]['graph'])
#         c.drawImage(vector,50,70, width=480, height=130,preserveAspectRatio=False, mask='auto')
#         c.drawImage(graph, 350, 175, width=250, height=250,preserveAspectRatio=True)
#         pagen+=1
#         c.showPage()
#c.save()








###################
CAZy_CBM_correlation={}
for x in list(df_sum_mean.index):
    CAZy_CBM_correlation[x]={}
    for y in list(df_sum_mean_CBM.index):
        CAZy_CBM_correlation[x][y]=0
        CAZy_CBM_correlation[x]['None']=0
        
for CAZy_family in CAZy_CBM_correlation:
    for ORF_index in list(df_mean.index):
        if bool(re.search(CAZy_family, ORF_index))==True:
            for CBM in CAZy_CBM_correlation[CAZy_family]:
                 if bool(re.search(str(CBM), ORF_index))==True:
                      CAZy_CBM_correlation[CAZy_family][CBM]=CAZy_CBM_correlation[CAZy_family][CBM]+len(re.findall(CBM, ORF_index))
                      
df_CAZy_CBM_correlation=pd.DataFrame.from_dict(CAZy_CBM_correlation, orient='index')
df_CAZy_CBM_correlation=df_CAZy_CBM_correlation.T

df_CAZy_CBM_correlation.loc['None'][df_CAZy_CBM_correlation.sum(axis=0)==0]=1
                           
cell_color='light magenta'
hemicell_color="bright sky blue"
chitin_color='seafoam'
    
cellulose_associated_colours=sns.xkcd_palette(['white',cell_color])
hemicellulose_associated_colours=sns.xkcd_palette(['white',hemicell_color])
chitin_associated_colours=sns.xkcd_palette(['white',chitin_color])

cbm_cell_list=[]
cbm_hemi_list=[]
cbm_chitin_list=[]              
              
for x in df_CAZy_CBM_correlation.index:
     if x == 'None':
         cbm_cell_list.append(cellulose_associated_colours[0])
         cbm_hemi_list.append(hemicellulose_associated_colours[0])
         cbm_chitin_list.append(chitin_associated_colours[0])
     if not x == 'None':
          if cbm_pickle_cellulose[x]==1:
              cbm_cell_list.append(cellulose_associated_colours[1])
          if cbm_pickle_cellulose[x]==0:
              cbm_cell_list.append(cellulose_associated_colours[0])
          if cbm_pickle_hemicellulose[x]==1:
              cbm_hemi_list.append(hemicellulose_associated_colours[1])
          if cbm_pickle_hemicellulose[x]==0:
              cbm_hemi_list.append(hemicellulose_associated_colours[0])
          if cbm_pickle_chitin[x]==1:
              cbm_chitin_list.append(chitin_associated_colours[1])
          if cbm_pickle_chitin[x]==0:
             cbm_chitin_list.append(chitin_associated_colours[0])
 
cell_row_series=pd.Series(cbm_cell_list, index=df_CAZy_CBM_correlation.index, name='') #Cellulose
hemicell_row_series=pd.Series(cbm_hemi_list, index=df_CAZy_CBM_correlation.index, name='') #Hemicellulose
chitin_row_series=pd.Series(cbm_chitin_list, index=df_CAZy_CBM_correlation.index, name='') #Chitin
  
#CBM_master_row_col_dict=pd.DataFrame(dict(Cellulose=cell_row_series, Hemicellulose=hemicell_row_series, Chitin=chitin_row_series))

CBM_master_row_col_handles=[]
CBM_master_row_col_handles.append(patches.Patch(color=sns.xkcd_rgb[cell_color], label='Cellulose'))
CBM_master_row_col_handles.append(patches.Patch(color=sns.xkcd_rgb[hemicell_color], label='Hemicellulose'))
CBM_master_row_col_handles.append(patches.Patch(color=sns.xkcd_rgb[chitin_color], label='Chitin'))

sns.set(font_scale=1.3)
#0 = normalised across the CBMs
row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_CAZy_CBM_correlation, df_CAZy_CBM_correlation.T))
cbm_cor = sns.clustermap(df_CAZy_CBM_correlation.T, row_linkage=row_linkage, col_linkage=col_linkage,col_cluster=False, row_cluster=False, figsize=(7,12), standard_scale=0, cmap="Blues", col_colors=[cell_row_series, hemicell_row_series, chitin_row_series])
plt.setp(cbm_cor.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=12)
cbm_cor.ax_heatmap.set(xlabel = 'Family', ylabel='Normalised count')
cbm_cor.ax_heatmap.legend(loc=(0.1,1.075), handles=CBM_master_row_col_handles, title='Known binding associations', ncol=5)
cbm_cor.cax.set_position([1.05, .2, .03, .45]) 
sns.clustermap
cbmcorclustermapstr='OUT_Clustermaps/CAZy_CBM_corelation_NOTnormalised'+time.strftime('%H%M%S') + '_500_dpi.png' 
cbm_cor.savefig(cbmcorclustermapstr, dpi=500)    

sns.set(font_scale=1.5)
#0 = normalised across the CBMs
#row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (df_CAZy_CBM_correlation, df_CAZy_CBM_correlation.T))
cbm_cor = sns.clustermap(df_CAZy_CBM_correlation.T, col_cluster=False, row_cluster=False, figsize=(7,12), cmap="Blues", col_colors=[cell_row_series, hemicell_row_series, chitin_row_series])
plt.setp(cbm_cor.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, size=12)
cbm_cor.ax_heatmap.set(xlabel = 'Family', ylabel='Associations (n)')
cbm_cor.ax_heatmap.legend(loc=(0.1,1.075), handles=CBM_master_row_col_handles, title='Known binding associations', ncol=5)
cbm_cor.cax.set_position([1.05, .2, .03, .45]) 
sns.clustermap 
cbmcorclustermapstr='OUT_Clustermaps/CAZy_CBM_corelation_normalised'+time.strftime('%H%M%S') + '_500_dpi.png' 
cbm_cor.savefig(cbmcorclustermapstr, dpi=500)   
plt.clf()
plt.close()   
####################GENERAL SEARCH QUERIES##########################
# for x in df_master_mean[df_master_mean['rank']<500].index:
#    for y in crossref_molarpct:
#        if bool(re.search(x, y))==True:
#             print(y)     
################################### 
#for k in crossref_molarpct:
#    ki=((''.join(reversed((''.join(reversed(k))).split(' ', 1)[0])).strip('()')))
#    for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
#          temp_seq=str(seq.seq)
#          if re.search(r'SSSSSSS',temp_seq):
#            if ki == seq.id:
#                print(k+'\n'+temp_seq)
##########################################
#temp={}
#temple=[]
#for x in df_master_mean[df_master_mean['rank']<2000].index:
#    for y in e_value_filter_keys_named:
#         if bool(re.search(x, y))==True:
#              temp[x]=y
#for x in temp:
#    if x not in [x.split(' ',x.count(' '))[x.count(' ')].strip('(').strip(')') for x in crossref_molarpct]:
#            for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
#                            if seq.id == x:
#                                print(temp[x] +'\n' + seq.seq +'\n')
##############################################################
#listx=[]
#list_w_interesting_DUFs={'DUF1100':'hydrolase','DUF1573':'PL3/locatonase','DUF4382':'PKD','DUF4957':'Aliginate lyase/Pentin esterase','DUF772':'GH2'}
#for x in pfam_DUF_dict:
#    if x in df_master_mean[df_master_mean['rank']<8000].index:
#        for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
#                        if seq.id == x and len(seq.seq)>=300 and str(pfam_DUF_dict[x]) in list_w_interesting_DUFs:
#                               print(x +' : ' + str(pfam_DUF_dict[x]) + ' associated with ' + list_w_interesting_DUFs[pfam_DUF_dict[x]] + '\n' + seq.seq +'\n')
#                               listx.append(pfam_DUF_dict[x])
###################################################################
#listx=[]
##list_w_interesting_DUFs={'DUF1100':'hydrolase','DUF1573':'PL3/locatonase','DUF4382':'PKD','DUF4957':'Aliginate lyase/Pentin esterase','DUF772':'GH2'}
#for x in pfam_DUF_dict:
#    if x in df_master_mean[df_master_mean['rank']<8000].index:
#        for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
#                        if seq.id == x and len(seq.seq)>=300 and pfam_DUF_dict[x]=='peroxidase': #and str(pfam_DUF_dict[x]) in list_w_interesting_DUFs:
#                               print(x +' : ' + str(pfam_DUF_dict[x]) +  '\n' + seq.seq +'\n')
#                               listx.append(pfam_DUF_dict[x])
##########################################
############TOP 1000 SEQUENCES####
for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
    seq.description='placeholder'
    seqs[seq.id]=seq
   
df_master_top_1000=df_master_mean[df_master_mean['rank']<1001]
df_master_top_1000=df_master_top_1000.sort_values(by=['rank'])
description_dict={}
top_1000=[]
for x in df_master_top_1000.index:
    description_dict[x]={}
    description_dict[x]['rank']=df_master_top_1000.loc[x]['rank']
    description_dict[x]['full']='rank: ' + str(df_master_top_1000.loc[x]['rank'])
    if x in annotation_dict:
           keylist = [x.split('_', 1)[0] for x in list(annotation_dict[x].keys())]
           cazy=', '.join(list(reversed(sorted(keylist))))
           description_dict[x]['cazy']=cazy
           description_dict[x]['full']=description_dict[x]['full']+', CAZY: ' + str(cazy)                
    if x not in annotation_dict:
        description_dict[x]['cazy']='NA'
        description_dict[x]['full']=description_dict[x]['full']+', CAZY: NA'
    if x in pfam_master_dict:
        for k, v in pfam_master_dict[x].items():
            description_dict[x]['pfam']=str(k.rsplit('_',1)[0])
            description_dict[x]['full']=description_dict[x]['full']+', Pfam: ' + str(k.rsplit('_',1)[0])
    if x not in pfam_master_dict:
        description_dict[x]['pfam']='NA'
        description_dict[x]['full']=description_dict[x]['full']+', Pfam: NA'
    if x in taxa_master_dict:
        description_dict[x]['phylum']=taxa_master_dict[x]['Taxonomy']['phylum']
        description_dict[x]['full']=description_dict[x]['full']+', Phylum: ' + str(taxa_master_dict[x]['Taxonomy']['phylum'])
    if x not in taxa_master_dict:
        description_dict[x]['pfam']='NA'
        description_dict[x]['full']=description_dict[x]['full']+', Phylum: NA'
    print(description_dict[x]['full'])
    if x in seqs:
        seqs[x].description=str(description_dict[x]['full'])
        top_1000.append(seqs[x])
        print(x)

SeqIO.write(top_1000, 'ISDE_top_1000_proteins_translated.fna', 'fasta')
    
for x in description_dict:
    if not description_dict[x]['cazy'] == 'NA':
        print(x, description_dict[x]['full'])




templist=[]
for x in crossref_molarpct:
    if re.search (r'AA2 ', x) or re.search(r'AA2,',x):
        templist.append(x)        
for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
    for x in list(range(0,len(templist))):
          if seq.id==((''.join(reversed((''.join(reversed(templist[x]))).split(' ', 1)[0])).strip('()'))):
             if len(seq.seq)>=300:
                 print(templist[x]+', '+seq.id + ':Length= '+str(len(seq.seq))+'\n'+seq.seq)
        
def give_me_target_seqs(input):
    for seq in SeqIO.parse('R:/rsrch/ncb/lab/Dan/RNAseq/MASCOT_seaches/Dereplicated_ORf_files/Derep_longest_100k_blocks/cat_ORFs_translated_gt300_from_gt500.fna', 'fasta'):
        seqname='_'.join(str(input).split('_',-1)[0:len(str(input).split('_',-1))-1])
        if seq.id==seqname:
            print(seq.id + ':Length= '+str(len(seq.seq))+'\n'+seq.seq)
            global tran
            tran=str(seq.seq)
    for seq in SeqIO.parse('R:/rsrch/ncb/lab/Dan/RNAseq/MASCOT_seaches/Dereplicated_ORf_files/Derep_longest_100k_blocks/Trinity_DL_gt500_t.fasta', 'fasta'):
        if seq.id==input[:-2]:
            print(seq.id + ':Length= '+str(len(seq.seq))+'\n'+seq.seq)
            global nuc
            nuc=str(seq.seq)
            break
    return tran, nuc             

def get_me_start(startregex, endregex, nuc, arg):
    replacements={'A':'T', 'T':'A', 'C':'G', 'G': 'C'}
    if arg.upper()=='YES':
        nuc=''.join(list(reversed(nuc)))
        nuc2=''
        for x in nuc:
            l=replacements[x.upper()]
            nuc2=nuc2+l
        nuc=nuc2
        nuc2=''
    global sequence
    sequence=nuc.upper().split(startregex.upper(),1)[1]
    sequence=startregex.upper()+sequence
    sequence=sequence.upper().split(endregex.upper(),1)[0]
    sequence=sequence+endregex.upper()
    print(str(len(sequence))+'\n')                   
    return sequence

#NB/if call returns list index not in range use nuc.reversed()
#import xlrd
#orflist=[]
#orfdict={}
#for k in trans_annotation_dict_named:
#    if re.search (r'GH9 ', k) or re.search(r'GH9,',k):
#        name=((''.join(reversed((''.join(reversed(k))).split(' ', 1)[0])).strip('()')))
#        name2='_'.join(name.split('_',-1)[0:len(name.split('_',-1))-1])
#        orflist.append(name2) 
#        orfdict[name2]={}
#for seq in SeqIO.parse('R:/rsrch/ncb/lab/Dan/RNAseq/MASCOT_seaches/Dereplicated_ORf_files/Derep_longest_100k_blocks/cat_ORFs_translated_gt300_from_gt500.fna', 'fasta'):
#       seqname='_'.join(str(seq.id).split('_',-1)[0:len(str(seq.id).split('_',-1))-1])
#       if seqname in orflist:
#            orfdict[seqname]['tran']=str(seq.seq)
#print('stage_1_complete')
##orflist2=[x[:-2] for x in orflist]
#countx=0
#for seq in SeqIO.parse('R:/rsrch/ncb/lab/Dan/RNAseq/MASCOT_seaches/Dereplicated_ORf_files/Derep_longest_100k_blocks/Trinity_DL_gt500_t.fasta', 'fasta'):
#       if seq.id in orflist:
#             countx+=1
#             orfdict[seq.id]['nuc']=str(seq.seq)
#             if countx%5==0:print(str(countx) + ' sequences stripped')
#             matched.append(seq.id)
#print('stage_2_complete')
#import xlwt
#book=xlwt.Workbook(encoding="utf-8")
#sheetx=book.add_sheet("sheet1")
#countk=0
#for x in orfdict:
#     countk+=1
#     sheetx.write(countk-1,0,str(x))
#     sheetx.write(countk-1,1,orfdict[x]['tran'])
#     sheetx.write(countk-1,2,orfdict[x]['nuc']) 
#book.save('sequences_GH9_tran_nuc_all.xls')



###superkingdom stuff############# TOTAL PROTEOME ########################
superkingdom_list=[]
superkingdom_dict={}
for x in taxa_master_dict:
    if taxa_master_dict[x]['Taxonomy']['superkingdom'] not in superkingdom_list:
        superkingdom_list.append(taxa_master_dict[x]['Taxonomy']['superkingdom'])
#Add superkingdoms to dict
for x in superkingdom_list:
    if x not in superkingdom_dict:
        superkingdom_dict[x]={}
for x in master_nested_dict:
     for y in master_nested_dict[x]:
         for z in superkingdom_dict:
             if y not in superkingdom_dict[z]:
                  superkingdom_dict[z][y]=0
for x in master_nested_dict:
    for x2 in master_nested_dict[x]:
         if x in taxa_master_dict:
             superkingdom_dict[taxa_master_dict[x]['Taxonomy']['superkingdom']][x2]=superkingdom_dict[taxa_master_dict[x]['Taxonomy']['superkingdom']][x2]+master_nested_dict[x][x2]
         else:   
             superkingdom_dict['NA'][x2]=superkingdom_dict['NA'][x2]+master_nested_dict[x][x2]
 
df_superkingdom=pd.DataFrame.from_dict(superkingdom_dict, orient='index')
df_mean_superkingdom=pd.DataFrame()
df_mean_superkingdom['One']=(df_superkingdom['W1_C201'] + df_superkingdom['W1_3_C429'] + df_superkingdom['W1_4_C429'])/3
df_mean_superkingdom['Three']=(df_superkingdom['W3_C201'] + df_superkingdom['W3_2_C429'] + df_superkingdom['W3_1_C429'])/3
df_mean_superkingdom['Five']=(df_superkingdom['W5_C201'] + df_superkingdom['W5_25_C429'] + df_superkingdom['W5_100_C429'])/3
df_mean_superkingdom['Ten']=(df_superkingdom['W10_C201'] + df_superkingdom['W10_100_i_C429'] + df_superkingdom['W10_100_ii_C429'])/3
 
df_mean_superkingdom['One_err']=df_superkingdom[['W1_C201','W1_3_C429', 'W1_4_C429']].std(axis=1)
df_mean_superkingdom['Three_err']=df_superkingdom[['W3_C201','W3_2_C429','W3_1_C429']].std(axis=1)
df_mean_superkingdom['Five_err']=df_superkingdom[['W5_C201','W5_25_C429','W5_100_C429']].std(axis=1)
df_mean_superkingdom['Ten_err']=df_superkingdom[['W10_C201','W10_100_i_C429','W10_100_ii_C429']].std(axis=1) 

fig, ax = plt.subplots(figsize=(8,6))
sns.set(style="white")
sns.despine(top=True, right=True)
N = 4
a = df_mean_superkingdom.loc['Bacteria'][:4]
a_err = df_mean_superkingdom.loc['Bacteria'][4:]
b = df_mean_superkingdom.loc['Eukaryota'][:4]
b_err = df_mean_superkingdom.loc['Eukaryota'][4:]
ab=a+b
c = df_mean_superkingdom.loc['NA'][:4]
c_err = df_mean_superkingdom.loc['NA'][4:]
abc=ab+c
d = df_mean_superkingdom.loc['Viruses'][:4]
d_err = df_mean_superkingdom.loc['Viruses'][4:]
abcd=abc+d
e = df_mean_superkingdom.loc['Archaea'][:4]
e_err = df_mean_superkingdom.loc['Archaea'][4:]
abcde=abcd+e
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

errorbar_format=dict(ecolor='black', lw=1, capsize=3, capthick=1)
col=4
colo=sns.xkcd_palette(["clear blue"])+sns.hls_palette(col)
p1 = ax.bar(ind, a, width, yerr=a_err,error_kw=errorbar_format, color=colo[0])
p2 = ax.bar(ind, b, width, bottom=a, yerr=b_err,error_kw=errorbar_format, color=colo[1])
p3=ax.bar(ind, c, width, bottom=ab, yerr=c_err,error_kw=errorbar_format, color=colo[2])
p4=ax.bar(ind, d, width, bottom=abc, yerr=d_err,error_kw=errorbar_format, color=colo[3])
p5=ax.bar(ind, e, width, bottom=abcd, yerr=c_err,error_kw=errorbar_format, color=colo[4])
plt.ylabel('Secretome molar percentage', fontsize=16)
#plt.title('Scores by group and gender')
ax.tick_params(direction='out', length=3, width=2)
plt.yticks(fontsize=16)
plt.xticks(ind, ('One', 'Three', 'Five', 'Ten'), fontsize=16)
plt.xlabel('Week', fontsize=18)
#plt.yticks(np.arange(0, 81, 10))
plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), (['Bacteria', 'Eukaryota', 'NA', 'Viruses', 'Archaea']), ncol=5, loc=(-0.15,-.25), fontsize=14)
fig.subplots_adjust(bottom=0.2)
plt.show()   
fig.savefig('OUT_Taxonomy_superkingdom_graphs/Stacked_barplot_superkingdom_stderror_500_dpi.png', dpi=500)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all') 
#####################CAZY FRACTION ONLY ######################
superkingdom_list=[]
superkingdom_dict={}
for x in taxa_master_dict:
    if taxa_master_dict[x]['Taxonomy']['superkingdom'] not in superkingdom_list:
        superkingdom_list.append(taxa_master_dict[x]['Taxonomy']['superkingdom'])
#Add superkingdoms to dict
for x in superkingdom_list:
    if x not in superkingdom_dict:
        superkingdom_dict[x]={}
for x in crossref_nested_dict:
     for y in crossref_nested_dict[x]:
         for z in superkingdom_dict:
             if y not in superkingdom_dict[z]:
                  superkingdom_dict[z][y]=0
for x in crossref_nested_dict:
    for x2 in crossref_nested_dict[x]:
         if x in taxa_master_dict:
             superkingdom_dict[taxa_master_dict[x]['Taxonomy']['superkingdom']][x2]=superkingdom_dict[taxa_master_dict[x]['Taxonomy']['superkingdom']][x2]+crossref_nested_dict[x][x2]
         else:   
             superkingdom_dict['NA'][x2]=superkingdom_dict['NA'][x2]+crossref_nested_dict[x][x2]
 
df_superkingdom=pd.DataFrame.from_dict(superkingdom_dict, orient='index')
df_mean_superkingdom=pd.DataFrame()

#####Normalise to percentages
for x in list(df_superkingdom.columns):
    df_superkingdom[x]=df_superkingdom[x]/(df_superkingdom[x].sum())*100


df_mean_superkingdom['One']=(df_superkingdom['W1_C201'] + df_superkingdom['W1_3_C429'] + df_superkingdom['W1_4_C429'])/3
df_mean_superkingdom['Three']=(df_superkingdom['W3_C201'] + df_superkingdom['W3_2_C429'] + df_superkingdom['W3_1_C429'])/3
df_mean_superkingdom['Five']=(df_superkingdom['W5_C201'] + df_superkingdom['W5_25_C429'] + df_superkingdom['W5_100_C429'])/3
df_mean_superkingdom['Ten']=(df_superkingdom['W10_C201'] + df_superkingdom['W10_100_i_C429'] + df_superkingdom['W10_100_ii_C429'])/3
 
                  
df_mean_superkingdom['One_err']=df_superkingdom[['W1_C201','W1_3_C429', 'W1_4_C429']].std(axis=1)
df_mean_superkingdom['Three_err']=df_superkingdom[['W3_C201','W3_2_C429','W3_1_C429']].std(axis=1)
df_mean_superkingdom['Five_err']=df_superkingdom[['W5_C201','W5_25_C429','W5_100_C429']].std(axis=1)
df_mean_superkingdom['Ten_err']=df_superkingdom[['W10_C201','W10_100_i_C429','W10_100_ii_C429']].std(axis=1) 

fig, ax = plt.subplots(figsize=(8,6))
sns.set(style="white")
sns.despine(top=True, right=True)
N = 4
a = df_mean_superkingdom.loc['Bacteria'][:4]
a_err = df_mean_superkingdom.loc['Bacteria'][4:]
b = df_mean_superkingdom.loc['Eukaryota'][:4]
b_err = df_mean_superkingdom.loc['Eukaryota'][4:]
ab=a+b
c = df_mean_superkingdom.loc['NA'][:4]
c_err = df_mean_superkingdom.loc['NA'][4:]
abc=ab+c
d = df_mean_superkingdom.loc['Viruses'][:4]
d_err = df_mean_superkingdom.loc['Viruses'][4:]
abcd=abc+d
e = df_mean_superkingdom.loc['Archaea'][:4]
e_err = df_mean_superkingdom.loc['Archaea'][4:]
abcde=abcd+e
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

errorbar_format=dict(ecolor='black', lw=1, capsize=3, capthick=1)
col=4
colo=sns.xkcd_palette(["clear blue", "coral"])+sns.hls_palette(col)
p1 = ax.bar(ind, a, width, yerr=a_err,error_kw=errorbar_format, color=colo[0])
p2 = ax.bar(ind, b, width, bottom=a, yerr=b_err,error_kw=errorbar_format, color=colo[1])
p3=ax.bar(ind, c, width, bottom=ab, yerr=c_err,error_kw=errorbar_format, color=colo[2])
p4=ax.bar(ind, d, width, bottom=abc, yerr=d_err,error_kw=errorbar_format, color=colo[3])
p5=ax.bar(ind, e, width, bottom=abcd, yerr=c_err,error_kw=errorbar_format, color=colo[4])
plt.ylabel('Total CAZy (%)', fontsize=16)
#plt.title('Scores by group and gender')
ax.tick_params(direction='out', length=3, width=2)
plt.yticks(fontsize=16)
plt.xticks(ind, ('One', 'Three', 'Five', 'Ten'), fontsize=16)
plt.xlabel('Week', fontsize=18)
#plt.yticks(np.arange(0, 81, 10))
plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), (['Bacteria', 'Eukaryota', 'NA', 'Viruses', 'Archaea']), ncol=5, loc=(-0.15,-.25), fontsize=14)
fig.subplots_adjust(bottom=0.2)
plt.show()  
fig.savefig('OUT_Taxonomy_superkingdom_graphs/Stacked_barplot_superkingdom_CAZy_stderror_500_dpi.png', dpi=500)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all') 



w1plot=list(df_sum_activity_taxa_dict_molpct_w1.astype(bool).sum(axis=0))
w3plot=list(df_sum_activity_taxa_dict_molpct_w3.astype(bool).sum(axis=0))
w5plot=list(df_sum_activity_taxa_dict_molpct_w5.astype(bool).sum(axis=0))
w10plot=list(df_sum_activity_taxa_dict_molpct_w10.astype(bool).sum(axis=0))

classcount=w1plot+w3plot+w5plot+w10plot
namlist=['One']*int(len(w1plot))+['Three']*int(len(w3plot))+['Five']*int(len(w5plot))+['Ten']*int(len(w10plot))
df_violoin=pd.DataFrame({'Diversity':classcount,'Week': namlist})         
columns=list(df_sum_activity_taxa_dict_molpct_w1.columns)
starlabels=['Proteobacteria', 'Bacteroidetes']
starcol=['blue', 'red', 'cyan', 'yellow']
lencol=len(columns) 
viocol=sns.xkcd_palette(["light grey","clear blue","light grey","clear blue"])       
clcol=sns.hls_palette(lencol, l=.5, s=.8)  
        
fig, ax = plt.subplots(figsize=(8,6))
sns.set(style="white")
ax=sns.violinplot(x="Week", y="Diversity", data=df_violoin, palette=viocol, inner='points')  

starcount=0 
for x in list(range(0, lencol)):
    if columns[x] not in starlabels:
         ax.plot([0,1,2,3], [w1plot[x], w3plot[x], w5plot[x], w10plot[x]], 'o', color=clcol[x], markersize=5, label=str(columns[x]))
    else:
         ax.plot([0,1,2,3], [w1plot[x], w3plot[x], w5plot[x], w10plot[x]], '*', color=starcol[starcount], markersize=15, label=str(columns[x]))
         starcount+=1

means=[np.mean(w1plot),np.mean(w3plot),np.mean(w5plot),np.mean(w10plot)]
ax.plot([0,1,2,3], means, '-', color='black', label='Mean')

plt.legend(loc=(0, .9), fontsize=12, ncol=3)
sns.despine(top=True, right=True)  
ax.tick_params(direction='out', length=3, width=2)
plt.yticks(fontsize=16)
plt.xticks(ind, ('One', 'Three', 'Five', 'Ten'), fontsize=16)
plt.xlabel('Week', fontsize=18)
plt.ylabel('CAZy (n)', fontsize=18)
fig.subplots_adjust(top=0.825)
fig.savefig("OUT_Taxonomy_phylum_graphs/violoinplot_cazyn.png", dpi=500)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')  
#######class level

w1plot=list(df_sum_activity_taxa_dict_molpct_w1_class.astype(bool).sum(axis=0))
w3plot=list(df_sum_activity_taxa_dict_molpct_w3_class.astype(bool).sum(axis=0))
w5plot=list(df_sum_activity_taxa_dict_molpct_w5_class.astype(bool).sum(axis=0))
w10plot=list(df_sum_activity_taxa_dict_molpct_w10_class.astype(bool).sum(axis=0))
classcount=w1plot+w3plot+w5plot+w10plot
namlist=['One']*len(w1plot)+['Three']*len(w1plot)+['Five']*len(w1plot)+['Ten']*len(w1plot)

columns=list(df_sum_activity_taxa_dict_molpct_w1_class.columns)
starlabels=['Gammaproteobacteria', 'Bacteroidia', 'Cytophagia', 'Deltaproteobacteria', 'Flavobacteria']
starcol=['blue', 'red', 'cyan', 'yellow', 'magenta']
lencol=len(columns) 
viocol=sns.xkcd_palette(["light grey","clear blue","light grey","clear blue"])       
clcol=sns.hls_palette(lencol, l=.5, s=.8)       
df_violoin=pd.DataFrame({'Diversity':classcount,'Week': namlist}) 
fig, ax = plt.subplots(figsize=(8,6))
sns.set(style="white")
ax=sns.violinplot(x="Week", y="Diversity", data=df_violoin, palette=viocol, inner='points') 

means=[np.mean(w1plot),np.mean(w3plot),np.mean(w5plot),np.mean(w10plot)]
ax.plot([0,1,2,3], means, '-', color='black', label='Mean')

starcount=0 
for x in list(range(0, lencol)):
    if columns[x] not in starlabels:
         ax.plot([0,1,2,3], [w1plot[x], w3plot[x], w5plot[x], w10plot[x]], 'o', color=clcol[x], markersize=5, label=str(columns[x]))
    else:
         ax.plot([0,1,2,3], [w1plot[x], w3plot[x], w5plot[x], w10plot[x]], '*', color=starcol[starcount], markersize=15, label=str(columns[x]))
         starcount+=1
#ax.plot([0,1,2,3],[19,25,28,27],'*', color='blue', markersize=15, label='Proteobacteria')
#ax.plot([0,1,2,3],[20,18,20,24],'*', color='red', markersize=15, label='Bacteroidetes')
plt.legend(loc=(0, .9), fontsize=12, ncol=3)
sns.despine(top=True, right=True)  
ax.tick_params(direction='out', length=3, width=2)
plt.yticks(fontsize=16)
plt.xticks(ind, ('One', 'Three', 'Five', 'Ten'), fontsize=16)
plt.xlabel('Week', fontsize=18)
plt.ylabel('CAZy (n)', fontsize=18)
fig.subplots_adjust(top=0.78)
fig.savefig("OUT_Taxonomy_phylum_graphs/violinplot_class_level.png", dpi=500)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all') 







###############FUNGAL ORFS

fungal_ORF_count=0
fungal_ORFs=[]
for x in taxa_master_dict:
     for y in taxa_master_dict[x]['Taxonomy']: 
          if taxa_master_dict[x]['Taxonomy'][y]=='Fungi':
                   fungal_ORF_count+=1
                   fungal_ORFs.append(x)
                   print(x, y)
                   

fungal_df=df_master_mean.ix[fungal_ORFs]
fungal_df.sum(axis=1)

#################################TEST network graph###########
#################################TEST network graph###########
#################################TEST network graph###########
#################################TEST network graph###########
import networkx as nx


###select main df here
subjectdf_phylum=(df_sum_activity_taxa_dict_molpct_w1+df_sum_activity_taxa_dict_molpct_w3+df_sum_activity_taxa_dict_molpct_w5+df_sum_activity_taxa_dict_molpct_w10)/4
subjectdf_class=(df_sum_activity_taxa_dict_molpct_w1_class+df_sum_activity_taxa_dict_molpct_w3_class+df_sum_activity_taxa_dict_molpct_w5_class+df_sum_activity_taxa_dict_molpct_w10_class)/4
subjectdf_3=(df_sum_activity_taxa_dict_molpct_w3_class+df_sum_activity_taxa_dict_molpct_w5_class+df_sum_activity_taxa_dict_molpct_w10_class)/3

#subjectdf=df_sum_activity_taxa_dict_molpct_w1
#subjectdf=df_sum_activity_taxa_dict_total_format
###SELECT INPUT DF HERE###
lvl='class'
subjectdf=subjectdf_class
########################################


#filter raw dataframe by taxa with < X connections or weight (mol%)
taxfiltervalue=0.025
enzfiltervalue=0.00125
#filter tax
df_sum_activity_networkx_filtered=subjectdf[subjectdf.columns[subjectdf.sum(axis=0)>taxfiltervalue]]
#filter enz
df_sum_activity_networkx_filtered=df_sum_activity_networkx_filtered[df_sum_activity_networkx_filtered.sum(axis=1)>enzfiltervalue]

#df_sum_activity_taxa_dict_total_format is test df
#UNSTACK the dataframe, level 0 = taxa, level 1 = enzyme, third and final column should = weight
dfnetwork=df_sum_activity_networkx_filtered.unstack().reset_index()
dfnetwork=dfnetwork.rename(columns={'level_0':'Taxa', 'level_1':'Enzyme', 0:'Weight'})
#Filter non_zero values
dfnetwork=dfnetwork[dfnetwork['Weight']>0]
namemaps={}
changedict={'Gammaproteobacteria':'γ-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria'}
for x in dfnetwork['Taxa']:
    if x in changedict:
        namemaps[x]=changedict[x]
    if x not in changedict:
        namemaps[x]=x
dfnetwork['Taxa']=dfnetwork['Taxa'].map(namemaps)

#delete troublesome nodes
troublesome_node=[]#['GH4']
for x in troublesome_node:
    dfnetwork=dfnetwork[~(dfnetwork['Enzyme']==x)]


#all unqiue taxa
list_taxa=list(dfnetwork['Taxa'].unique())
list_enzyme=list(dfnetwork['Enzyme'].unique())

###enzyme colors same as heatmaps
colors=sns.xkcd_palette(['pale orange', 'coral', "baby blue", 'pastel yellow', 'lilac', 'pale lilac', 'light green'])    
color_dict={'AA':colors[0], 'GH':colors[2], 'CE':colors[1], 'PL':colors[3], 'SLH':colors[5], 'CBM':colors[4]}

#taxa colors use same as circos plots
taxa_color_dict={'Proteobacteria': (36, 122, 253),
 'Bacteroidetes': (254, 44, 84),
 'Verrucomicrobia': (255, 255, 45),
 'Ignavibacteriae': (117, 253, 99),
 'NA': (255, 255, 45),
 'Planctomycetes': (178, 91, 178),
 'Spirochaetes': (150, 150, 150),
 'Chlorophyta': (229, 25, 127),
 'Firmicutes': (229, 25, 229),
 'Actinobacteria': (25, 127, 229),
 'Candidatus Riflebacteria': (229, 25, 25),
 'Tenericutes': (25, 127, 229),
 'Annelida': (25, 229, 229),
 'Lentisphaerae': (127, 229, 25),
 'Candidatus Moranbacteria': (25, 229, 127),
 'γ-proteobacteria': (36, 122, 253),
 'Bacteroidia': (254, 44, 84),
 'Flavobacteriia': (117, 253, 99),
 'δ-proteobacteria': (255, 165, 0),
 'Ignavibacteria': (178, 91, 178),
 'Verrucomicrobiae': (150, 150, 150),
 'Cytophagia': (185, 66, 244),
 'Trebouxiophyceae': (229, 25, 229),
 'Clostridia': (25, 25, 229),
 'Bacilli': (229, 25, 25),
 'Phycisphaerae': (25, 229, 229),
 'Mollicutes': (127, 229, 25),
 'Sphingobacteriia': (25, 229, 127),
 'Saprospiria': (229, 25, 229),
 'Opitutae': (25, 25, 229),
 'Chitinophagia': (229, 25, 25),
 'Clitellata': (25, 127, 229)}

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap





fig=plt.figure(figsize=(8,8))
plt.axis('off')

###drawing taxa inside an inner circle
##assuming centre of circle is 0,0
filtered_taxa_list=['γ-proteobacteria', 'Cytophagia', 'Bacteroidia', 'Flavobacteriia', 'NA', 'Verrucomicrobiae','δ-proteobacteria']
filtered_enzyme_list=list_enzyme

num_taxa=len(filtered_taxa_list)
angles=[360/(num_taxa)*A for A in list(range(0,num_taxa))]
#convert_angles_to_radians
radians=[A*math.pi/180 for A in angles]
#radius =1
r=.75

coord_tuples=[(r*math.cos(A), r*math.sin(A)) for A in radians]

fixed_positions=dict(zip(filtered_taxa_list, coord_tuples))
fixed_nodes=fixed_positions.keys()

###inner circle####
r2=r*0.575
##call G2 as graph to extract nodes
G2=nx.from_pandas_edgelist(dfnetwork, source='Enzyme', target='Taxa', edge_attr=True)
enz_degree=3
enzdeg = G2.degree()
enz_nodes_circle= [n[0] for n in enzdeg if n[0] in filtered_enzyme_list if n[1] >= enz_degree]
manual_order=['CBM6','CBM60', 'GH16', 'GH23', 'CE6', 'GH3', 'GH10', 'CE10', 'GH13', 'GH5', 'GH11', 'GH30','CE11', 'CE1', 'CBM50']
enz_nodes_circle=manual_order+[A for A in enz_nodes_circle if A not in manual_order]
num_enz=len(enz_nodes_circle)
angles_enz=[360/(num_enz)*A for A in list(range(0,num_enz))]
radians_enz=[A*math.pi/180 for A in angles_enz]
coord_tuples_enz=[(r2*math.cos(A), r2*math.sin(A)) for A in radians_enz]
##set order of appearance
############
#add other fixed points for clarity
manual_append_list=['GH20', 'GH65', 'GH31', 'GH8', 'GH103', 'AA7', 'CBM30', 'GH29', 'CBM35', 'GH6', 'CBM2', 'CBM48', 'GH127', 'GH130', 'GH53', 'GH43', 'AA2', 'GH9', 'CBM4', 'CBM9', 'GH109']
manual_append_coords=[(-.7, 0.1),(0.5, 0.3), (-0.6, -0.2), (-0.4, 0.4), (.8, -.75), (.7, -.8), (.55, -.75),(-.4,.9),(.45,.9), (0.95,0.55), (.8, -.45), (0.125,0.65), (-0.4,1.1), (-.3,1.2), (-.1,1.15),(-0.5, 0.5),(0.6, 0.4),(0.5, -0.325),(0.65, -0.35), (0.1,0.55),(0.18722709, 1.01182814)]
enz_nodes_circle=enz_nodes_circle+manual_append_list
coord_tuples_enz=coord_tuples_enz+manual_append_coords

fixed_positions_enz=dict(zip(enz_nodes_circle, coord_tuples_enz))
fixed_nodes_enz=fixed_positions_enz.keys()
####
#append the data from both circles together
aggregated_coord_tuples=coord_tuples+coord_tuples_enz
aggregated_fixed_nodes=list(fixed_nodes)+list(fixed_nodes_enz)
aggregated_fixed_positions=dict(zip(aggregated_fixed_nodes, aggregated_coord_tuples))

#create graph
G=nx.from_pandas_edgelist(dfnetwork, source='Enzyme', target='Taxa', edge_attr=True)
#filter nodes with X or less outgoing edges
tax_edge_filter=5
outdeg = G.degree()
to_remove= [n[0] for n in outdeg if n[0] in list_taxa if n[1] <= tax_edge_filter]
G.remove_nodes_from(to_remove)
##now remove filtered enzyme groups
to_remove_enz= [n[0] for n in outdeg if n[1] == 0]
to_remove_enz=to_remove_enz#+['GH6', 'CBM2']
G.remove_nodes_from(to_remove_enz)
#generate layout for graph
#edit k value for spreading out nodes
filtered_taxa=[tax for tax in list_taxa if tax in G.nodes()]
if lvl=='class':
   # filtered_taxa_list=['Gammaproteobacteria', 'Bacteroidia', 'Flavobacteriia', 'Cytophagia', 'NA', 'Deltaproteobacteria']
     filtered_taxa_list=['γ-proteobacteria', 'Cytophagia', 'Bacteroidia', 'Flavobacteriia', 'NA', 'Verrucomicrobiae','δ-proteobacteria']









[filtered_taxa_list.append(x) for x in filtered_taxa if x not in filtered_taxa_list]
#random.shuffle(filtered_taxa_list)
filtered_enzyme_list=[enz for enz in list_enzyme if enz in G.nodes()]
#make dictionary to label taxa
taxa_label_dict=dict(zip(filtered_taxa_list, filtered_taxa_list))
enzyme_label_dict=dict(zip(filtered_enzyme_list, filtered_enzyme_list))


#layout = nx.spring_layout(G,iterations=50,k=0.2, pos=fixed_positions, fixed=fixed_nodes)
layout = nx.spring_layout(G,iterations=50,k=0.2, pos=aggregated_fixed_positions, fixed=aggregated_fixed_nodes)
#generate weights
#weight of taxa
#dfnetwork.loc[dfnetwork['Taxa']=='Chlorophyta']['Weight'].sum()
tax_multiplier=1000*10
taxa_size=[dfnetwork.loc[dfnetwork['Taxa']==tax]['Weight'].sum()*tax_multiplier for tax in filtered_taxa_list]
taxa_node_size_dict=dict(zip(filtered_taxa_list, taxa_size))
###enzyme stuff
enz_multiplier=1600*5
enzyme_size=[dfnetwork.loc[dfnetwork['Enzyme']==enz]['Weight'].sum()*enz_multiplier for enz in filtered_enzyme_list]
enzyme_node_size_dict=dict(zip(filtered_enzyme_list, enzyme_size))

####
#node colors and shapes
tax_size=[]
tax_shape='h'
tax_node_size=[]
tax_color=[]
enz_size=[]
enz_shape='o'
enz_node_size=[]
enz_color=[]
nodecolors=[]
nodesizes=[]
for node in G:
     if node in taxa_color_dict:
         nodecolors.append(tuple([float(x)/255 for x in taxa_color_dict[node]]))
         nodesizes.append(taxa_node_size_dict[node])
         ####
         tax_size.append(taxa_node_size_dict[node])
         tax_node_size.append(taxa_node_size_dict[node])
         tax_color.append(tuple([float(x)/255 for x in taxa_color_dict[node]]))
     else:
         for start in color_dict:
             if node.startswith(start):
                 nodecolors.append((color_dict[start]))
                 nodesizes.append(enzyme_node_size_dict[node])
                 #####
                 enz_size.append(enzyme_node_size_dict[node])
                 enz_node_size.append(enzyme_node_size_dict[node])
                 enz_color.append((color_dict[start]))



linewidths=[1]*len(G.nodes())
edgecols=['black']*len(G.nodes())

nx.draw_networkx_nodes(G,layout,nodelist=filtered_taxa,node_shape='H', node_color=tax_color, node_size=tax_node_size, linewidths=[2]*len(filtered_taxa), edgecolors=['black']*len(filtered_taxa))
nx.draw_networkx_nodes(G,layout,nodelist=filtered_enzyme_list,node_shape='o', node_color=enz_color, node_size=enz_node_size, linewidths=[1]*len(filtered_enzyme_list), edgecolors=['black']*len(filtered_enzyme_list))
#draw second smalled hex marker infront of other to get white inner outline and make more bold
tax_node_size_inner=[(t-125)*0.9 for t in tax_node_size]
nx.draw_networkx_nodes(G,layout,nodelist=filtered_taxa,node_shape='H', node_color=tax_color, node_size=tax_node_size_inner, linewidths=[2.5]*len(filtered_taxa), edgecolors=['white']*len(filtered_taxa))
tax_node_size_inner2=[(sx-50)*0.9 for sx in tax_node_size_inner]
nx.draw_networkx_nodes(G,layout,nodelist=filtered_taxa,node_shape='H', node_color=tax_color, node_size=tax_node_size_inner2, linewidths=[1]*len(filtered_taxa), edgecolors=['black']*len(filtered_taxa))

#nx.draw_networkx_labels(G, layout, labels=taxa_label_dict, font_size=10, font_weight='bold')
##offset Y coordinates to sit above the node
#duplicate dictionary so nodes remain in position, but text coords can be changed
offset_layout=deepcopy(layout)
divby=100
offset=0.02
ignore_list=['CBM2','GH6']
offset_enzymes=[]
for node in offset_layout:
     if node not in manual_order and node not in filtered_taxa_list and node not in ignore_list:
         size=enzyme_node_size_dict[node]
         #size = cirlce area, therefore radius = 
         radius=math.sqrt(size)/math.pi
         offset_layout[node]=np.array([offset_layout[node][0], offset_layout[node][1]+(radius/divby)+offset])
         offset_enzymes.append(node)

####make labels for inner circle on the outside
#r3=0.28
ignore_list=['CE1','GH5','GH3']
#coord_tuples_enz=[(r3*math.cos(A), r3*math.sin(A)) for A in radians_enz]
innerlabeldict=dict(zip(manual_order, coord_tuples_enz))
for node in innerlabeldict:
    if node not in ignore_list:
       size=enzyme_node_size_dict[node]
       #size = cirlce area, therefore radius = 
       radius=math.sqrt(size)/math.pi
       offset_layout[node]=np.array([offset_layout[node][0], offset_layout[node][1]+(radius/divby)+offset])
       offset_enzymes.append(node)


#nx.draw_networkx_labels(G, layout, labels=enzyme_label_dict, font_size=10, font_weight='bold')
nx.draw_networkx_labels(G, offset_layout, labels=enzyme_label_dict, font_size=10, font_weight='bold')

#scale edges
actual_edge_weights=[]
edge_weights=[]
edgeweightmultiplier=15*4
edgelist=G.edges()
for x in edgelist:
    if x[0] in filtered_taxa_list:
        tax=x[0]
        enz=x[1]
    if x[1] in filtered_taxa_list:
        tax=x[1]
        enz=x[0]
    weight=dfnetwork[(dfnetwork['Taxa']==tax) & (dfnetwork['Enzyme']==enz)]['Weight'].values[0]
    edge_weights.append(weight*edgeweightmultiplier)
    actual_edge_weights.append(weight)



minimum_edge_weight=0.4
edge_weights=[x if x>=minimum_edge_weight else minimum_edge_weight for x in edge_weights]
####custom cmap####
orig_cmap = matplotlib.cm.coolwarm
orig_cmap= matplotlib.colors.LinearSegmentedColormap.from_list("", ["blue", 'violet','red'])
shifted_cmap = shiftedColorMap(orig_cmap, midpoint=0.1, name='shifted')
###make edges conform to heatmap
#cmap = matplotlib.cm.get_cmap('coolwarm')#'bwr'. 'seismic', 'RdBu'
cmap = shifted_cmap
norm = matplotlib.colors.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))  
norm_im=  matplotlib.colors.Normalize(vmin=min(actual_edge_weights), vmax=max(actual_edge_weights))  
#normalise the value (norm(x)), then convert to rgb with cmap, cmap(norm(x))  
edge_weighted_colors=[cmap(norm(w)) for w in edge_weights]

nx.draw_networkx_edges(G, layout,edgelist=G.edges(), width=edge_weights, edge_color=edge_weighted_colors)
plt.text(-.875, 1.1, ''u"\u2211"r' $\bar{x}$'" mol%", fontsize=12,weight='bold', horizontalalignment='center', verticalalignment='center')
##circle
alphaval=0.2
circle_ax=fig.gca()
circle1 = plt.Circle((0, 0), r, color='#cccccc', alpha=alphaval, zorder=0)
circle2 = plt.Circle((0, 0), r2, color='#cccccc', alpha=alphaval, zorder=0)
circle_ax.add_artist(circle1)
circle_ax.add_artist(circle2)
######
#####draw elipses, pull taxa, pull edges, if edges connections only = 1, rip coords, draw cloud
cloud_dict={}
cloud_dict_formatted={}
for tax in filtered_taxa_list:
    cloud_dict[tax]={}
    cloud_dict_formatted[tax]={}
    nodezzz=[]
    coordzzz=[]
    nodezzz.append(tax)
    coordzzz.append(layout[tax])
    for node in G.edges(tax):
        if G.degree(node[1])<=3:
            cloud_dict_formatted[tax][str(node[1])]=layout[node[1]]
            nodezzz.append(str(node[1]))
            coordzzz.append(layout[node[1]])
    cloud_dict_formatted[tax][str(tax)]=layout[tax]
    cloud_dict[tax]['node']=nodezzz
    cloud_dict[tax]['coord']=coordzzz



draw_order={'γ-proteobacteria':['γ-proteobacteria','GH9','CBM4','CBM2','CE4','CBM5','PL4','GH62','CBM3','CBM10','PL1','CBM44','CBM13','GH6', 'γ-proteobacteria'],\
            'δ-proteobacteria':['δ-proteobacteria','GH103', 'AA7', 'CBM30','δ-proteobacteria'],\
            'Cytophagia':['Cytophagia','CBM37','GH74','CBM48','CBM9','GH109','Cytophagia'],\
            'Bacteroidia':['Bacteroidia','GH109', 'CBM9','CBM48','GH53','GH130','GH127','GH29','Bacteroidia'],\
            'Flavobacteriia':['Flavobacteriia','GH115','GH78', 'GH26','Flavobacteriia'],\
            'NA':['NA','CBM56','GH133', 'GH27','GH18','CE9','SLH','NA']}

draw_order_coords={}
for tax in draw_order:
    templist=[]
    for node in draw_order[tax]:
        templist.append(list(cloud_dict_formatted[tax][node]))
    draw_order_coords[tax]=templist



draw_order_coords={'γ-proteobacteria': [[0.75, 0.0],  [0.8, -0.45],
  [1.0549701529936888, -0.532163394573219],
  [1.2125696430063262, -0.47829814882348864],
  [1.333339362194685, -0.36472957765966796],
  [1.388147497352032, -0.19528240296518512],
  [1.4244153640728336, -0.059095258814673805],
  [1.422447706010992, 0.09186325228924613],
  [1.3485625629515028, 0.21645022545290143],
  [1.304263487657538, 0.34724169998398224],
  [1.1889012553592713, 0.46392495425605723],
  [0.95, 0.55],
  [0.75, 0.0]],
 'δ-proteobacteria': [[0.46761735139405003, -0.5863736118510224],
  [0.8, -0.75],
  [0.7, -0.8],
  [0.55, -0.75],
  [0.46761735139405003, -0.5863736118510224]],
 'Cytophagia': [[0.4676173513940502, 0.5863736118510223],
  [0.7977806951651315, 1.0784102794248869],
  [0.6368174684485154, 1.1725314944346223],
  [0.18722709, 1.01182814],
  [0.4676173513940502, 0.5863736118510223]],
 'Bacteroidia': [[-0.16689070046723575, 0.7311959341363677],
  [0.18722709, 1.01182814],
  [-0.1, 1.15],
  [-0.3, 1.2],
  [-0.4, 1.1],
  [-0.4, 0.9],
  [-0.16689070046723575, 0.7311959341363677]],
 'Flavobacteriia': [[-0.6757266509268143, 0.32541280433816866],
  [-1.2323100285746866, 0.3779247860512955],
  [-1.2159771594672946, 0.5467739993752234],
  [-1.1196832583847025, 0.6907830678127567],
  [-0.6757266509268143, 0.32541280433816866]],
 'NA': [[-0.6757266509268144, -0.3254128043381685],
  [-0.8525286609044485, -0.8460530586000379],
  [-1.025669993230176, -0.8222257794359518],
  [-1.1907433498686864, -0.6873064574843211],
  [-1.151260595837927, -0.5248075363596111],
  [-1.2823374024299041, -0.4470844466167486],
  [-1.2172150700835507, -0.2917447945318761],
  [-0.6757266509268144, -0.3254128043381685]]}
ax=plt.gca()
ax.set_ylim([0,0])
ax.set_xlim([0,0])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticks([])
ax.set_yticklabels([])
Path=mpath.Path    
for tax in draw_order_coords:
     path_data=[]
     start=(Path.MOVETO, cloud_dict_formatted[tax][tax])
     path_data.append(start)
     curvecount=0
     for coord in draw_order_coords[tax][1:-1]:
         if curvecount==0:
              path_data.append((Path.LINETO, coord))
              print('l0',curvecount)
              curvecount+=1
         if curvecount>0 and not curvecount==len(draw_order_coords[tax][1:-1])+1:
              path_data.append((Path.CURVE3, coord))
              print('c',curvecount)
              curvecount+=1
         if not curvecount==0 and curvecount==len(draw_order_coords[tax][1:-1])+1:
              path_data.append((Path.LINETO, coord))
              print('lend',curvecount)
              curvecount+=1
     path_data.append((Path.MOVETO, cloud_dict_formatted[tax][tax]))
     codes, verts = zip(*path_data)
     path = mpath.Path(verts, codes)
     col=taxa_color_dict[tax]
     patch = patches.PathPatch(path, facecolor=[tx/255 for tx in taxa_color_dict[tax]], alpha=0.2, zorder=0)
     ax.add_patch(patch) 
     print(path_data)
     #x, y = zip(*path.vertices)
     #line, = ax.plot(x, y, 'go-')



####plot legend#####
####plot legend#####
#x=[-1.25,-.55,0.15, .85,-1.05,-0.2,0.65]
#y=[-1,-1,-1,-1,-1.25,-1.25,-1.25]
x=[0.95]*7
y=[1.24-(0.091*P) for P in list(range(0,7))]

s=[250]*7
plt.scatter(x, y, s, c='w', edgecolor='black', marker='H', linewidths=1.5)
s2=np.array([(sx-125)*0.85 for sx in s])
plt.scatter(x, y, s2, c='g', edgecolor='white', marker='H', linewidths=1)
s3=np.array([(sx-5)*0.95 for sx in s2])
#
labels=filtered_taxa_list
cols=[(taxa_color_dict[x][0]/255,taxa_color_dict[x][1]/255,taxa_color_dict[x][2]/255) for x in labels]
plt.scatter(x, y, s3, c=cols, edgecolor='black', marker='H', linewidths=1)
xtextpad=0.065
#get coords for text
font_dict= {'family': 'arial','color':  'black','weight': 'bold','size': 12}
for name, xc, yc in zip(labels,x, y):
    plt.text(xc+xtextpad, yc, name, fontdict=font_dict, verticalalignment='center')




ax.axis('equal')
# [left, bottom, width, height]
 ####HEATMAP#########################
ax = fig.add_axes([0.1, 0.675, 0.2, 0.6])
im=ax.imshow([actual_edge_weights], cmap=cmap, norm=norm_im, interpolation='none')
cb=plt.colorbar(im, orientation='horizontal', aspect=5, ticks=[0.001, 0.06, 0.125])#, cax=ax)
cb.ax.set_xticklabels([0, 0.06,0.125], weight='bold')
cb.ax.tick_params(labelsize=12) 
cb.outline.set_edgecolor('black')
cb.outline.set_linewidth(2)
ax.set_ylim([0,0])
ax.set_xlim([0,0])
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_yticks([])
ax.set_yticklabels([])
plt.gca().set_visible(False)
plt.axis('off')
plt.tight_layout()
#fig.savefig("NETWORKx_CURRENT_8_by_8_1200dpi.png",dpi=1200, bbox_inches='tight', pad_inches = 0)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')


##test
fig=plt.figure(figsize=(8,8))
ax=fig.add_subplot(111)
x=np.array(np.random.rand(10))
y=np.array(np.random.rand(10))
s=np.array(np.random.rand(10)*5000)
s2=np.array([(sx-125)*0.85 for sx in s])
ax.scatter(x, y, s, c='g', edgecolor='black', marker='H', linewidths=2)
ax.scatter(x, y, s2, c='g', edgecolor='white', marker='H', linewidths=2)
s3=np.array([(sx-50)*0.9 for sx in s2])
ax.scatter(x, y, s3, c='g', edgecolor='black', marker='H', linewidths=1)








































######Find sequences for gene cynthesis

#
#
#df_synthesis=deepcopy(df)
#df_synthesis['sum']=df_synthesis[['W1_C201','W1_3_C429','W1_4_C429','W3_C201','W3_2_C429','W3_1_C429','W5_C201','W5_25_C429','W5_100_C429','W10_C201','W10_100_i_C429','W10_100_ii_C429']].sum(axis=1)
#df_synthesis['rank']=df_synthesis['sum'].rank(ascending=False)
#df_synthesis=df_synthesis.sort_values('rank', axis=0)
#
#contig_seqs={}
#loop=0
#
#for seq in SeqIO.parse('Raw_contig_files/Trinity_DL_gt500_t.fasta', 'fasta'):
#        contig_seqs[seq.id]=seq.seq
#        if loop%100000==0:
#            print('Iterated over ' + str(loop) + 'seqs')
#        loop+=1
#
#nuc=[]
#trans=[]
#contig=[]
#
#for x in df_synthesis.index:
#    orf=x.split(' ')[-1].strip('()')
#    nuc.append(seqs_nuc[orf].seq)
#    trans.append(seqs[orf].seq)
#    #split on the ORF delimiter
#    contignam=orf.rsplit('_',1)[0]
#    contig.append(str(contig_seqs[contignam]))
#df_synthesis['trans']=trans
#df_synthesis['nuc']=nuc
#df_synthesis['contig_nuc']=contig
#
#def get_stats(indf, orf):
#    for x in indf.index:
#        #if bool(re.search(orf, x))==True:
#        if x == orf:
#            rank=str(indf.loc[x]['rank'])
#            nuc=str(indf.loc[x]['nuc'])
#            trans=str(indf.loc[x]['trans'])
#            contig=indf.loc[x]['contig_nuc']
#            print('ORF: ' + str(orf)+'\nProtein: '+str(trans)+'\nNuc: '+str(nuc)+'\nContig: '+str(contig)+'\nRank: ' +str(rank))
#
#get_stats(df_synthesis, 'GH6 (a_c489891_g1_i1_1)')
#get_stats(df_synthesis, 'GH5, CBM2 (a_c566843_g4_i1_1)')
#get_stats(df_synthesis, 'GH6 (b_TRINITY_DN712366_c0_g2_i1_1)') 
#get_stats(df_synthesis, 'GH5 (c_c998769_g1_i1_2)')
#get_stats(df_synthesis, 'GH5 (b_TRINITY_DN883935_c0_g1_i1_1)')
#get_stats(df_synthesis, 'GH5 (c_c957016_g1_i1_1)')
#get_stats(df_synthesis, 'GH5 (b_TRINITY_DN525671_c0_g1_i1_1)')
#get_stats(df_synthesis, 'GH5 (c_c1012251_g5_i1_1)')
#get_stats(df_synthesis, 'GH5 (c_c991288_g2_i1_1)')
#get_stats(df_synthesis, 'GH5 (a_c561927_g2_i2_1)')
#get_stats(df_synthesis, 'GH6 (d_TRINITY_DN1318669_c0_g6_i1_2)')
#
#get_stats(df_synthesis, 'GH6, CBM2, CBM10 (d_TRINITY_DN1318669_c0_g5_i1_2)')
#get_stats(df_synthesis, 'GH6, CBM44 (a_c544148_g1_i1_2)')
#get_stats(df_synthesis, 'GH5, CBM4, CBM2 (a_c580607_g3_i1_2)')
#get_stats(df_synthesis, 'GH5 (a_c577269_g3_i1_1)')
#get_stats(df_synthesis, 'GH5 (c_c987261_g2_i1_1)')
#
#get_stats(df_synthesis, 'GH3 (b_TRINITY_DN1479168_c0_g2_i2_3)')
#get_stats(df_synthesis, 'GH3 (b_TRINITY_DN183641_c0_g2_i1_1)')
#get_stats(df_synthesis, 'GH3 (c_c892952_g1_i2_2)')
#get_stats(df_synthesis, 'GH3 (b_TRINITY_DN1480675_c4_g1_i1_1)')
#get_stats(df_synthesis, 'GH3 (c_c219158_g1_i1_1)')
#get_stats(df_synthesis, 'GH3 (d_TRINITY_DN1285524_c0_g2_i2_1)')
#get_stats(df_synthesis, 'GH3 (b_TRINITY_DN915448_c0_g1_i1_1)')
#get_stats(df_synthesis, 'GH3 (c_c2259399_g1_i1_2)')
#get_stats(df_synthesis, 'GH3 (b_TRINITY_DN1462741_c0_g1_i3_1)')
#get_stats(df_synthesis, 'GH3 (c_c988814_g6_i1_1)')
#get_stats(df_synthesis, 'GH3 (c_c977743_g2_i1_1)')
#get_stats(df_synthesis, 'GH3 (a_c51752_g1_i1_1)')
#get_stats(df_synthesis, 'GH3 (b_TRINITY_DN1426143_c0_g2_i3_1)')
#
#GH3 (c_c988814_g6_i1_1)
#GH3 (c_c977743_g2_i1_1)
#GH3 (a_c51752_g1_i1_1)
#GH3 (a_c554422_g1_i1_1)
#GH3 (c_c1007998_g3_i1_1)
#GH5, GH30, CBM35 (c_c470125_g1_i1_1)
#GH3 (c_c977534_g2_i1_1)
#GH30 (a_c525221_g1_i1_1)
#GH3 (c_c977534_g3_i1_1)
#GH30 (b_TRINITY_DN1496319_c1_g1_i6_6)
#GH3 (c_c948380_g1_i1_1)
#GH3 (b_TRINITY_DN1426143_c0_g2_i3_1)
#GH3 (c_c977743_g1_i1_1)
#GH31 (a_c581716_g1_i1_2)
#sample_to_week_dict={'W1_C201':'One','W1_3_C429':'One','W1_4_C429':'One','W3_C201':'Three','W3_2_C429':'Three','W3_1_C429':'Three','W5_C201':'Five','W5_25_C429':'Five','W5_100_C429':'Five','W10_C201':'Ten','W10_100_i_C429':'Ten','W10_100_ii_C429':'Ten'}
#key=str(key)
#target_index=df.index.get_loc(key)
#temp_df=df.iloc[target_index][['W1_C201','W1_3_C429', 'W1_4_C429','W3_C201','W3_2_C429','W3_1_C429', 'W5_C201','W5_25_C429','W5_100_C429','W10_C201','W10_100_i_C429','W10_100_ii_C429']]
#temp_df=temp_df.reset_index()
#temp_df['Week']=temp_df['index'].map(sample_to_week_dict)
plt.close('all')






























#pull out CE1 sequences##
read_and_writece1='no'
CE1_list=[]
for x in annotation_dict:
    for y in annotation_dict[x]:
        if bool(re.search('CE1_',y))==True:
            CE1_list.append(x)
            print(x, y)
#            
#CE1_file_nuc='CE1_targets/CE1_ISDE_targets_nuc.fna'
#CE1_file_prot='CE1_targets/CE1_ISDE_targets_protein.fna'
#######WRITE_all_dbCANN_annotated_seqs_to_file###########
#if read_and_writece1 == 'yes':
#     ce1seqs=[]
#     ORFcount=0
#     for seq in SeqIO.parse(all_peptide_matching_ORFs_file, 'fasta'):
#         if seq.id in CE1_list:
#             ce1seqs.append(seq)
#             ORFcount+=1
#SeqIO.write(ce1seqs, CE1_file_prot, 'fasta')
#print("Written " + str(ORFcount) + " sequences in total out of a total of " + str(len(CE1_list)) + "\nWriting complete\n") 
#  
#######WRITE_all_dbCANN_annotated_seqs_to_file_NUCLEOTIDE###########
#if read_and_writece1 == 'yes':
#     ce1seqsprot=[]
#     ORFcount=0
#     for seq in SeqIO.parse(all_peptide_matching_ORFs_file_nuc, 'fasta'):
#         if seq.id in CE1_list:
#             ce1seqsprot.append(seq)
#             ORFcount+=1
#SeqIO.write(ce1seqsprot, CE1_file_nuc, 'fasta')
#print("Written " + str(ORFcount) + " sequences in total out of a total of " + str(len(CE1_list)) + "\nWriting complete\n")           

#######Extract_FUll_length_CONTIGS###########
contig_seqs={}
loop=0
for seq in SeqIO.parse('Raw_contig_files/Trinity_DL_gt500_t.fasta', 'fasta'):
        contig_seqs[seq.id]=seq
        if loop%100000==0:
            print('Iterated over ' + str(loop) + ' seqs')
        loop+=1

CE1_contigs={}
CE1_contigslist=[]
c=0
for orf in CE1_list:
    if orf[:-2] in contig_seqs:
        CE1_contigs[orf]=contig_seqs[orf[:-2]]
        c+=1
        print(c)

for contig in CE1_contigs:
    print(contig+'\n'+CE1_contigs[contig].seq)




targets=['TRINITY_DN258850_c0_g1_i1_2','TRINITY_DN219547_c0_g3_i1_2','TRINITY_DN247265_c0_g1_i5_1','TRINITY_DN254144_c0_g2_i4_1']
c=0
for x in contig_seqs:
    if bool(re.search(target[:-5],x))==True:
        print(x, contig_seqs[x].seq)
    if c%100000==0:
        print('Iterated over ' + str(c) + ' seqs')    
    c+=1

top_30_transcriptome_CE1s=['b_TRINITY_DN1505956_c4_g3_i1_1',\
'b_TRINITY_DN1499015_c1_g2_i1_1',\
'b_TRINITY_DN1505956_c4_g2_i3_3',\
'b_TRINITY_DN1505301_c3_g4_i8_2',\
'b_TRINITY_DN1474402_c1_g1_i2_1',\
'b_TRINITY_DN1495014_c1_g10_i1_1',\
'd_TRINITY_DN1322502_c0_g4_i1_2',\
'c_c1015929_g2_i1_3',\
'c_c998217_g1_i4_1',\
'd_TRINITY_DN1261140_c0_g4_i1_1',\
'c_c1013537_g2_i3_3',\
'b_TRINITY_DN1506518_c10_g10_i6_2',\
'c_c996741_g2_i1_3',\
'c_c999949_g1_i1_1',\
'a_c564951_g1_i1_1',\
'b_TRINITY_DN1484296_c1_g2_i1_1',\
'b_TRINITY_DN1484296_c1_g2_i5_2',\
'c_c998217_g1_i3_3',\
'c_c997460_g1_i1_1',\
'a_c489878_g1_i2_1',\
'd_TRINITY_DN1319199_c1_g4_i3_2',\
'a_c583258_g1_i3_3',\
'c_c1012355_g1_i1_2',\
'c_c998217_g1_i2_3',\
'b_TRINITY_DN913134_c0_g1_i1_2',\
'c_c1006875_g1_i1_1',\
'd_TRINITY_DN1308080_c5_g3_i1_1',\
'c_c1016573_g1_i1_2',\
'b_TRINITY_DN1502492_c3_g5_i4_1']

top_30_transcriptome_CE1s_2=[x[:-2] for x in top_30_transcriptome_CE1s]
    
target=['b_TRINITY_DN1505956_c4_g3_i1_1']    
targets=['TRINITY_DN258850_c0_g1_i1_2','TRINITY_DN219547_c0_g3_i1_2','TRINITY_DN247265_c0_g1_i5_1','TRINITY_DN254144_c0_g2_i4_1']
targets2=['TRINITY_DN258850_c0_g1_i1','TRINITY_DN219547_c0_g3_i1','TRINITY_DN247265_c0_g1_i5','TRINITY_DN254144_c0_g2_i4']
loop=0
#for seq in SeqIO.parse('100k_fasta_files/Dereplicated_ORFs_cat_nucleotide_gt300_from_gt50008_02__11_28_30.fna','fasta'):
for seq in SeqIO.parse('Raw_contig_files/Trinity_DL_gt500_t.fasta', 'fasta'):
    if seq.id in top_30_transcriptome_CE1s_2:
        print('>'+seq.id+'\n'+seq.seq)
    if loop%100000==0:
            print('Iterated over ' + str(loop) + ' seqs')
    loop+=1
    
from Bio.Seq import Seq   
from Bio.Alphabet import generic_dna    
contig=Seq(t, generic_dna)
print(start+contig.lower().split(start)[1].split(end)[0]+end)


for x in df_mean_CAZy_ORF.index:
    if bool(re.search('CE1,', x))==True or bool(re.search('CE1 ', x))==True:
        print(x, df_mean_CAZy_ORF.loc[x]['sum_rank'])
        
        
        
        
        
        
        
        
#######Percentage variability for paper ##############
        
        
tx=df3[['Taxonomy','Week','Molar_percentage']]
txgb=tx.groupby(['Taxonomy'])['Molar_percentage'].sum()
txgb/sum(txgb)*100

#groupby week to get sum of molar percentage
weeksums=tx.groupby(['Week'])['Molar_percentage'].sum().reset_index(drop=True)

a=tx[tx['Taxonomy']=='Sphingobacteriia']['Molar_percentage'].reset_index(drop=True)
weekvals=a/weeksums*100
weekvals.std()



#############Construct final table for 
#use df_master_mean and tax_master_dict
##Protein\tTophit\tTophit_eval\tAccession\tTaxid\tGenus\tspecies


#
f=open('blastP\BlastP_result_files\Taxonomy_from_blastP\Feature_table_taxonomy_no_filter_FULL_LINEAGE.tsv','w')
evallist=[]
linedict={}
for orf in annotation_dict:
    orfstr=str(orf)
    tophit=ORF_taxa_dict[orf]['Tophit']
    topeval=ORF_taxa_dict[orf]['Tophit_eval']
    accession=ORF_taxa_dict[orf]['Accession']
    taxid=taxa_master_dict[orf]['Taxid']
    phyl=taxa_master_dict[orf]['Taxonomy']['phylum']
    classx=taxa_master_dict[orf]['Taxonomy']['class']
    order=taxa_master_dict[orf]['Taxonomy']['order']
    fam=order=taxa_master_dict[orf]['Taxonomy']['family']
    genus=taxa_master_dict[orf]['Taxonomy']['genus']
    species=taxa_master_dict[orf]['Taxonomy']['species']
    evallist.append(float(topeval))
    line=orfstr+'\t'+tophit+'\t'+topeval+'\t'+accession+'\t'+taxid+'\t'+phyl+'\t'+classx+'\t'+order+'\t'+fam+'\t'+genus+'\t'+species+'\n'
    linedict[orf]=line
    f.write(line)
    print(orfstr+'\t'+tophit+'\t'+topeval+'\t'+accession+'\t'+taxid+'\t'+phyl+'\t'+classx+'\t'+order+'\t'+fam+'\t'+genus+'\t'+species)
f.close()
  
#plt.violinplot(evallist,showextrema=False)

  
orflist=[x for x in annotation_dict]
molpct_df=df_master.loc[orflist]

#molpct_df.to_csv('blastP\BlastP_result_files\Taxonomy_from_blastP\Feature_table_molpct_no_filter.csv')