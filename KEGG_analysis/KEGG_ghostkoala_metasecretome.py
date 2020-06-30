# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 11:17:18 2019

@author: dl923
"""
import xlrd
import pandas as pd
import numpy as np
import seaborn as sns
import random
import matplotlib.pyplot as plt
import copy
import pickle
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from Bio.KEGG import REST
#####Ghostkoala total proteome

#### FILTERING PARAMETERS###
e_value_filter=1e-0
############################


###files###
metasecretome_kegg_file='Proteomics/KEGG_GHOSTKoala/user_ko.txt'
file_master='Proteomics/mascot_100k_comp_searches_with_db_0pt05_dbCAN/DB0pt05_vs_0pt1SIG/ISDE_Master_mascot_COMPOSITE_100k_search_0pt05DB_with_0pt1sig_11268U.xlsx'
file_dbCAN='Proteomics/mascot_100k_comp_searches_with_db_0pt05_dbCAN/DB0pt05_vs_0pt1SIG/dbCAN_OUTPUT_DB0pt05_vs_0pt1SIG.xlsx'

############ METASECRETOME KEGG ORTHOLOGS ##############
f=open(metasecretome_kegg_file, 'r')
metasecretome_kegg_dict={}
linecount=0
unclassified_count=0
unclassified_list=[]
for line in f:
    if line.count('\t')>0:
        ORF=line.split('\t')[0].strip('\n')
        kegg=line.split('\t')[1].strip('\n')
        if ORF not in metasecretome_kegg_dict:
            metasecretome_kegg_dict[ORF]=kegg
            linecount+=1
    elif line.count('\t')==0:
        ORF=line.split('\t')[0].strip('\n')
        kegg='Unassigned'  
        if ORF not in metasecretome_kegg_dict:
            metasecretome_kegg_dict[ORF]=kegg
            linecount+=1
            unclassified_count+=1
            unclassified_list.append(ORF)
print(str(linecount) + ' ORFs and associated keggs written to dict\n'+ \
      str(len(metasecretome_kegg_dict.keys())) + ' ORFs with keggs in dict\n'+ \
      str(len(unclassified_list)) + ' ORFs had unclassified or no homologous kegg orthologs (' + str(len(unclassified_list)/len(metasecretome_kegg_dict.keys())*100) + '%)')
###########
#mol% accounted for by unclassified KO proteins
df_master_mean[['One','Three','Five','Ten']].loc[unclassified_list].sum(axis=0)
np.std(df_master_mean[['One','Three','Five','Ten']].loc[unclassified_list].sum(axis=0))
np.mean(df_master_mean[['One','Three','Five','Ten']].loc[unclassified_list].sum(axis=0))
############ RAW METASECRETOME DICT ##############
workbook_master=xlrd.open_workbook(file_master)
#Write create a dictionary of sheetnames (key) and indices (create two lists then merge into dictionary using dict(zip(a+b)))
sheetlist=workbook_master.sheet_names()
indexlist=list(range(0,workbook_master.nsheets))
xl_dict=dict(zip(sheetlist, indexlist))
#filter only the sheets containing contigs and empai scores using dict comprehension
filter_xl_dict= {k: v for k, v in xl_dict.items() if k.startswith('ISDE')}
#####Write peptide matching contigs to list then dictionary

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
    
master_nested_dict={}
print("Reading spreadsheet " + str(file_master) + "\n")
print("Opening and extracting protein matching ORFs from spreadsheet")
for k in filter_xl_dict:
    wb=workbook_master.sheet_by_index(filter_xl_dict[k])
    create_dict(wb)   

#####FILTER UNWANTED KEYS################
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
   
###BUILD DF ####
df_master=pd.DataFrame.from_dict(master_nested_dict, orient='index')
df_master_mean=pd.DataFrame(index=df_master.index)
df_master_mean['One']=(df_master['W1_C201'] + df_master['W1_3_C429'] + df_master['W1_4_C429'])/3
df_master_mean['Three']=(df_master['W3_C201'] + df_master['W3_2_C429'] + df_master['W3_1_C429'])/3
df_master_mean['Five']=(df_master['W5_C201'] + df_master['W5_25_C429'] + df_master['W5_100_C429'])/3
df_master_mean['Ten']=(df_master['W10_C201'] + df_master['W10_100_i_C429'] + df_master['W10_100_ii_C429'])/3    


#### MAP ORFS TO KEGG ORTHOLOGS FROM GHOSTKOLA ####
df_master_mean['KO']=df_master_mean.index.map(metasecretome_kegg_dict)



######################GENERATE LINK DICTS#########################
#Get LINKS between pathway and KO numbers E.G...
#'ko:K00001\tpath:map00010\n',
#'ko:K00001\tpath:ko00010\n',
strat_pathway_links=REST.kegg_link(source_db="ko", target_db="pathway")
strat_pathway_link_dict={}
count=0
#format = #KO : link
for x in list(strat_pathway_links):
    key=x.split('\t')[0][3:]
    if x.split('\t')[1][5:].strip('\n').startswith('map'):
         link=x.split('\t')[1][5:].strip('\n')
         if key not in strat_pathway_link_dict:
             strat_pathway_link_dict[key]=[]
         if link not in strat_pathway_link_dict[key]:
             strat_pathway_link_dict[key].append(link)
         count+=1
    strat_pathway_link_dict['Unassigned']=['Unassigned']
print(str(count) + ' pathways iterated into dictionary')
##Unassigned means ORF has no KO
##Unclassified means KO has no pathway classifications

#################### QUALITY CHECK NUMBER OF KO THAT ARE MAPPABLE##########################

unmappable_KOn=[]
mappable_KOn=[]
for KOn in df_master_mean['KO']:
    if KOn not in strat_pathway_link_dict:
        if KOn not in unmappable_KOn:
            unmappable_KOn.append(KOn)
    if KOn in strat_pathway_link_dict:
        if KOn not in mappable_KOn:
            mappable_KOn.append(KOn)
            
print(str(len(unmappable_KOn)) + ' unmappable KO numbers')
print(str(len(mappable_KOn)) + ' mappable KO numbers')

#Append unmappable KO numbers to strat_pathway_link_dict as KO number : ['Unclassified'] 
#unclassified as string in list to be unpacked at next processing stage
for x in unmappable_KOn:
    if x not in strat_pathway_link_dict:
        strat_pathway_link_dict[x]=['Unclassified']



######################GENERATE PATHWAY DICTS######################### 
#get list of maps:pathway DEFINITION from pathway database E.G...
# 'path:map00010\tGlycolysis / Gluconeogenesis\n',
# 'path:map00020\tCitrate cycle (TCA cycle)\n',
# 'path:map00030\tPentose phosphate pathway\n'
strat_pathway_definition_list=REST.kegg_list("pathway")
strat_pathway_definition_dict={}
pcount=0  
#format = map number : definition
for y in list(strat_pathway_definition_list):
    mapstr=y.split('\t')[0][5:]
    definition=y.split('\t')[1].strip('\n')
    if key not in strat_pathway_definition_dict:
        strat_pathway_definition_dict[mapstr]=definition
    pcount+=1
strat_pathway_definition_dict['Unclassified']='Unclassified'
strat_pathway_definition_dict['Unassigned']='Unassigned'
#above - add unclassified : unclassified to deifnition dict
print(str(pcount) + ' definition appended to dictionary')
    


########## map KO to KO maps #############
df_master_mean['maps']=df_master_mean['KO'].map(strat_pathway_link_dict)
##drop index
df_master_mean=df_master_mean.reset_index()
df_master_mean_unmapped=pd.DataFrame({col:np.repeat(df_master_mean[col].values, df_master_mean['maps'].str.len())for col in df_master_mean.columns.difference(['maps'])}).assign(**{'maps':np.concatenate(df_master_mean['maps'].values)})[df_master_mean.columns.tolist()]
###map MAPS to pathways
df_master_mean_unmapped['pathway']=df_master_mean_unmapped['maps'].map(strat_pathway_definition_dict)


################## KEGG ATLAS ###########################
kegg_atlas={'Metabolism':{'Global and overview maps':['Metabolic pathways', 'Biosynthesis of secondary metabolites','Microbial metabolism in diverse environments', 'Biosynthesis of antibiotics', 'Carbon metabolism','2-Oxocarboxylic acid metabolism','Fatty acid metabolism','Biosynthesis of amino acids','Degradation of aromatic compounds'],
                          'Carbohydrate metabolism':['Glycolysis / Gluconeogenesis','Citrate cycle (TCA cycle)','Pentose phosphate pathway','Pentose and glucuronate interconversions','Fructose and mannose metabolism','Galactose metabolism','Ascorbate and aldarate metabolism','Starch and sucrose metabolism','Amino sugar and nucleotide sugar metabolism','Pyruvate metabolism','Glyoxylate and dicarboxylate metabolism','Propanoate metabolism','Butanoate metabolism','C5-Branched dibasic acid metabolism','Inositol phosphate metabolism'],
                          'Energy metabolism':['Oxidative phosphorylation','Photosynthesis','Photosynthesis - antenna proteins','Carbon fixation in photosynthetic organisms','Carbon fixation pathways in prokaryotes','Methane metabolism','Nitrogen metabolism','Sulfur metabolism'],
                          'Lipid metabolism':['Fatty acid biosynthesis','Fatty acid elongation','Fatty acid degradation','Synthesis and degradation of ketone bodies','Cutin, suberine and wax biosynthesis','Steroid biosynthesis','Primary bile acid biosynthesis','Secondary bile acid biosynthesis','Steroid hormone biosynthesis','Glycerolipid metabolism','Glycerophospholipid metabolism','Ether lipid metabolism','Sphingolipid metabolism','Arachidonic acid metabolism','Linoleic acid metabolism','alpha-Linolenic acid metabolism','Biosynthesis of unsaturated fatty acids'],
                          'Nucleotide metabolism':['Purine metabolism','Pyrimidine metabolism'],
                          'Amino acid metabolism':['Alanine, aspartate and glutamate metabolism', 'Glycine, serine and threonine metabolism', 'Cysteine and methionine metabolism', 'Valine, leucine and isoleucine degradation', 'Valine, leucine and isoleucine biosynthesis', 'Lysine biosynthesis', 'Lysine degradation', 'Arginine biosynthesis', 'Arginine and proline metabolism', 'Histidine metabolism', 'Tyrosine metabolism', 'Phenylalanine metabolism', 'Tryptophan metabolism', 'Phenylalanine, tyrosine and tryptophan biosynthesis'],
                          'Metabolism of other amino acids':['beta-Alanine metabolism', 'Taurine and hypotaurine metabolism', 'Phosphonate and phosphinate metabolism', 'Selenocompound metabolism', 'Cyanoamino acid metabolism', 'D-Glutamine and D-glutamate metabolism', 'D-Arginine and D-ornithine metabolism', 'D-Alanine metabolism', 'Glutathione metabolism'],
                          'Glycan biosynthesis and metabolism':['N-Glycan biosynthesis', 'Various types of N-glycan biosynthesis', 'Mucin type O-glycan biosynthesis', 'Mannose type O-glycan biosynthesis', 'Other types of O-glycan biosynthesis', 'Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate', 'Glycosaminoglycan biosynthesis - heparan sulfate / heparin', 'Glycosaminoglycan biosynthesis - keratan sulfate', 'Glycosaminoglycan degradation', 'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis', 'Glycosphingolipid biosynthesis - lacto and neolacto series', 'Glycosphingolipid biosynthesis - globo and isoglobo series', 'Glycosphingolipid biosynthesis - ganglio series', 'Lipopolysaccharide biosynthesis', 'Peptidoglycan biosynthesis', 'Other glycan degradation', 'Lipoarabinomannan (LAM) biosynthesis', 'Arabinogalactan biosynthesis - Mycobacterium'],
                          'Metabolism of cofactors and vitamins':['Thiamine metabolism', 'Riboflavin metabolism', 'Vitamin B6 metabolism', 'Nicotinate and nicotinamide metabolism', 'Pantothenate and CoA biosynthesis', 'Biotin metabolism', 'Lipoic acid metabolism', 'Folate biosynthesis', 'One carbon pool by folate', 'Retinol metabolism', 'Porphyrin and chlorophyll metabolism', 'Ubiquinone and other terpenoid-quinone biosynthesis'],
                          'Metabolism of terpenoids and polyketides':['Terpenoid backbone biosynthesis', 'Monoterpenoid biosynthesis', 'Sesquiterpenoid and triterpenoid biosynthesis', 'Diterpenoid biosynthesis', 'Carotenoid biosynthesis', 'Brassinosteroid biosynthesis', 'Insect hormone biosynthesis', 'Zeatin biosynthesis', 'Limonene and pinene degradation', 'Geraniol degradation', 'Type I polyketide structures', 'Biosynthesis of 12-, 14- and 16-membered macrolides', 'Biosynthesis of ansamycins', 'Biosynthesis of enediyne antibiotics', 'Biosynthesis of type II polyketide backbone', 'Biosynthesis of type II polyketide products', 'Tetracycline biosynthesis', 'Polyketide sugar unit biosynthesis', 'Nonribosomal peptide structures', 'Biosynthesis of siderophore group nonribosomal peptides', 'Biosynthesis of vancomycin group antibiotics'],
                          'Biosynthesis of other secondary metabolites':['Phenylpropanoid biosynthesis', 'Stilbenoid, diarylheptanoid and gingerol biosynthesis', 'Flavonoid biosynthesis', 'Flavone and flavonol biosynthesis', 'Anthocyanin biosynthesis', 'Isoflavonoid biosynthesis', 'Indole alkaloid biosynthesis', 'Indole diterpene alkaloid biosynthesis', 'Isoquinoline alkaloid biosynthesis', 'Tropane, piperidine and pyridine alkaloid biosynthesis', 'Acridone alkaloid biosynthesis', 'Caffeine metabolism', 'Betalain biosynthesis', 'Glucosinolate biosynthesis', 'Benzoxazinoid biosynthesis', 'Penicillin and cephalosporin biosynthesis', 'Carbapenem biosynthesis', 'Monobactam biosynthesis', 'Clavulanic acid biosynthesis', 'Streptomycin biosynthesis', 'Neomycin, kanamycin and gentamicin biosynthesis', 'Acarbose and validamycin biosynthesis', 'Puromycin biosynthesis', 'Novobiocin biosynthesis', 'Staurosporine biosynthesis', 'Phenazine biosynthesis', 'Prodigiosin biosynthesis', 'Aflatoxin biosynthesis', 'Biosynthesis of secondary metabolites - unclassified'],
                          'Xenobiotics biodegradation and metabolism':['Benzoate degradation', 'Aminobenzoate degradation', 'Fluorobenzoate degradation', 'Chloroalkane and chloroalkene degradation', 'Chlorocyclohexane and chlorobenzene degradation', 'Toluene degradation', 'Xylene degradation', 'Nitrotoluene degradation', 'Ethylbenzene degradation', 'Styrene degradation', 'Atrazine degradation', 'Caprolactam degradation', 'Bisphenol degradation', 'Dioxin degradation', 'Naphthalene degradation', 'Polycyclic aromatic hydrocarbon degradation', 'Furfural degradation', 'Steroid degradation', 'Metabolism of xenobiotics by cytochrome P450', 'Drug metabolism - cytochrome P450', 'Drug metabolism - other enzymes'],
                          'Chemical structure transformation maps':['Overview of biosynthetic pathways', 'Biosynthesis of plant secondary metabolites', 'Biosynthesis of phenylpropanoids', 'Biosynthesis of terpenoids and steroids', 'Biosynthesis of alkaloids derived from shikimate pathway', 'Biosynthesis of alkaloids derived from ornithine, lysine and nicotinic acid', 'Biosynthesis of alkaloids derived from histidine and purine', 'Biosynthesis of alkaloids derived from terpenoid and polyketide', 'Biosynthesis of plant hormones']},
            'Genetic Information Processing':{'Transcription':['RNA polymerase','Basal transcription factors','Spliceosome'],
                          'Translation':['Ribosome', 'Aminoacyl-tRNA biosynthesis', 'RNA transport', 'mRNA surveillance pathway', 'Ribosome biogenesis in eukaryotes'],
                          'Folding, sorting and degradation':['Protein export', 'Protein processing in endoplasmic reticulum', 'SNARE interactions in vesicular transport', 'Ubiquitin mediated proteolysis', 'Sulfur relay system', 'Proteasome', 'RNA degradation'],
                          'Replication and repair':['DNA replication', 'Base excision repair', 'Nucleotide excision repair', 'Mismatch repair', 'Homologous recombination', 'Non-homologous end-joining', 'Fanconi anemia pathway']},
           'Environmental Information Processing':{'Membrane transport':['ABC transporters','Phosphotransferase system (PTS)','Bacterial secretion system'],
                          'Signal transduction':['Two-component system', 'Ras signaling pathway', 'Rap1 signaling pathway', 'MAPK signaling pathway', 'MAPK signaling pathway - fly', 'MAPK signaling pathway - plant', 'MAPK signaling pathway - yeast', 'ErbB signaling pathway', 'Wnt signaling pathway', 'Notch signaling pathway', 'Hedgehog signaling pathway', 'Hedgehog signaling pathway - fly', 'TGF-beta signaling pathway', 'Hippo signaling pathway', 'Hippo signaling pathway - fly', 'Hippo signaling pathway - multiple species', 'VEGF signaling pathway', 'Apelin signaling pathway', 'JAK-STAT signaling pathway', 'NF-kappa B signaling pathway', 'TNF signaling pathway', 'HIF-1 signaling pathway', 'FoxO signaling pathway', 'Calcium signaling pathway', 'Phosphatidylinositol signaling system', 'Phospholipase D signaling pathway', 'Sphingolipid signaling pathway', 'cAMP signaling pathway', 'cGMP-PKG signaling pathway', 'PI3K-Akt signaling pathway', 'AMPK signaling pathway', 'mTOR signaling pathway', 'Plant hormone signal transduction'],
                          'Signaling molecules and interaction':['Neuroactive ligand-receptor interaction', 'Cytokine-cytokine receptor interaction', 'ECM-receptor interaction', 'Cell adhesion molecules (CAMs)']},
           'Cellular Processes':{'Transport and catabolism':['Endocytosis', 'Phagosome', 'Lysosome', 'Peroxisome', 'Autophagy - animal', 'Autophagy - yeast', 'Autophagy - other', 'Mitophagy - animal', 'Mitophagy - yeast'],
                          'Cell growth and death':['Cell cycle', 'Cell cycle - yeast', 'Cell cycle - Caulobacter', 'Meiosis - yeast', 'Oocyte meiosis', 'Apoptosis', 'Apoptosis - fly', 'Apoptosis - multiple species', 'Ferroptosis', 'Necroptosis', 'p53 signaling pathway', 'Cellular senescence'],
                          'Cellular community - eukaryotes':['Focal adhesion', 'Adherens junction', 'Tight junction', 'Gap junction', 'Signaling pathways regulating pluripotency of stem cells'],
                          'Cellular community - prokaryotes':['Quorum sensing', 'Biofilm formation - Vibrio cholerae', 'Biofilm formation - Pseudomonas aeruginosa', 'Biofilm formation - Escherichia coli'],
                          'Cell motility':['Bacterial chemotaxis', 'Flagellar assembly', 'Regulation of actin cytoskeleton']},
           'Organismal Systems':{'Immune system':['Hematopoietic cell lineage', 'Complement and coagulation cascades', 'Platelet activation', 'Toll-like receptor signaling pathway', 'Toll and Imd signaling pathway', 'NOD-like receptor signaling pathway', 'RIG-I-like receptor signaling pathway', 'Cytosolic DNA-sensing pathway', 'C-type lectin receptor signaling pathway', 'Natural killer cell mediated cytotoxicity', 'Antigen processing and presentation', 'T cell receptor signaling pathway', 'Th1 and Th2 cell differentiation', 'Th17 cell differentiation', 'IL-17 signaling pathway', 'B cell receptor signaling pathway', 'Fc epsilon RI signaling pathway', 'Fc gamma R-mediated phagocytosis', 'Leukocyte transendothelial migration', 'Intestinal immune network for IgA production', 'Chemokine signaling pathway'],
                          'Endocrine system':['Insulin secretion', 'Insulin signaling pathway', 'Glucagon signaling pathway', 'Regulation of lipolysis in adipocytes', 'Adipocytokine signaling pathway', 'PPAR signaling pathway', 'GnRH signaling pathway', 'Ovarian steroidogenesis', 'Estrogen signaling pathway', 'Progesterone-mediated oocyte maturation', 'Prolactin signaling pathway', 'Oxytocin signaling pathway', 'Relaxin signaling pathway', 'Thyroid hormone synthesis', 'Thyroid hormone signaling pathway', 'Parathyroid hormone synthesis, secretion and action', 'Melanogenesis', 'Renin secretion', 'Renin-angiotensin system', 'Aldosterone synthesis and secretion', 'Cortisol synthesis and secretion'],
                          'Circulatory system':['Cardiac muscle contraction', 'Adrenergic signaling in cardiomyocytes', 'Vascular smooth muscle contraction'],
                          'Digestive system':['Salivary secretion', 'Gastric acid secretion', 'Pancreatic secretion', 'Bile secretion', 'Carbohydrate digestion and absorption', 'Protein digestion and absorption', 'Fat digestion and absorption', 'Cholesterol metabolism', 'Vitamin digestion and absorption', 'Mineral absorption'],
                          'Excretory system':['Vasopressin-regulated water reabsorption', 'Aldosterone-regulated sodium reabsorption', 'Endocrine and other factor-regulated calcium reabsorption', 'Proximal tubule bicarbonate reclamation', 'Collecting duct acid secretion'],
                          'Nervous system':['Glutamatergic synapse', 'GABAergic synapse', 'Cholinergic synapse', 'Dopaminergic synapse', 'Serotonergic synapse', 'Long-term potentiation', 'Long-term depression', 'Retrograde endocannabinoid signaling', 'Synaptic vesicle cycle', 'Neurotrophin signaling pathway'],
                          'Sensory system':['Phototransduction', 'Phototransduction - fly', 'Olfactory transduction', 'Taste transduction', 'Inflammatory mediator regulation of TRP channels'],
                          'Development':['Dorso-ventral axis formation', 'Axon guidance', 'Osteoclast differentiation'],
                          'Aging':['Longevity regulating pathway', 'Longevity regulating pathway - worm', 'Longevity regulating pathway - multiple species'],
                          'Environmental adaptation':['Circadian rhythm', 'Circadian entrainment', 'Circadian rhythm - fly', 'Circadian rhythm - plant', 'Thermogenesis', 'Plant-pathogen interaction']},
           'Human Diseases':{'Cancers: Overview':['Pathways in cancer', 'Central carbon metabolism in cancer', 'Choline metabolism in cancer', 'Transcriptional misregulation in cancer', 'MicroRNAs in cancer', 'Proteoglycans in cancer', 'Chemical carcinogenesis', 'Viral carcinogenesis'],
                          'Cancers: Specific types':['Colorectal cancer', 'Pancreatic cancer', 'Hepatocellular carcinoma', 'Gastric cancer', 'Glioma', 'Thyroid cancer', 'Acute myeloid leukemia', 'Chronic myeloid leukemia', 'Basal cell carcinoma', 'Melanoma', 'Renal cell carcinoma', 'Bladder cancer', 'Prostate cancer', 'Endometrial cancer', 'Breast cancer', 'Small cell lung cancer', 'Non-small cell lung cancer'],
                          'Immune diseases':['Asthma', 'Systemic lupus erythematosus', 'Rheumatoid arthritis', 'Autoimmune thyroid disease', 'Inflammatory bowel disease (IBD)', 'Allograft rejection', 'Graft-versus-host disease', 'Primary immunodeficiency'],
                          'Neurodegenerative diseases':['Alzheimer disease', 'Parkinson disease', 'Amyotrophic lateral sclerosis (ALS)', 'Huntington disease', 'Prion diseases'],
                          'Substance dependence':['Cocaine addiction', 'Amphetamine addiction', 'Morphine addiction', 'Nicotine addiction', 'Alcoholism'],
                          'Cardiovascular diseases':['Fluid shear stress and atherosclerosis', 'Hypertrophic cardiomyopathy (HCM)', 'Arrhythmogenic right ventricular cardiomyopathy (ARVC)', 'Dilated cardiomyopathy (DCM)', 'Viral myocarditis'],
                          'Endocrine and metabolic diseases':['Type II diabetes mellitus', 'Type I diabetes mellitus', 'Maturity onset diabetes of the young', 'Non-alcoholic fatty liver disease (NAFLD)', 'Insulin resistance', 'AGE-RAGE signaling pathway in diabetic complications', 'Cushing syndrome'],
                          'Infectious diseases: Bacterial':['Vibrio cholerae infection', 'Epithelial cell signaling in Helicobacter pylori infection', 'Pathogenic Escherichia coli infection', 'Salmonella infection', 'Shigellosis', 'Pertussis', 'Legionellosis', 'Staphylococcus aureus infection', 'Tuberculosis', 'Bacterial invasion of epithelial cells'],
                          'Infectious diseases: Viral':['Human T-cell leukemia virus 1 infection', 'Human immunodeficiency virus 1 infection', 'Measles', 'Influenza A', 'Hepatitis B', 'Hepatitis C', 'Herpes simplex infection', 'Human cytomegalovirus infection', 'Kaposi sarcoma-associated herpesvirus infection', 'Epstein-Barr virus infection', 'Human papillomavirus infection'],
                          'Infectious diseases: Parasitic':['Amoebiasis', 'Malaria', 'Toxoplasmosis', 'Leishmaniasis', 'Chagas disease (American trypanosomiasis)', 'African trypanosomiasis'],
                          'Drug resistance: Antimicrobial':['beta-Lactam resistance', 'Vancomycin resistance', 'Cationic antimicrobial peptide (CAMP) resistance'],
                          'Drug resistance: Antineoplastic':['EGFR tyrosine kinase inhibitor resistance', 'Platinum drug resistance', 'Antifolate resistance', 'Endocrine resistance']},
           'Drug Development':{'Chronology: Antiinfectives':['Penicillins', 'Cephalosporins - parenteral agents', 'Cephalosporins - oral agents', 'Aminoglycosides', 'Tetracyclines', 'Macrolides and ketolides', 'Quinolones', 'Rifamycins', 'Antifungal agents', 'Antiviral agents', 'Anti-HIV agents'],
                          'Chronology: Antineoplastics':['Antineoplastics - alkylating agents', 'Antineoplastics - antimetabolic agents', 'Antineoplastics - agents from natural products', 'Antineoplastics - hormones', 'Antineoplastics - protein kinases inhibitors'],
                          'Chronology: Nervous system agents':['Hypnotics', 'Anxiolytics', 'Anticonvulsants', 'Local analgesics', 'Opioid analgesics', 'Antipsychotics', 'Antipsychotics - phenothiazines', 'Antipsychotics - butyrophenones', 'Antidepressants', 'Agents for Alzheimer-type dementia', 'Antiparkinsonian agents'],
                          'Chronology: Other drugs':['Sulfonamide derivatives - overview', 'Sulfonamide derivatives - sulfa drugs', 'Sulfonamide derivatives - diuretics', 'Sulfonamide derivatives - hypoglycemic agents', 'Antiarrhythmic drugs', 'Antiulcer drugs', 'Immunosuppressive agents', 'Osteoporosis drugs', 'Antimigraines', 'Antithrombosis agents', 'Antirheumatics - DMARDs and biological agents', 'Antidiabetics', 'Antidyslipidemic agents', 'Antiglaucoma agents'],
                          'Target-based classification: G protein-coupled receptors':['Cholinergic and anticholinergic drugs', 'alpha-Adrenergic receptor agonists/antagonists', 'beta-Adrenergic receptor agonists/antagonists', 'Dopamine receptor agonists/antagonists', 'Histamine H receptor antagonists', 'Histamine H/H receptor agonists/antagonists', 'Serotonin receptor agonists/antagonists', 'Eicosanoid receptor agonists/antagonists', 'Opioid receptor agonists/antagonists', 'Angiotensin receptor and endothelin receptor antagonists'],
                          'Target-based classification: Nuclear receptors':['Glucocorticoid and mineralocorticoid receptor agonists/antagonists', 'Progesterone, androgen and estrogen receptor agonists/antagonists', 'Retinoic acid receptor (RAR) and retinoid X receptor (RXR) agonists/antagonists', 'Peroxisome proliferator-activated receptor (PPAR) agonists'],
                          'Target-based classification: Ion channels':['Nicotinic cholinergic receptor antagonists', 'GABA-A receptor agonists/antagonists', 'Calcium channel blocking drugs', 'Sodium channel blocking drugs', 'Potassium channel blocking and opening drugs', 'N-Metyl-D-aspartic acid receptor antagonists'],
                          'Target-based classification: Transporters':['Ion transporter inhibitors', 'Neurotransmitter transporter inhibitors'],
                          'Target-based classification: Enzymes':['Catecholamine transferase inhibitors', 'Cyclooxygenase inhibitors', 'HMG-CoA reductase inhibitors', 'Renin-angiotensin system inhibitors', 'HIV protease inhibitors'],
                          'Structure-based classification':['Quinolines', 'Eicosanoids', 'Prostaglandins'],
                          'Skeleton-based classification':['Benzoic acid family', ',-Diphenyl substitution family', 'Naphthalene family', 'Benzodiazepine family']}}

kegg_atlas_mappable_level_1={} # {l2:l1}
kegg_atlas_mappable_level_2={} # {l3:l2}
for l1 in kegg_atlas:
    for l2 in kegg_atlas[l1]:
        if l2 not in kegg_atlas_mappable_level_1:
            kegg_atlas_mappable_level_1[l2]=l1
        for l3 in kegg_atlas[l1][l2]:
            if l3 not in kegg_atlas_mappable_level_2:
                kegg_atlas_mappable_level_2[l3]=l2
kegg_atlas_mappable_level_1['Unclassified']='Unclassified'
kegg_atlas_mappable_level_2['Unclassified']='Unclassified'
kegg_atlas_mappable_level_1['Unassigned']='Unassigned'
kegg_atlas_mappable_level_2['Unassigned']='Unassigned'

#####################Map to different levels #######################
df_master_mean_unmapped['L2']=df_master_mean_unmapped['pathway'].map(kegg_atlas_mappable_level_2)
df_master_mean_unmapped['L1']=df_master_mean_unmapped['L2'].map(kegg_atlas_mappable_level_1)

##########################################################ANALYSIS#####################################################################
##########################################################ANALYSIS#####################################################################
##########################################################ANALYSIS#####################################################################
##########################################################ANALYSIS#####################################################################


target_level='L2'


df_target=df_master_mean_unmapped[['One','Three', 'Five', 'Ten', target_level]].groupby(target_level).sum()

#add L1 targets in exception list
exception_list=['Global and overview maps', 'Cancers: Overview', 'Unassigned']
# drop global exceptions
for ex in exception_list:
        if ex in df_target.index:
            df_target=df_target.drop(index=ex)
        for l1 in kegg_atlas:
            if ex == l1:
                for name in kegg_atlas[l1]:
                    if name in df_target.index:
                         df_target=df_target.drop(index=name)
            for l2 in kegg_atlas[l1]:
                if ex == l2:
                    for name in kegg_atlas[l1][l2]:
                        if name in df_target.index:
                             df_target=df_target.drop(index=name)

#reset index
df_target=df_target.reset_index()

def color_jenga(df):
    newpal=sns.xkcd_palette(["clear blue", "coral", "light aqua", "pale yellow", "light blue", "pinkish red"])
    length=int(len(df.index)/2)+1
    #plus 1 so as to not get caught out when int round 0.45 and 0.49 to 0...
    lightpal=sns.hls_palette(length, l=.5, s=.8, h=.5)
    darkpal=sns.hls_palette(length, l=.75, s=.9, h=0)
    for x in list(range(0, length)):
        newpal.append(lightpal[x])
        newpal.append(darkpal[x])
    random.shuffle(newpal)
    return newpal

#Plot as STACK ##
def plot_my_df_STACK(indf):
    colorz=color_jenga(indf)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(5,4))
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ####Plot stack
    SPLOT=ax.stackplot(indf.columns, [indf.loc[x] for x in indf.index], colors=colorz[:len(indf.index)], labels=indf.index)
   # for x in list(range(0, len(indf.index))):
    #    ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor=None, color=colors[x], label=indf.iloc[x].name)
    #    yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_xlim([0, max(position)])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xticklabels(list(indf.columns))
    sns.despine(top=True, right=True)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Function abundance (%)", fontsize=24)
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=1,title="Function", fontsize=10)
    legend.get_title().set_fontsize(20)
    plt.show()
    return fig, legend

collapse_list_L2_to_L1=['Genetic Information Processing', 'Environmental Information Processing', 'Organismal Systems']#, 'Human Diseases',]
dont_collapse_list=['Membrane transport']
##collapse different levels to show different hiercharies on the same plot, e.g. merge multiple L2 into L1, or single L2 into multiple L3
def RE_collapse_df(indf, kegg_atlas, collapse_list_L2_to_L1, dont_collapse_list, target_level):
    merge_dict={}
    if target_level=='L2':
        for L1 in collapse_list_L2_to_L1:
            for L2 in kegg_atlas[L1]:
                if L2 not in dont_collapse_list:
                     merge_dict[L2]=L1
    ##make sure other values map to itself
    for x in indf[target_level]:
        if x not in merge_dict:
            merge_dict[x]=x
    return merge_dict
        
merge_dict=RE_collapse_df(df_target, kegg_atlas, collapse_list_L2_to_L1, dont_collapse_list, target_level)
###MANUALLY appended collapsed stuff
manual_append_dict={'Amino acid metabolism': 'Amino acid metabolism (+ other)', \
                    'Metabolism of other amino acids': 'Amino acid metabolism (+ other)',\
                    'Metabolism of cofactors and vitamins' : 'Cofactor, vitamin, terpenoid and polyketide metabolism',\
                    'Metabolism of terpenoids and polyketides': 'Cofactor, vitamin, terpenoid and polyketide metabolism',\
                    'Carbon fixation pathways in prokaryotes': 'Carbon fixation: photosynthetic and prokaryotic', \
                    'Carbon fixation in photosynthetic organisms':'Carbon fixation: photosynthetic and prokaryotic'}
merge_dict.update(manual_append_dict)

def filter_df(indf, filter_pct):
    dff=copy.deepcopy(indf)
    #filter by percentage
    taxastr='Functions < ' + str(filter_pct)+"%"
    otherlist=[]
    for x in indf.columns:
        dff[x]=indf[x][((indf[x]/indf[x].sum())*100)>=filter_pct]    
        otherlist.append(indf[x][~(((indf[x]/indf[x].sum())*100)>=filter_pct)].sum(axis=0))
    row=pd.DataFrame([otherlist], columns=list(dff.columns), index=[taxastr])
    dff=dff.append(row)
    #remove empty index by sum
    dff=dff[dff.sum(axis=1)>0]
    return dff

def filter_df_taxa(indf, filter_pct):
    dff=copy.deepcopy(indf)
    #filter by percentage
    taxastr='Taxa < ' + str(filter_pct)+"%"
    otherlist=[]
    for x in indf.columns:
        dff[x]=indf[x][((indf[x]/indf[x].sum())*100)>=filter_pct]    
        otherlist.append(indf[x][~(((indf[x]/indf[x].sum())*100)>=filter_pct)].sum(axis=0))
    row=pd.DataFrame([otherlist], columns=list(dff.columns), index=[taxastr])
    dff=dff.append(row)
    #remove empty index by sum
    dff=dff[dff.sum(axis=1)>0]
    return dff

def normalise(indf):
    for col in indf.columns:
        indf[col]=indf[col]/(indf[col].sum())*100 
    return indf

#collapse high level into lower level
pathways_to_collapse=['Energy metabolism']#, 'Carbohydrate metabolism']#, 'Membrane transport']
def collapse_pathways(indf, pathways_to_collapse, reference_df, target_level):
    for path in pathways_to_collapse:
        tempdf=reference_df[reference_df[target_level]==path]
        tempdf2=tempdf[['One','Three', 'Five', 'Ten', 'pathway']].groupby('pathway').sum()
        indf=indf.drop(path)
        indf=indf.append(tempdf2)
    return indf
        

#merged low level pathways into higher levels    
df_target[target_level]=df_target[target_level].map(merge_dict, na_action=None)

##groupby high level first
df_target_x=df_target[['One','Three', 'Five', 'Ten', target_level]].groupby(target_level).sum()
##then append collapsed, grouped pathways
df_target_x_col=collapse_pathways(df_target_x, pathways_to_collapse, df_master_mean_unmapped, target_level)
###RE MAP then.... RE groupby for collapsed pathways
for x in df_target_x_col.index:
    if x not in merge_dict:
        merge_dict[x]=x
df_target_x_col=df_target_x_col.reset_index()
df_target_x_col[target_level]=df_target_x_col['index'].map(merge_dict, na_action=None)
df_target_x_col=df_target_x_col[['One','Three', 'Five', 'Ten', target_level]].groupby(target_level).sum()
#normalise
df_target_x_col=normalise(df_target_x_col)
#filter
testdf=filter_df(df_target_x_col, 1)
f, l = plot_my_df_STACK(testdf)
plt.show()
plt.clf()
plt.cla()
plt.close(f)  
plt.close('all')
##family



##########################################################ANALYSIS TAXA#####################################################################
##########################################################ANALYSIS TAXA#####################################################################
##########################################################ANALYSIS TAXA#####################################################################
##########################################################ANALYSIS TAXA#####################################################################
#load in taxonomy
tax_level='class'
with open('Proteomics/blastP/BlastP_result_files/Taxonomy_from_blastP/ISDE_blastp_taxonomy_entez_ete3.pickle', 'rb') as fp:
    taxa_master_dict = pickle.load(fp) 
tax_map_dict={}
for x in taxa_master_dict:
    tax_map_dict[x]=taxa_master_dict[x]['Taxonomy'][tax_level]
#add unanotated ORFs as 'NA' 
for x in df_master.index:
    if x not in tax_map_dict:
        tax_map_dict[x]='NA'
print(str(len(df_master)) +' total ORFs - ' + str(len(tax_map_dict)) + ' mappable taxa (including NAs)')

tax_level2='family'
tax_map_dict_family={}
for x in taxa_master_dict:
    tax_map_dict_family[x]=taxa_master_dict[x]['Taxonomy'][tax_level2]
#add unanotated ORFs as 'NA' 
for x in df_master.index:
    if x not in tax_map_dict_family:
        tax_map_dict_family[x]='NA'
print(str(len(df_master)) +' total ORFs - ' + str(len(tax_map_dict_family)) + ' mappable taxa (including NAs)')
   


    
df_taxa_master=copy.deepcopy(df_master_mean_unmapped)
df_taxa_master['taxa']=df_master_mean_unmapped['index'].map(tax_map_dict)
df_taxa_master['taxa_family']=df_master_mean_unmapped['index'].map(tax_map_dict_family)  

#######

df_taxa_grouped=df_taxa_master[['One','Three', 'Five', 'Ten', target_level, 'taxa']].groupby([target_level, 'taxa']).sum().reset_index()
df_taxa_grouped=df_taxa_grouped.set_index(target_level)

df_taxa_grouped_family=df_taxa_master[['One','Three', 'Five', 'Ten', target_level, 'taxa_family']].groupby([target_level, 'taxa_family']).sum().reset_index()
df_taxa_grouped_family=df_taxa_grouped_family.set_index(target_level)
#add L1 targets in exception list
#exception_list IS DEFINED ABOVE
# drop global exceptions
for ex in exception_list:
        if ex in df_taxa_grouped.index:
            df_taxa_grouped=df_taxa_grouped.drop(index=ex)
        for l1 in kegg_atlas:
            if ex == l1:
                for name in kegg_atlas[l1]:
                    if name in df_taxa_grouped.index:
                         df_taxa_grouped=df_taxa_grouped.drop(index=name)
            for l2 in kegg_atlas[l1]:
                if ex == l2:
                    for name in kegg_atlas[l1][l2]:
                        if name in df_taxa_grouped.index:
                             df_taxa_grouped=df_taxa_grouped.drop(index=name)
                             
###make new function to handle taxa for collapsable pathways
def collapse_pathways_taxa(indf, pathways_to_collapse, reference_df, target_level, taxa_level):
    for path in pathways_to_collapse:
        tempdf=reference_df[reference_df[target_level]==path]
        tempdf2=tempdf[['One','Three', 'Five', 'Ten', 'pathway', taxa_level]].groupby(['pathway', taxa_level]).sum().reset_index().set_index('pathway')
        indf=indf.drop(path)
        indf=indf.append(tempdf2)
    return indf

#Plot as STACK ##
def plot_my_df_STACK_non_norm(indf):
    colorz=color_jenga(indf)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(5,4))
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ####Plot stack
    SPLOT=ax.stackplot(indf.columns, [indf.loc[x] for x in indf.index], colors=colorz[:len(indf.index)], labels=indf.index)
   # for x in list(range(0, len(indf.index))):
    #    ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor=None, color=colors[x], label=indf.iloc[x].name)
    #    yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_xlim([0, max(position)])
    #ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xticklabels(list(indf.columns))
    sns.despine(top=True, right=True)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel(''u"\u2211"r' $\bar{x}$'" mol%", fontsize=24)
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=1,title="Function", fontsize=10)
    legend.get_title().set_fontsize(20)
    plt.show()
    return fig, legend


#reset index
df_taxa_grouped=df_taxa_grouped.reset_index()
#merged low level pathways into higher levels    
df_taxa_grouped[target_level]=df_taxa_grouped[target_level].map(merge_dict, na_action=None)
##groupby high level first
df_taxa_grouped_x=df_taxa_grouped[['One','Three', 'Five', 'Ten', target_level, 'taxa']].groupby([target_level, 'taxa']).sum().reset_index()
#set index for collapse function
df_taxa_grouped_x=df_taxa_grouped_x.set_index('L2')
##then append collapsed, grouped pathways
df_taxa_grouped_x_col=collapse_pathways_taxa(df_taxa_grouped_x, pathways_to_collapse, df_taxa_master, target_level, 'taxa')
###RE MAP then.... RE groupby for collapsed pathways
df_taxa_grouped_x_col=df_taxa_grouped_x_col.reset_index()
df_taxa_grouped_x_col['index']=df_taxa_grouped_x_col['index'].map(merge_dict, na_action=None)
df_taxa_grouped_x_col_master=df_taxa_grouped_x_col[['One','Three', 'Five', 'Ten', 'taxa', 'index']].groupby(['index', 'taxa']).sum().reset_index()

###family level df####
df_taxa_grouped_family=df_taxa_grouped_family.reset_index()
#merged low level pathways into higher levels    
df_taxa_grouped_family[target_level]=df_taxa_grouped_family[target_level].map(merge_dict, na_action=None)
##groupby high level first
df_taxa_grouped_x_family=df_taxa_grouped_family[['One','Three', 'Five', 'Ten', target_level, 'taxa_family']].groupby([target_level, 'taxa_family']).sum().reset_index()
#set index for collapse function
df_taxa_grouped_x_family=df_taxa_grouped_x_family.set_index('L2')
##then append collapsed, grouped pathways
df_taxa_grouped_x_col_family=collapse_pathways_taxa(df_taxa_grouped_x_family, pathways_to_collapse, df_taxa_master, target_level, 'taxa_family')
###RE MAP then.... RE groupby for collapsed pathways
df_taxa_grouped_x_col_family=df_taxa_grouped_x_col_family.reset_index()
df_taxa_grouped_x_col_family['index']=df_taxa_grouped_x_col_family['index'].map(merge_dict, na_action=None)
df_taxa_grouped_x_col_master_family=df_taxa_grouped_x_col_family[['One','Three', 'Five', 'Ten', 'taxa_family', 'index']].groupby(['index', 'taxa_family']).sum().reset_index()




####TEST this DF is the same as above ###########
testdf2_for_comparison=df_taxa_grouped_x_col_master[['One','Three', 'Five', 'Ten', 'index']].groupby(['index']).sum()
#normalise
testdf2_for_comparison_norm=normalise(testdf2_for_comparison)
#filter
testdf2_for_comparison_norm_filt=filter_df(testdf2_for_comparison_norm, 1)
f, l = plot_my_df_STACK(testdf2_for_comparison_norm_filt)
plt.show()
plt.clf()
plt.cla()
plt.close(f)  
plt.close('all')


focus_groups=['Gammaproteobacteria', 'Deltaproteobacteria', 'Bacteroidia', 'Alphaproteobacteria', 'Flavobacteriia', 'Cytophagia']#, 'NA']

for x in focus_groups:
    tempdf=df_taxa_grouped_x_col_master[df_taxa_grouped_x_col_master['taxa']==x].set_index('index')[['One', 'Three', 'Five', 'Ten']]
    plot_my_df_STACK_non_norm(tempdf)


##############################GLOBAL PLOTS#######################################
##############################GLOBAL PLOTS#######################################
##############################GLOBAL PLOTS#######################################    

##global color dict
def color_jenga_from_list(df):
    newpal=sns.xkcd_palette(["light blue","grey", "light grey", "clear blue", "coral", "light aqua", "pale yellow", "pinkish red"])
    length=int(len(df)/2)+2
    #plus 1 so as to not get caught out when int round 0.45 and 0.49 to 0...
    lightpal=sns.hls_palette(length, l=.5, s=.8, h=.5)
    darkpal=sns.hls_palette(length, l=.75, s=.9, h=0)
    for x in list(range(0, length)):
        newpal.append(darkpal[length-1-x])
        newpal.append(lightpal[x])
    #random.shuffle(newpal)
    return newpal

master_color_dict={}
index_for_colors=set(list(df_taxa_grouped_x_col_master['index']))
colist=color_jenga_from_list(index_for_colors)
count=0
for x in index_for_colors:
    master_color_dict[x]=colist[count]
    count+=1

    
def plot_my_df_strat_stack_stacked(indf, coldict, groups, merge_dict, normalize, normalize_total, filterx_total):
    #split the figure into two gridspec boxes, top and bottom
    fig = plt.figure(figsize=(8,10.5))
    gs_top = GridSpec(4, 2, wspace=0.25, hspace=0.05)

    #(gridsize, row, col), (row pos, col pos)
    ax1 = fig.add_subplot(gs_top[0, 0])
    ax2 = fig.add_subplot(gs_top[0, 1])
    ax3 = fig.add_subplot(gs_top[1, 0])
    ax4 = fig.add_subplot(gs_top[1, 1])
    ax5 = fig.add_subplot(gs_top[2, 0])
    ax6 = fig.add_subplot(gs_top[2, 1])
    
    gs_bottom=GridSpec(4, 2, wspace=0, hspace=0.3)
    ax7 = fig.add_subplot(gs_bottom[3, :])
    
    ax=[ax1, ax2, ax3, ax4, ax5, ax6]
    #plt.subplot2grid((2,2),(1,0), colspan=cols, rowspan=1)
    sns.set(style="white")
    ### Plot the amount of data plotted
    ####Plot stack
    group_count=0
    #sns.despine(top=True, right=True)

    for a in ax:
         #Make sure every kegg is in the index or it will cause plotting errors, loop over and form new list
        # pathway_list=[]
         #for kegg in kegg_atlas[L1][L2]:
          #   if kegg in list(indf[indf['sequence']==groups[group_count]]['function']):
           #       pathway_list.append(kegg)
           
         focus_group_df=indf[indf['taxa']==groups[group_count]].set_index('index')[['One', 'Three', 'Five', 'Ten']]
         #dff=copy.deepcopy(focus_group_df)
         if normalize=='normalize':
              for x in focus_group_df.columns:
                  focus_group_df[x]=focus_group_df[x]/focus_group_df[x].sum()*100
    #reorder the dataframe by abundance (make it easier to plot)
         focus_group_df=focus_group_df.reindex(index=focus_group_df.sum(axis=1).rank(ascending=0).sort_values().index)
         colorz=[]
         patches=[]
         for x in focus_group_df.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                 # patches.append(mpatches.Patch(color=coldict[x], label=x))
         a.stackplot(focus_group_df.columns, [focus_group_df.loc[x] for x in focus_group_df.index], colors=colorz, edgecolor='black')   
         if normalize=='normalize':
              a.set_ylim([0,100])
              a.set_yticks([50, 100])
         #if not normalize == 'normalize':
             #a.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
         changedict={'Gammaproteobacteria':'γ-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria', 'Alphaproteobacteria':'α-proteobacteria'}
         if groups[group_count] in changedict:
                 newname=changedict[groups[group_count]]
         if groups[group_count] not in changedict:
                 newname=groups[group_count]
         a.set_ylabel(newname, fontsize=12)
         a.set_xlim([0,3])
         a.set_xticks([])
         group_count+=1 
    ax5.set_xticks(['One','Three','Five','Ten'])
    #ax5.set_xlabel('Week', fontsize=10)
    ax6.set_xticks(['One','Three','Five','Ten'])
    #ax6.set_xlabel('Week', fontsize=10)
    totaldf=indf[['One','Three', 'Five', 'Ten', 'index']].groupby(['index']).sum()
    if normalize_total =='normalize_total':
         for x in totaldf.columns:
              totaldf[x]=totaldf[x]/totaldf[x].sum()*100
         ax7.set_ylim([0,100])
         ax7.set_yticks([0, 50, 100])
    totaldf=filter_df(totaldf, filterx_total)
    coldict['Functions < ' + str(filterx_total)+"%"]=(0,0,0)
    totaldf=totaldf.reindex(index=totaldf.sum(axis=1).rank(ascending=0).sort_values().index)
    colorz=[]
    patches=[]
    label_change_dict={'Cofactor, vitamin, terpenoid and polyketide metabolism':'Cofactor, vitamin, terpenoid and\npolyketide metabolism',\
                       'Xenobiotics biodegradation and metabolism':'Xenobiotics biodegradation and\nmetabolism',\
                       'Carbon fixation: photosynthetic and prokaryotic':'Carbon fixation: photosynthetic and\nprokaryotic'}
    for x in totaldf.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                  if x in label_change_dict:
                      patches.append(mpatches.Patch(color=coldict[x], label=label_change_dict[x]))
                  elif x not in label_change_dict:
                    patches.append(mpatches.Patch(color=coldict[x], label=x))
    totaldf=totaldf.fillna(0)
    ax7.stackplot(totaldf.columns, [totaldf.loc[x] for x in totaldf.index], colors=colorz,edgecolor='black')     
    ax7.set_xlim([0,3])
    ax7.set_xticks([0,1,2,3])
    ax7.set_xlabel('Week', fontsize=14)
    #ax7.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
    ax7.yaxis.set_label_text("Normalised total secretome")
    #plt.subplots_adjust(hspace=0.05, wspace=0.0)
    legend=plt.legend(handles=patches[::-1],loc='upper left', bbox_to_anchor=(0,-.35), ncol=2,title="", fontsize=10, frameon=False)
    return fig
    ####
  
    
fig=plot_my_df_strat_stack_stacked(df_taxa_grouped_x_col_master, master_color_dict, focus_groups, merge_dict, 'normalize', 'normalize_total', 1)    
#fig.savefig("Proteomics/KEGG_GHOSTKoala/Normalised_metasecretome_KEGG.png", bbox_inches='tight', dpi=1200)

fig=plot_my_df_strat_stack_stacked(df_taxa_grouped_x_col_master, master_color_dict, focus_groups, merge_dict, 'no', 'no', 1)    
#fig.savefig("Proteomics/KEGG_GHOSTKoala/NON_Normalised_metasecretome_KEGG.png", bbox_inches='tight', dpi=1200)

fig=plot_my_df_strat_stack_stacked(df_taxa_grouped_x_col_master, master_color_dict, focus_groups, merge_dict, 'normalize', 'no', 1)  
#fig.savefig("Proteomics/KEGG_GHOSTKoala/Normalised_taxa_NON_norm_total_metasecretome_KEGG.png", bbox_inches='tight', dpi=1200)
  
fig=plot_my_df_strat_stack_stacked(df_taxa_grouped_x_col_master, master_color_dict, focus_groups, merge_dict, 'no', 'normalize_total', 1)  
#fig.savefig("Proteomics/KEGG_GHOSTKoala/NON_Norm_taxa_NORM_TOTAL_metasecretome_KEGG.png", bbox_inches='tight', dpi=1200)

####################################################################################################
focus_groups_large=['Deltaproteobacteria','Gammaproteobacteria', 'Bacteroidia', 'Alphaproteobacteria', 'Flavobacteriia', 'Cytophagia' ,'Epsilonproteobacteria', 'Clostridia']#, 'NA']

def plot_my_df_strat_stack_stacked_LARGE(indf, coldict, groups, merge_dict, normalize, normalize_total, filterx_total, otu_table, otu_filter):
    ylabelsize=13
    #split the figure into two gridspec boxes, top and bottom
    fig = plt.figure(figsize=(12,10))
    gs_top = GridSpec(4, 4, wspace=0.3, hspace=0.05)

    #(gridsize, row, col), (row pos, col pos)
    ax1 = fig.add_subplot(gs_top[0, 0])
    ax2 = fig.add_subplot(gs_top[0, 1])
    ax3 = fig.add_subplot(gs_top[0, 2])
    ax4 = fig.add_subplot(gs_top[0, 3])
    ax5 = fig.add_subplot(gs_top[1, 0])
    ax6 = fig.add_subplot(gs_top[1, 1])
    ax7 = fig.add_subplot(gs_top[1, 2])
    ax8 = fig.add_subplot(gs_top[1, 3])
    #ax9 = fig.add_subplot(gs_top[2, 1])
    
    gs_bottom=GridSpec(4, 4, wspace=0.325, hspace=0.325)
    ax_bot = fig.add_subplot(gs_bottom[2:, :2])
    ax_bar = fig.add_subplot(gs_bottom[2:, 2:])
 #   gs_bar=GridSpec(4, 2, wspace=0, hspace=0.2)
 #   ax_bar=fig.add_subplot(gs_bar[2:, 2:])
    
    ax=[ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    #plt.subplot2grid((2,2),(1,0), colspan=cols, rowspan=1)
    sns.set(style="white")
    ### Plot the amount of data plotted
    ####Plot stack
    group_count=0
    #sns.despine(top=True, right=True)

    for a in ax:
         #Make sure every kegg is in the index or it will cause plotting errors, loop over and form new list
        # pathway_list=[]
         #for kegg in kegg_atlas[L1][L2]:
          #   if kegg in list(indf[indf['sequence']==groups[group_count]]['function']):
           #       pathway_list.append(kegg)
           
         focus_group_df=indf[indf['taxa']==groups[group_count]].set_index('index')[['One', 'Three', 'Five', 'Ten']]
         #dff=copy.deepcopy(focus_group_df)
         if normalize=='normalize':
              for x in focus_group_df.columns:
                  focus_group_df[x]=focus_group_df[x]/focus_group_df[x].sum()*100
    #reorder the dataframe by abundance (make it easier to plot)
         focus_group_df=focus_group_df.reindex(index=focus_group_df.sum(axis=1).rank(ascending=0).sort_values().index)
         colorz=[]
         patches=[]
         for x in focus_group_df.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                 # patches.append(mpatches.Patch(color=coldict[x], label=x))
         a.stackplot(focus_group_df.columns, [focus_group_df.loc[x] for x in focus_group_df.index], colors=colorz,edgecolor='white')   
         if normalize=='normalize':
              a.set_ylim([0,100])
              a.set_yticks([20,40,60,80])
         #if not normalize == 'normalize':
             #a.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
         changedict={'Gammaproteobacteria':'γ-proteobacteria','Epsilonproteobacteria':'ε-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria', 'Alphaproteobacteria':'α-proteobacteria'}
         if groups[group_count] in changedict:
                 newname=changedict[groups[group_count]]
         if groups[group_count] not in changedict:
                 newname=groups[group_count]
         a.set_ylabel(newname, fontsize=ylabelsize)
         a.set_xlim([0,3])
         a.set_xticks([])
         group_count+=1 
    ax5.set_xticks(['One','Three','Five','Ten'])
    ax5.set_xticklabels(['1','3','5','10'])
    ax6.set_xticks(['One','Three','Five','Ten'])
    ax6.set_xticklabels(['1','3','5','10'])
    ax7.set_xticks(['One','Three','Five','Ten'])
    ax7.set_xticklabels(['1','3','5','10'])
    ax8.set_xticks(['One','Three','Five','Ten'])
    ax8.set_xticklabels(['1','3','5','10'])
    #change clostridia graph to integer ytick_labels
    ax8.set_yticks([0,1,2])
    totaldf=indf[['One','Three', 'Five', 'Ten', 'index']].groupby(['index']).sum()
    if normalize_total =='normalize_total':
         for x in totaldf.columns:
              totaldf[x]=totaldf[x]/totaldf[x].sum()*100
         ax_bot.set_ylim([0,100])
         ax_bot.set_yticks([0, 50, 100])
    totaldf=filter_df(totaldf, filterx_total)
    coldict['Functions < ' + str(filterx_total)+"%"]=(0,0,0)
    totaldf=totaldf.reindex(index=totaldf.sum(axis=1).rank(ascending=0).sort_values().index)
    colorz=[]
    patches=[]
    label_change_dict={'Cofactor, vitamin, terpenoid and polyketide metabolism':'Cofactor, vitamin, terpenoid and\npolyketide metabolism',\
                       'Xenobiotics biodegradation and metabolism':'Xenobiotics biodegradation\nand metabolism',\
                       'Carbon fixation: photosynthetic and prokaryotic':'Carbon fixation: photosynthetic\nand prokaryotic',\
                       'Biosynthesis of other secondary metabolites':'Biosynthesis of other secondary\nmetabolites',\
                       'Environmental Information Processing':'Environmental Information\nProcessing'}
    for x in totaldf.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                  if x in label_change_dict:
                      patches.append(mpatches.Patch(color=coldict[x], label=label_change_dict[x]))
                  elif x not in label_change_dict:
                    patches.append(mpatches.Patch(color=coldict[x], label=x))
    totaldf=totaldf.fillna(0)
    ax_bot.stackplot(totaldf.columns, [totaldf.loc[x] for x in totaldf.index], colors=colorz,edgecolor='white')     
    ax_bot.set_xlim([0,3])
    ax_bot.set_xticks([0,1,2,3])
    ax_bot.set_xlabel('Week', fontsize=14)
    ax_bot.set_yticks([20,40,60,80])
    ax_bot.set_xticks(['One','Three','Five','Ten'])
    ax_bot.set_xticklabels(['1','3','5','10'])
    #ax7.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
    ax_bot.set_ylabel('Total secretome', fontsize=ylabelsize)
    #plt.subplots_adjust(hspace=0.05, wspace=0.0)
    legend=ax_bot.legend(handles=patches[::-1],loc='upper left', bbox_to_anchor=(-0.1,-.15), ncol=2, fontsize=10, frameon=False, title='Function (Ascending abundance)')
    legend.get_title().set_fontsize(12)
    ####BARPLOT
    colorsbar=color_jenga(otu_table)
    #reorder the dataframe by abundance (make it easier to plot)
    otu_table=otu_table.reindex(index=otu_table.sum(axis=1).rank(ascending=0).sort_values().index)
    otu_table=otu_table.fillna(value=0)
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(otu_table.columns)))
    position=list(range(0, len(otu_table.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(otu_table.index))):
        labeln=otu_table.iloc[x].name
        if labeln in changedict:
            labeln=changedict[labeln]
        ax_bar.bar(position, otu_table.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor='black', color=colorsbar[x], label=labeln)
        yvalue_cumulative=np.array(yvalue_cumulative+otu_table.iloc[x])
    ax_bar.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax_bar.set_ylim([0,100.5])
    ax_bar.set_xticks(position)
    ax_bar.set_xticklabels(list(otu_table.columns))
    ax_bar.tick_params(direction='out', length=6, width=4)  
    ax_bar.tick_params(axis='both', labelsize=12)
    ax_bar.set_ylabel("Relative abundance (%)", fontsize=ylabelsize)
    ax_bar.set_xlabel('Week', fontsize=14)
    ax_bar.set_yticks([20,40,60,80])
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax_bar.get_legend_handles_labels()
    legendbar=ax_bar.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.05,-0.15), ncol=2,title="Taxonomy (≥"+str(filterx)+"%) (Ascending abundance)", fontsize=12, frameon=False)
    legendbar.get_title().set_fontsize(12)
    
    
    return fig, legend, legendbar
    ####

fig, l1, l2=plot_my_df_strat_stack_stacked_LARGE(df_taxa_grouped_x_col_master, master_color_dict, focus_groups_large, merge_dict, 'no', 'normalize_total', 1, df_filtered_normalised, 2)    
#AVOID callin bbox_inches='tight' as this rescales the figsize!!
#fig.savefig("Proteomics/KEGG_GHOSTKoala/Metasecretome_KEGG_Taxa_nonnorm_TOTAL_norm_BAR.png", bbox_inches='tight', dpi=1200)

#########################PAPER VERSION#############################################################
def plot_my_df_strat_stack_stacked_LARGE_paper(indf, coldict, groups, merge_dict, normalize, normalize_total, filterx_total, otu_table, otu_filter):
    ylabelsize=11
    #split the figure into two gridspec boxes, top and bottom
    fig = plt.figure(figsize=(8,7*1.5))
    gs_top = GridSpec(6, 4, wspace=0.3, hspace=0.05)

    #(gridsize, row, col), (row pos, col pos)
    ax1 = fig.add_subplot(gs_top[0, 0])
    ax2 = fig.add_subplot(gs_top[0, 1])
    ax3 = fig.add_subplot(gs_top[0, 2])
    ax4 = fig.add_subplot(gs_top[0, 3])
    ax5 = fig.add_subplot(gs_top[1, 0])
    ax6 = fig.add_subplot(gs_top[1, 1])
    ax7 = fig.add_subplot(gs_top[1, 2])
    ax8 = fig.add_subplot(gs_top[1, 3])
    #ax9 = fig.add_subplot(gs_top[2, 1])
    
    gs_bottom=GridSpec(6, 4, wspace=0.325, hspace=0.325)
    ax_bot = fig.add_subplot(gs_bottom[2:4, 0:2])
    ax_bar = fig.add_subplot(gs_bottom[2:4, 2:4])
    ax_leg1_bot = fig.add_subplot(gs_bottom[4:6, 0:2])
    ax_leg2_bar = fig.add_subplot(gs_bottom[4:6, 2:])
    
    ax=[ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    #plt.subplot2grid((2,2),(1,0), colspan=cols, rowspan=1)
    sns.set(style="white")
    ### Plot the amount of data plotted
    ####Plot stack
    group_count=0
    #sns.despine(top=True, right=True)

    for a in ax:
         focus_group_df=indf[indf['taxa']==groups[group_count]].set_index('index')[['One', 'Three', 'Five', 'Ten']]
         #dff=copy.deepcopy(focus_group_df)
         if normalize=='normalize':
              for x in focus_group_df.columns:
                  focus_group_df[x]=focus_group_df[x]/focus_group_df[x].sum()*100
    #reorder the dataframe by abundance (make it easier to plot)
         focus_group_df=focus_group_df.reindex(index=focus_group_df.sum(axis=1).rank(ascending=0).sort_values().index)
         colorz=[]
         patches=[]
         for x in focus_group_df.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                 # patches.append(mpatches.Patch(color=coldict[x], label=x))
         a.stackplot(focus_group_df.columns, [focus_group_df.loc[x] for x in focus_group_df.index], colors=colorz,edgecolor='white')   
         if normalize=='normalize':
              a.set_ylim([0,100])
              a.set_yticks([20,40,60,80])
         #if not normalize == 'normalize':
             #a.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
         changedict={'Gammaproteobacteria':'γ-proteobacteria','Epsilonproteobacteria':'ε-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria', 'Alphaproteobacteria':'α-proteobacteria'}
         if groups[group_count] in changedict:
                 newname=changedict[groups[group_count]]
         if groups[group_count] not in changedict:
                 newname=groups[group_count]
         a.set_ylabel(newname, fontsize=ylabelsize)
         a.set_xlim([0,3])
         a.set_xticks([])
         a.yaxis.labelpad=0
         a.tick_params(axis='both', which='major', labelsize=10, pad=-3)
         group_count+=1 
    #change clostridia graph to integer ytick_labels
    ax4.set_yticks([0,2,4,6,8])
    ax5.set_xticks(['One','Three','Five','Ten'])
    ax5.set_xticklabels(['1','3','5','10'])
    ax6.set_xticks(['One','Three','Five','Ten'])
    ax6.set_xticklabels(['1','3','5','10'])
    ax7.set_xticks(['One','Three','Five','Ten'])
    ax7.set_xticklabels(['1','3','5','10'])
    ax8.set_xticks(['One','Three','Five','Ten'])
    ax8.set_xticklabels(['1','3','5','10'])
    #change clostridia graph to integer ytick_labels
    ax8.set_yticks([0,1,2])
    totaldf=indf[['One','Three', 'Five', 'Ten', 'index']].groupby(['index']).sum()
    if normalize_total =='normalize_total':
         for x in totaldf.columns:
              totaldf[x]=totaldf[x]/totaldf[x].sum()*100
         ax_bot.set_ylim([0,100])
         ax_bot.set_yticks([0, 50, 100])
    totaldf=filter_df(totaldf, filterx_total)
    coldict['Functions < ' + str(filterx_total)+"%"]=(0,0,0)
    totaldf=totaldf.reindex(index=totaldf.sum(axis=1).rank(ascending=0).sort_values().index)
    colorz=[]
    patches=[]
    label_change_dict={'Cofactor, vitamin, terpenoid and polyketide metabolism':'Cofactor, vitamin, terpenoid and\npolyketide metabolism',\
                       'Xenobiotics biodegradation and metabolism':'Xenobiotics biodegradation\nand metabolism',\
                       'Carbon fixation: photosynthetic and prokaryotic':'Carbon fixation: photosynthetic\nand prokaryotic',\
                       'Biosynthesis of other secondary metabolites':'Biosynthesis of other secondary\nmetabolites',\
                       'Environmental Information Processing':'Environmental Information\nProcessing'}
    for x in totaldf.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                  if x in label_change_dict:
                      patches.append(mpatches.Patch(color=coldict[x], label=label_change_dict[x]))
                  elif x not in label_change_dict:
                    patches.append(mpatches.Patch(color=coldict[x], label=x))
    totaldf=totaldf.fillna(0)
    ax_bot.stackplot(totaldf.columns, [totaldf.loc[x] for x in totaldf.index], colors=colorz,edgecolor='white')     
    ax_bot.set_xlim([0,3])
    ax_bot.set_xticks([0,1,2,3])
    ax_bot.set_xlabel('Week', fontsize=12)
    ax_bot.set_yticks([20,40,60,80])
    ax_bot.set_xticks(['One','Three','Five','Ten'])
    ax_bot.set_xticklabels(['1','3','5','10'])
    ax_bot.tick_params(axis='both', which='major', labelsize=10, pad=2)
    ax_bot.tick_params(axis='both', which='major', labelsize=10, pad=-3)
    ax_bot.yaxis.labelpad=0
    #ax7.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
    ax_bot.set_ylabel('Total secretome', fontsize=ylabelsize)
    #plt.subplots_adjust(hspace=0.05, wspace=0.0)
    legend=ax_bot.legend(handles=patches[::-1],loc='upper left', bbox_to_anchor=(-0.1,-.15), ncol=2, fontsize=6, frameon=False, title='Function (Ascending abundance)')
    legend.get_title().set_fontsize(10)
    ####BARPLOT
    colorsbar=color_jenga(otu_table)
    #reorder the dataframe by abundance (make it easier to plot)
    otu_table=otu_table.reindex(index=otu_table.sum(axis=1).rank(ascending=0).sort_values().index)
    otu_table=otu_table.fillna(value=0)
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(otu_table.columns)))
    position=list(range(0, len(otu_table.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(otu_table.index))):
        labeln=otu_table.iloc[x].name
        if labeln in changedict:
            labeln=changedict[labeln]
        ax_bar.bar(position, otu_table.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor='black', color=colorsbar[x], label=labeln)
        yvalue_cumulative=np.array(yvalue_cumulative+otu_table.iloc[x])
    ax_bar.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax_bar.set_ylim([0,100.5])
    ax_bar.set_xticks(position)
    ax_bar.set_xticklabels(list(otu_table.columns))
    ax_bar.tick_params(direction='out', length=6, width=4)  
    ax_bar.tick_params(axis='both', labelsize=10, pad=2)
    ax_bar.set_ylabel("Relative abundance (%)", fontsize=ylabelsize)
    ax_bar.set_xlabel('Week', fontsize=12)
    ax_bar.set_yticks([20,40,60,80])
    ax_bar.yaxis.labelpad=0
    ax_bar.tick_params(axis='both', which='major', labelsize=10, pad=-3)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax_bar.get_legend_handles_labels()
    legendbar=ax_bar.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(-0.1,-0.135), ncol=2,title="Taxonomy (≥"+str(filterx)+"%) (Ascending abundance)", fontsize=9, frameon=False)
    legendbar.get_title().set_fontsize(10)
    #set alpha for subplots to 0 to make legend show through and remove axis
    ax_leg1_bot.patch.set_alpha(0)
    ax_leg2_bar.patch.set_alpha(0)
    ax_leg1_bot.axis('off')
    ax_leg2_bar.axis('off')
    return fig

fig=plot_my_df_strat_stack_stacked_LARGE_paper(df_taxa_grouped_x_col_master, master_color_dict, focus_groups_large, merge_dict, 'no', 'normalize_total', 1, df_filtered_normalised, 2)    
#fig.savefig("Proteomics/KEGG_GHOSTKoala/Metasecretome_KEGG_Taxa_nonnorm_TOTAL_norm_BAR_8_7_1200dpi_PAPER.png", dpi=1200)

family_level_global_tax_col_dict={'Flavobacteriaceae': (0.1411764705882353,
  0.47843137254901963,
  0.9921568627450981),
 'Desulfobulbaceae': (0.4588235294117647,
  0.7333333333333333,
  0.9921568627450981),
 'Marinilabiaceae': (0.8705882352941177,
  0.047058823529411764,
  0.3843137254901961),
 'Saccharospirillaceae': (1.0, 0.8196078431372549, 0.8745098039215686),
 'Alteromonadaceae': (1.0, 0.996078431372549, 0.47843137254901963),
 'Vibrionaceae': (0.5882352941176471, 0.9764705882352941, 0.4823529411764706),
 'SB-1': (0.984313725490196, 0.1607843137254902, 0.2627450980392157),
 'NA': (0.0, 0.0, 0.0),
 'Oceanospirillaceae': (0.7843137254901961,
  0.4588235294117647,
  0.7686274509803922),
 'Desulfuromonadaceae': (1.0, 1.0, 1.0),
 'Rhodobacteraceae': (0.6784313725490196,
  0.5058823529411764,
  0.3137254901960784),
 'OM60': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
 'Flammeovirgaceae': (0.8470588235294118,
  0.8627450980392157,
  0.8392156862745098),
 '[Marinicellaceae]': (0.3176470588235294,
  0.396078431372549,
  0.4470588235294118),
 'Halomonadaceae': (1.0, 0.3568627450980392, 0.0),
 '[Acidaminobacteraceae]': (0.03137254901960784, 1.0, 0.03137254901960784),
 'Peptococcaceae': (0.996078431372549,
  0.00392156862745098,
  0.6941176470588235),
 'Chromatiaceae': (0.01568627450980392, 0.8509803921568627, 1.0),
 'Pseudomonadaceae': (1.0, 0.0, 0.050980392156862744),
 'Desulfobacteraceae': (0.8, 0.6784313725490196, 0.3764705882352941),
 'Shewanellaceae': (0.9749999999999999,
  0.6361764705882356,
  0.5250000000000001),
 'Cryomorphaceae': (0.2899999999999995, 0.15000000000000002, 0.85),
 'Acholeplasmataceae': (0.9749999999999999,
  0.5250000000000001,
  0.5726470588235296),
 'Opitutaceae': (0.5370588235294114, 0.15000000000000002, 0.85),
 'Desulfovibrionaceae': (0.9749999999999999,
  0.5250000000000001,
  0.7314705882352943),
 'Taxa < 1%': (0.7841176470588231, 0.15000000000000002, 0.85),
 'Piscirickettsiaceae': (0.9749999999999999,
  0.5250000000000001,
  0.8902941176470589),
 'Bacillaceae': (0.85, 0.15000000000000002, 0.6688235294117643),
 'Oleiphilaceae': (0.9008823529411762, 0.5250000000000001, 0.9749999999999999),
 'Saprospiraceae': (0.85, 0.15000000000000002, 0.42176470588235243),
 'Cytophagaceae': (0.7420588235294117, 0.5250000000000001, 0.9749999999999999),
 'Colwelliaceae': (0.85, 0.15000000000000002, 0.17470588235294066),
 'Verrucomicrobiaceae': (0.583235294117647,
  0.5250000000000001,
  0.9749999999999999),
 'Pirellulaceae': (0.85, 0.3723529411764698, 0.15000000000000002),
 '211ds20': (0.5250000000000001, 0.6255882352941177, 0.9749999999999999),
 'Erythrobacteraceae': (0.85, 0.6194117647058826, 0.15000000000000002),
 'Campylobacteraceae': (0.5250000000000001,
  0.7844117647058824,
  0.9749999999999999),
 'Paenibacillaceae': (0.8335294117647066, 0.85, 0.15000000000000002),
 'Cohaesibacteraceae': (0.5250000000000001,
  0.943235294117647,
  0.9749999999999999),
 'Spirochaetaceae': (0.586470588235294, 0.85, 0.15000000000000002),
 'Helicobacteraceae': (0.5250000000000001,
  0.9749999999999999,
  0.8479411764705881),
 'Planococcaceae': (0.3394117647058831, 0.85, 0.15000000000000002),
 'Sphingobacteriaceae': (0.5250000000000001,
  0.9749999999999999,
  0.6891176470588237),
 'Bacteroidaceae': (0.15000000000000002, 0.85, 0.2076470588235299),
 'Chitinophagaceae': (0.5250000000000001,
  0.9749999999999999,
  0.5302941176470591),
 'Phyllobacteriaceae': (0.15000000000000002, 0.85, 0.4547058823529408),
 'Christensenellaceae': (0.6785294117647057,
  0.9749999999999999,
  0.5250000000000001),
 'Cyclobacteriaceae': (0.15000000000000002, 0.85, 0.7017647058823525),
 'Pseudoalteromonadaceae': (0.8373529411764703,
  0.9749999999999999,
  0.5250000000000001),
 'Ruminococcaceae': (0.15000000000000002, 0.7511764705882354, 0.85),
 'Lachnospiraceae': (0.9749999999999999,
  0.9538235294117646,
  0.5250000000000001),
 'Cellulomonadaceae': (0.15000000000000002, 0.5041176470588237, 0.85)}

###check taxa as a function of keggs
focus_metabs=['Sulfur metabolism', 'Xenobiotics biodegradation and metabolism','Nitrogen metabolism', 'Methane metabolism']
fams=[]
for x in focus_metabs:
    tempdf2=df_taxa_grouped_x_col_master_family[df_taxa_grouped_x_col_master_family['index']==x].set_index('taxa_family')[['One', 'Three', 'Five', 'Ten']]
    tempdf3=filter_df(tempdf2, 3)
    for x in tempdf3.index:
        if x not in fams:
            fams.append(x)
    plot_my_df_STACK_non_norm(tempdf3)
    
##get taxa levels
for fam in ind:
     for x in taxa_master_dict:
         if fam == taxa_master_dict[x]['Taxonomy']['family']:
             if taxa_master_dict[x]['Taxonomy']['class']=='Alphaproteobacteria':
                  print(fam + ' : ' + taxa_master_dict[x]['Taxonomy']['class'])
                  break    
#########    

missing=[]
for x in fams:
    if x not in family_level_global_tax_col_dict:
        missing.append(x)
missing_pal=sns.hls_palette(len(missing), l=.3, s=1, h=.5)
misspalcount=0
for x in missing:
    family_level_global_tax_col_dict[x]=missing_pal[misspalcount]
    misspalcount+=1
#########################PAPER VERSION WITH ENERGY METAB#############################################################
#########################PAPER VERSION WITH ENERGY METAB#############################################################
#########################PAPER VERSION WITH ENERGY METAB#############################################################
#########################PAPER VERSION WITH ENERGY METAB#############################################################
metab_list_in=['Sulfur metabolism', 'Methane metabolism', 'Xenobiotics biodegradation and metabolism']
def plot_my_df_strat_stack_stacked_LARGE_paper(indf, coldict, groups, merge_dict, normalize, normalize_total, filterx_total, otu_table, otu_filter, indf_fam, global_colo_dict, metab_list, metab_filter):
    ylabelsize=11
    #split the figure into two gridspec boxes, top and bottom
    fig = plt.figure(figsize=(8*2,7*1.5))
    gs_top = GridSpec(6, 8, wspace=0.3, hspace=0.05)
    
    #(gridsize, row, col), (row pos, col pos)
    ax1 = fig.add_subplot(gs_top[0, 0])
    ax2 = fig.add_subplot(gs_top[0, 1])
    ax3 = fig.add_subplot(gs_top[0, 2])
    ax4 = fig.add_subplot(gs_top[0, 3])
    ax5 = fig.add_subplot(gs_top[1, 0])
    ax6 = fig.add_subplot(gs_top[1, 1])
    ax7 = fig.add_subplot(gs_top[1, 2])
    ax8 = fig.add_subplot(gs_top[1, 3])
    #ax9 = fig.add_subplot(gs_top[2, 1])
    
    gs_bottom=GridSpec(6, 8, wspace=0.325, hspace=0.325)
    ax_bot = fig.add_subplot(gs_bottom[2:4, 0:2])
    ax_bar = fig.add_subplot(gs_bottom[2:4, 2:4])
    ax_leg1_bot = fig.add_subplot(gs_bottom[4:6, 0:2])
    ax_leg2_bar = fig.add_subplot(gs_bottom[4:6, 2:4])
    
    gs_right=GridSpec(6, 8, wspace=0.5, hspace=0.15)
    ax_tl = fig.add_subplot(gs_right[0:2, 4:6])
    ax_tr = fig.add_subplot(gs_right[0:2, 6:8])
    ax_bl = fig.add_subplot(gs_right[2:4, 4:6])
    ax_br = fig.add_subplot(gs_right[2:4, 6:8])
    
    metab_ax=[ax_tl, ax_tr, ax_br]#ax_bl
    
    ax=[ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    #plt.subplot2grid((2,2),(1,0), colspan=cols, rowspan=1)
    sns.set(style="white")
    ### Plot the amount of data plotted
    ####Plot stack
    group_count=0
    #sns.despine(top=True, right=True)
    
    for a in ax:
         focus_group_df=indf[indf['taxa']==groups[group_count]].set_index('index')[['One', 'Three', 'Five', 'Ten']]
         #dff=copy.deepcopy(focus_group_df)
         if normalize=='normalize':
              for x in focus_group_df.columns:
                  focus_group_df[x]=focus_group_df[x]/focus_group_df[x].sum()*100
    #reorder the dataframe by abundance (make it easier to plot)
         focus_group_df=focus_group_df.reindex(index=focus_group_df.sum(axis=1).rank(ascending=0).sort_values().index)
         colorz=[]
         patches=[]
         for x in focus_group_df.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                 # patches.append(mpatches.Patch(color=coldict[x], label=x))
         a.stackplot(focus_group_df.columns, [focus_group_df.loc[x] for x in focus_group_df.index], colors=colorz,edgecolor='white')   
         if normalize=='normalize':
              a.set_ylim([0,100])
              a.set_yticks([20,40,60,80])
         #if not normalize == 'normalize':
             #a.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
         changedict={'Gammaproteobacteria':'γ-proteobacteria','Epsilonproteobacteria':'ε-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria', 'Alphaproteobacteria':'α-proteobacteria'}
         if groups[group_count] in changedict:
                 newname=changedict[groups[group_count]]
         if groups[group_count] not in changedict:
                 newname=groups[group_count]
         a.set_ylabel(newname, fontsize=ylabelsize)
         a.set_xlim([0,3])
         a.set_xticks([])
         a.yaxis.labelpad=0
         a.tick_params(axis='both', which='major', labelsize=10, pad=-3)
         group_count+=1 
    #change clostridia graph to integer ytick_labels
    ax4.set_yticks([0,2,4,6,8])
    ax5.set_xticks(['One','Three','Five','Ten'])
    ax5.set_xticklabels(['1','3','5','10'])
    ax6.set_xticks(['One','Three','Five','Ten'])
    ax6.set_xticklabels(['1','3','5','10'])
    ax7.set_xticks(['One','Three','Five','Ten'])
    ax7.set_xticklabels(['1','3','5','10'])
    ax8.set_xticks(['One','Three','Five','Ten'])
    ax8.set_xticklabels(['1','3','5','10'])
    #change clostridia graph to integer ytick_labels
    ax8.set_yticks([0,1,2])
    totaldf=indf[['One','Three', 'Five', 'Ten', 'index']].groupby(['index']).sum()
    if normalize_total =='normalize_total':
         for x in totaldf.columns:
              totaldf[x]=totaldf[x]/totaldf[x].sum()*100
         ax_bot.set_ylim([0,100])
         ax_bot.set_yticks([0, 50, 100])
    totaldf=filter_df(totaldf, filterx_total)
    coldict['Functions < ' + str(filterx_total)+"%"]=(0,0,0)
    totaldf=totaldf.reindex(index=totaldf.sum(axis=1).rank(ascending=0).sort_values().index)
    colorz=[]
    patches=[]
    label_change_dict={'Cofactor, vitamin, terpenoid and polyketide metabolism':'Cofactor, vitamin, terpenoid and\npolyketide metabolism',\
                       'Xenobiotics biodegradation and metabolism':'Xenobiotics biodegradation\nand metabolism',\
                       'Carbon fixation: photosynthetic and prokaryotic':'Carbon fixation: photosynthetic\nand prokaryotic',\
                       'Biosynthesis of other secondary metabolites':'Biosynthesis of other secondary\nmetabolites',\
                       'Environmental Information Processing':'Environmental Information\nProcessing'}
    for x in totaldf.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                  if x in label_change_dict:
                      patches.append(mpatches.Patch(facecolor=coldict[x], label=label_change_dict[x], edgecolor='black'))
                  elif x not in label_change_dict:
                    patches.append(mpatches.Patch(facecolor=coldict[x], label=x, edgecolor='black'))
    totaldf=totaldf.fillna(0)
    ax_bot.stackplot(totaldf.columns, [totaldf.loc[x] for x in totaldf.index], colors=colorz,edgecolor='white')     
    ax_bot.set_xlim([0,3])
    ax_bot.set_xticks([0,1,2,3])
    ax_bot.set_xlabel('Week', fontsize=12)
    ax_bot.set_yticks([20,40,60,80])
    ax_bot.set_xticks(['One','Three','Five','Ten'])
    ax_bot.set_xticklabels(['1','3','5','10'])
    ax_bot.tick_params(axis='both', which='major', labelsize=10, pad=2)
    ax_bot.tick_params(axis='both', which='major', labelsize=10, pad=-3)
    ax_bot.yaxis.labelpad=0
    #ax7.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
    ax_bot.set_ylabel('Total secretome', fontsize=ylabelsize)
    #plt.subplots_adjust(hspace=0.05, wspace=0.0)
    legend=ax_bot.legend(handles=patches[::-1],loc='upper left', bbox_to_anchor=(-0.1,-.15), ncol=2, fontsize=6, frameon=False, title='Function (Ascending abundance)')
    legend.get_title().set_fontsize(10)
    ####BARPLOT
    colorsbar=color_jenga(otu_table)
    #reorder the dataframe by abundance (make it easier to plot)
    otu_table=otu_table.reindex(index=otu_table.sum(axis=1).rank(ascending=0).sort_values().index)
    otu_table=otu_table.fillna(value=0)
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(otu_table.columns)))
    position=list(range(0, len(otu_table.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(otu_table.index))):
        labeln=otu_table.iloc[x].name
        if labeln in changedict:
            labeln=changedict[labeln]
        ax_bar.bar(position, otu_table.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor='black', color=colorsbar[x], label=labeln)
        yvalue_cumulative=np.array(yvalue_cumulative+otu_table.iloc[x])
    ax_bar.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax_bar.set_ylim([0,100.5])
    ax_bar.set_xticks(position)
    ax_bar.set_xticklabels(list(otu_table.columns))
    ax_bar.tick_params(direction='out', length=6, width=4)  
    ax_bar.tick_params(axis='both', labelsize=10, pad=2)
    ax_bar.set_ylabel("Relative abundance (%)", fontsize=ylabelsize)
    ax_bar.set_xlabel('Week', fontsize=12)
    ax_bar.set_yticks([20,40,60,80])
    ax_bar.yaxis.labelpad=0
    ax_bar.tick_params(axis='both', which='major', labelsize=10, pad=-3)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax_bar.get_legend_handles_labels()
    legendbar=ax_bar.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(-0.05,-0.135), ncol=2,title="Taxonomy (≥"+str(filterx)+"%) (Ascending abundance)", fontsize=9, frameon=False)
    legendbar.get_title().set_fontsize(10)
    #set alpha for subplots to 0 to make legend show through and remove axis
    ax_leg1_bot.patch.set_alpha(0)
    ax_leg2_bar.patch.set_alpha(0)
    ax_leg1_bot.axis('off')
    ax_leg2_bar.axis('off')
    ###############################################last 4 plot#########################################
    ####add extra colors to global color palette
    fams=[]
    for x in metab_list:
         tempdf2=indf_fam[indf_fam['index']==x].set_index('taxa_family')[['One', 'Three', 'Five', 'Ten']]
         tempdf3=filter_df_taxa(tempdf2, metab_filter)    
         for x in tempdf3.index:
             if x not in fams:
                 fams.append(x)
    missing=[]
    for x in fams:
        if x not in family_level_global_tax_col_dict:
            missing.append(x)
    missing_pal=sns.hls_palette(len(missing), l=.3, s=1, h=.5)
    misspalcount=0
    for x in missing:
        global_colo_dict[x]=missing_pal[misspalcount]
        misspalcount+=1
    ###generate graphs
    metab_ax_count=0
    for axiss in metab_ax:
        tempdf2=indf_fam[indf_fam['index']==metab_list[metab_ax_count]].set_index('taxa_family')[['One', 'Three', 'Five', 'Ten']]
        tempdf3=filter_df(tempdf2, metab_filter) 
        tempdf3=tempdf3.fillna(value=0) 
        color_order=[]
        for x in tempdf3.index:
             if global_colo_dict[x] not in color_order:
                  color_order.append(global_colo_dict[x])
        axiss.stackplot(tempdf3.columns, [tempdf3.loc[x] for x in tempdf3.index], colors=color_order,edgecolor='white') 
        axiss.set_ylabel(str(metab_list[metab_ax_count]), fontsize=ylabelsize)
        if str(metab_list[metab_ax_count]) =='Xenobiotics biodegradation and metabolism':
            axiss.set_ylabel('Xenobiotics biodegradation\n and metabolism', fontsize=ylabelsize)
        axiss.set_xlim([0,3])
        axiss.set_xticks([0,1,2,3])
        #ax_bot.set_yticks([20,40,60,80])
        #ax_bot.set_xticks(['One','Three','Five','Ten'])
        axiss.set_xticklabels(['1','3','5','10'])
        axiss.tick_params(axis='both', which='major', labelsize=10, pad=2)
        axiss.tick_params(axis='both', which='major', labelsize=10, pad=-3)
        axiss.yaxis.labelpad=0
        ##
        metab_ax_count+=1
    ax_br.set_xlabel('Week', fontsize=12)
    ax_bl.set_xlabel('Week', fontsize=12)
    ax_tr.set_yticks([0, 1, 2])
    ax_br.set_yticks([0, 1, 2])
    ###final legend
    leg3_patches=[]
    for x in fams:
        if x.startswith('Clostridiales Family'):
            leg3_patches.append(mpatches.Patch(facecolor=global_colo_dict[x], label='Clostridiales Family XII.\nIncertae Sedis', edgecolor='black'))
        else:
             leg3_patches.append(mpatches.Patch(facecolor=global_colo_dict[x], label=x, edgecolor='black'))
    legend3=ax_bl.legend(handles=leg3_patches[::-1],loc='upper left', bbox_to_anchor=(0.05,-0.11), ncol=3,title="Taxonomy (≥"+str(metab_filter)+"%) (Ascending abundance)", fontsize=9, frameon=False, labelspacing=0.1)
    legend3.get_title().set_fontsize(10)
    return fig

filterx=2
df=tax_extract(merged_otu_table, tax_table, 'Class')
df_filtered_otu, df_filtered_normalised_class=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised_class, filterx)

fig=plot_my_df_strat_stack_stacked_LARGE_paper(df_taxa_grouped_x_col_master,\
                                               master_color_dict, focus_groups_large, merge_dict,\
                                               'no', 'normalize_total', 1, df_filtered_normalised_class,\
                                               2, df_taxa_grouped_x_col_master_family, family_level_global_tax_col_dict,\
                                               metab_list_in, 3)    
#fig.savefig("KEGG_Ortholog_graphs_paper_metabs_PAPER.png", dpi=1200)

###############PAPER VERSION WITH BACKGROUND#######################
metab_list_in=['Sulfur metabolism', 'Methane metabolism', 'Xenobiotics biodegradation and metabolism']
def plot_my_df_strat_stack_stacked_LARGE_paper(indf, coldict, groups, merge_dict, normalize, normalize_total, filterx_total, otu_table, otu_filter, indf_fam, global_colo_dict, metab_list, metab_filter):
    ylabelsize=10
    legendfont=9.25
    columnssspacing=0.25
    #split the figure into two gridspec boxes, top and bottom
    fig = plt.figure(figsize=(14,10))
    gs_top = GridSpec(6, 8, wspace=0.3, hspace=0.05)
    
    #(gridsize, row, col), (row pos, col pos)
    ax1 = fig.add_subplot(gs_top[0, 0])
    ax2 = fig.add_subplot(gs_top[0, 1])
    ax3 = fig.add_subplot(gs_top[0, 2])
    ax4 = fig.add_subplot(gs_top[0, 3])
    ax5 = fig.add_subplot(gs_top[1, 0])
    ax6 = fig.add_subplot(gs_top[1, 1])
    ax7 = fig.add_subplot(gs_top[1, 2])
    ax8 = fig.add_subplot(gs_top[1, 3])
    #ax9 = fig.add_subplot(gs_top[2, 1])
    
    gs_bottom=GridSpec(6, 8, wspace=0.325, hspace=0.325)
    ax_bot = fig.add_subplot(gs_bottom[2:4, 0:4])
    #ax_bar = fig.add_subplot(gs_bottom[2:4, 2:4])
    ax_leg1_bot = fig.add_subplot(gs_bottom[4:6, 0:2])
    ax_leg2_bar = fig.add_subplot(gs_bottom[4:6, 2:4])
    
    gs_right=GridSpec(6, 8, wspace=0.5, hspace=0.35)
    ax_tl = fig.add_subplot(gs_right[0:2, 4:6])
    ax_tr = fig.add_subplot(gs_right[0:2, 6:8])
    ax_bar = fig.add_subplot(gs_right[2:4, 4:6])
    ax_br = fig.add_subplot(gs_right[2:4, 6:8])
    
    metab_ax=[ax_tl, ax_tr, ax_br]#, ax_bar]#, ax_br]
    
    ax=[ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    #plt.subplot2grid((2,2),(1,0), colspan=cols, rowspan=1)
    sns.set(style="white")
    ### Plot the amount of data plotted
    ####Plot stack
    group_count=0
    #sns.despine(top=True, right=True)
    
    for a in ax:
         focus_group_df=indf[indf['taxa']==groups[group_count]].set_index('index')[['One', 'Three', 'Five', 'Ten']]
         #dff=copy.deepcopy(focus_group_df)
         if normalize=='normalize':
              for x in focus_group_df.columns:
                  focus_group_df[x]=focus_group_df[x]/focus_group_df[x].sum()*100
    #reorder the dataframe by abundance (make it easier to plot)
         focus_group_df=focus_group_df.reindex(index=focus_group_df.sum(axis=1).rank(ascending=0).sort_values().index)
         colorz=[]
         patches=[]
         for x in focus_group_df.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                 # patches.append(mpatches.Patch(color=coldict[x], label=x))
         a.stackplot(focus_group_df.columns, [focus_group_df.loc[x] for x in focus_group_df.index], colors=colorz,edgecolor='white')   
         if normalize=='normalize':
              a.set_ylim([0,100])
              a.set_yticks([20,40,60,80])
         #if not normalize == 'normalize':
             #a.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
         changedict={'Gammaproteobacteria':'γ-proteobacteria','Epsilonproteobacteria':'ε-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria', 'Alphaproteobacteria':'α-proteobacteria'}
         if groups[group_count] in changedict:
                 newname=changedict[groups[group_count]]
         if groups[group_count] not in changedict:
                 newname=groups[group_count]
         a.set_ylabel(newname, fontsize=ylabelsize)
         a.set_xlim([0,3])
         a.set_xticks([])
         a.yaxis.labelpad=0
         a.tick_params(axis='both', which='major', labelsize=10, pad=-3)
         group_count+=1 
    #change clostridia graph to integer ytick_labels
    ax4.set_yticks([0,2,4,6,8])
    ax5.set_xticks(['One','Three','Five','Ten'])
    ax5.set_xticklabels(['1','3','5','10'])
    ax6.set_xticks(['One','Three','Five','Ten'])
    ax6.set_xticklabels(['1','3','5','10'])
    ax7.set_xticks(['One','Three','Five','Ten'])
    ax7.set_xticklabels(['1','3','5','10'])
    ax8.set_xticks(['One','Three','Five','Ten'])
    ax8.set_xticklabels(['1','3','5','10'])
    #change clostridia graph to integer ytick_labels
    ax8.set_yticks([0,1,2])
    totaldf=indf[['One','Three', 'Five', 'Ten', 'index']].groupby(['index']).sum()
    if normalize_total =='normalize_total':
         for x in totaldf.columns:
              totaldf[x]=totaldf[x]/totaldf[x].sum()*100
         ax_bot.set_ylim([0,100])
         ax_bot.set_yticks([0, 50, 100])
    totaldf=filter_df(totaldf, filterx_total)
    coldict['Functions < ' + str(filterx_total)+"%"]=(0,0,0)
    totaldf=totaldf.reindex(index=totaldf.sum(axis=1).rank(ascending=0).sort_values().index)
    colorz=[]
    patches=[]
    label_change_dict={'Cofactor, vitamin, terpenoid and polyketide metabolism':'Cofactor, vitamin, terpenoid and\npolyketide metabolism',\
                   'Xenobiotics biodegradation and metabolism':'Xenobiotics biodegradation\nand metabolism',\
                   'Carbon fixation: photosynthetic and prokaryotic':'Carbon fixation: photosynthetic\nand prokaryotic',\
                   'Biosynthesis of other secondary metabolites':'Biosynthesis of other secondary\nmetabolites',\
                   'Environmental Information Processing':'Environmental Information\nProcessing'}
    for x in totaldf.index:
             if coldict[x] not in colorz:
                  colorz.append(coldict[x])
                  if x in label_change_dict:
                      patches.append(mpatches.Patch(facecolor=coldict[x], label=label_change_dict[x], edgecolor='black'))
                  elif x not in label_change_dict:
                    patches.append(mpatches.Patch(facecolor=coldict[x], label=x, edgecolor='black'))
    totaldf=totaldf.fillna(0)
    ax_bot.stackplot(totaldf.columns, [totaldf.loc[x] for x in totaldf.index], colors=colorz,edgecolor='white')     
    ax_bot.set_xlim([0,3])
    ax_bot.set_xticks([0,1,2,3])
    ax_bot.set_xlabel('Week', fontsize=10, labelpad=0)
    ax_bot.set_yticks([20,40,60,80])
    ax_bot.set_xticks(['One','Three','Five','Ten'])
    ax_bot.set_xticklabels(['1','3','5','10'])
    ax_bot.tick_params(axis='both', which='major', labelsize=8, pad=2)
    ax_bot.tick_params(axis='both', which='major', labelsize=8, pad=-3)
    ax_bot.yaxis.labelpad=0
    #ax7.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
    ax_bot.set_ylabel('Total secretome', fontsize=ylabelsize)
    #plt.subplots_adjust(hspace=0.05, wspace=0.0)
    legend=ax_bot.legend(handles=patches[::-1],loc='upper left', bbox_to_anchor=(0,-0.1), ncol=2, fontsize=legendfont, frameon=False, title='Function (Ascending abundance)', columnspacing=columnssspacing+5, labelspacing=0.1, handletextpad=0.5)
    legend.get_title().set_fontsize(10)
    ####BARPLOT
    colorsbar=color_jenga(otu_table)
    #reorder the dataframe by abundance (make it easier to plot)
    otu_table=otu_table.reindex(index=otu_table.sum(axis=1).rank(ascending=0).sort_values().index)
    otu_table=otu_table.fillna(value=0)
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(otu_table.columns)))
    position=list(range(0, len(otu_table.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(otu_table.index))):
        labeln=otu_table.iloc[x].name
        if labeln in changedict:
            labeln=changedict[labeln]
        ax_bar.bar(position, otu_table.iloc[x],width, bottom=yvalue_cumulative, linewidth=2, edgecolor='black', color=colorsbar[x], label=labeln)
        yvalue_cumulative=np.array(yvalue_cumulative+otu_table.iloc[x])
    ax_bar.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax_bar.set_ylim([0,100.5])
    ax_bar.set_xticks(position)
    ax_bar.set_xticklabels(list(otu_table.columns))
    ax_bar.tick_params(direction='out', length=6, width=4)  
    ax_bar.tick_params(axis='both', labelsize=8, pad=2)
    ax_bar.set_ylabel("Relative abundance (%)", fontsize=ylabelsize)
    ax_bar.set_xlabel('Week', fontsize=10, labelpad=0)
    ax_bar.set_yticks([20,40,60,80])
    ax_bar.yaxis.labelpad=0
    ax_bar.tick_params(axis='both', which='major', labelsize=10, pad=-3)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax_bar.get_legend_handles_labels()
    legendbar=ax_bar.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(-0.125,-0.1), ncol=2,title="Taxonomy (≥"+str(filterx)+"%) (Ascending abundance)", fontsize=legendfont, frameon=False, columnspacing=columnssspacing, handletextpad=0.5, handlelength=1, handleheight=1)
    legendbar.get_title().set_fontsize(9)
    #set alpha for subplots to 0 to make legend show through and remove axis
    ax_leg1_bot.patch.set_alpha(0)
    ax_leg2_bar.patch.set_alpha(0)
    ax_leg1_bot.axis('off')
    ax_leg2_bar.axis('off')
    ###############################################last 4 plot#########################################
    ####add extra colors to global color palette
    fams=[]
    for x in metab_list:
         tempdf2=indf_fam[indf_fam['index']==x].set_index('taxa_family')[['One', 'Three', 'Five', 'Ten']]
         tempdf3=filter_df_taxa(tempdf2, metab_filter)    
         for x in tempdf3.index:
             if x not in fams:
                 fams.append(x)
    missing=[]
    for x in fams:
        if x not in family_level_global_tax_col_dict:
            missing.append(x)
    missing_pal=sns.hls_palette(len(missing), l=.3, s=1, h=.5)
    misspalcount=0
    for x in missing:
        global_colo_dict[x]=missing_pal[misspalcount]
        misspalcount+=1
    ###generate graphs
    metab_ax_count=0
    for axiss in metab_ax:
        tempdf2=indf_fam[indf_fam['index']==metab_list[metab_ax_count]].set_index('taxa_family')[['One', 'Three', 'Five', 'Ten']]
        tempdf3=filter_df(tempdf2, metab_filter) 
        tempdf3=tempdf3.fillna(value=0) 
        color_order=[]
        for x in tempdf3.index:
             if global_colo_dict[x] not in color_order:
                  color_order.append(global_colo_dict[x])
        axiss.stackplot(tempdf3.columns, [tempdf3.loc[x] for x in tempdf3.index], colors=color_order,edgecolor='white') 
        axiss.set_ylabel(str(metab_list[metab_ax_count]), fontsize=ylabelsize)
        if str(metab_list[metab_ax_count]) =='Xenobiotics biodegradation and metabolism':
            axiss.set_ylabel('Xenobiotics biodegradation\n and metabolism', fontsize=ylabelsize)
        axiss.set_xlim([0,3])
        axiss.set_xticks([0,1,2,3])
        #ax_bot.set_yticks([20,40,60,80])
        #ax_bot.set_xticks(['One','Three','Five','Ten'])
        axiss.set_xticklabels(['1','3','5','10'])
        axiss.tick_params(axis='both', which='major', labelsize=10, pad=2)
        axiss.tick_params(axis='both', which='major', labelsize=10, pad=-3)
        axiss.yaxis.labelpad=0
        axiss.set_xlabel('Week', fontsize=10, labelpad=0)
        ##
        metab_ax_count+=1
    ax_br.set_xlabel('Week', fontsize=10)
    #ax_bl.set_xlabel('Week', fontsize=12)
    ax_tr.set_yticks([0, 1, 2])
    #ax_br.set_yticks([0, 1, 2])
    ###final legend
    leg3_patches=[]
    for x in fams:
        if x.startswith('Clostridiales Family'):
            leg3_patches.append(mpatches.Patch(facecolor=global_colo_dict[x], label='Clostridiales Family XII.\nIncertae Sedis', edgecolor='black'))
        else:
             leg3_patches.append(mpatches.Patch(facecolor=global_colo_dict[x], label=x, edgecolor='black'))
    legend3=ax_br.legend(handles=leg3_patches[::-1],loc='upper left', bbox_to_anchor=(-0.185,-0.1), ncol=2,title="Taxonomy (≥"+str(metab_filter)+"%) (Ascending abundance)", fontsize=legendfont-0.5, frameon=False, labelspacing=0.1, columnspacing=columnssspacing, handletextpad=0.05)
    legend3.get_title().set_fontsize(10)
    ###letters
    font_dict= {'family': 'arial','color':  'black','weight': 'bold','size': 16}
    fig.text(0.1035,.8725, 'a', fontdict=font_dict)
    fig.text(0.1035,0.6, 'b', fontdict=font_dict)
    fig.text(0.515,.8725, 'c', fontdict=font_dict)
    fig.text(0.515,0.6, 'd', fontdict=font_dict)
    ####add colored backgrounds
    alphaz=0.1
    redv=(1,.95,.95)
    #(xy)z,m = (lowerleft), width, height
    bot=0.1
    fancybox1 = mpatches.FancyBboxPatch((0.12,bot), 0.365, 0.775,zorder=-1, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor='b',edgecolor='none',alpha=alphaz, transform=fig.transFigure)
    fancybox2 = mpatches.FancyBboxPatch((0.5275, bot), 0.155, 0.503,zorder=-1, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor='g',edgecolor='none',alpha=alphaz, transform=fig.transFigure)
    fancybox3 = mpatches.FancyBboxPatch((0.528,0.645), 0.361, 0.23,zorder=-1, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor=redv,edgecolor='none', transform=fig.transFigure)
    fancybox4 = mpatches.FancyBboxPatch((0.725,bot), 0.165, 0.55,zorder=-1, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor=redv,edgecolor='none', transform=fig.transFigure)
    fig.patches.append(fancybox1)
    fig.patches.append(fancybox2)
    fig.patches.append(fancybox3)
    fig.patches.append(fancybox4)
    return fig



filterx=1
df=tax_extract(merged_otu_table, tax_table, 'Class')
df_filtered_otu, df_filtered_normalised_class=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised_class, filterx)

fig=plot_my_df_strat_stack_stacked_LARGE_paper(df_taxa_grouped_x_col_master, master_color_dict, focus_groups_large, merge_dict, 'no', 'normalize_total', 1, df_filtered_normalised_class, 2, df_taxa_grouped_x_col_master_family, family_level_global_tax_col_dict, metab_list_in, 3)

#fig.savefig('KEGG_PLOT_BACKGROUNDS_750dpi_16_10pt5_SILVA.png', dpi=750, bbox_inches='tight', pad_inches = 0)



































####check what genes are involved
sulf_df=df_taxa_master[df_taxa_master['L2']=='Xenobiotics biodegradation and metabolism']
sulf_df['blastp']=sulf_df['index'].map(blastp_map_dict)
set(sulf_df['blastp'])
###sulfur enzymes are largely sulfite reductase, thiosulfate reductase, adenylyl-sulfate reductase
###nitrogen metabolism enzymes are largely glutamate dehydrogenase, glutamine synthetase, hydroxylamine reductase, nitrate reductase, nitrogenase


##############################OTU TABLE ################   
def tax_extract(otutable, taxtable, taxlevel):
    #Get an array of all unique taxa at specific level
    uniques_at_level=taxtable[taxlevel].unique()
    #Get the True False values for the index at these values, stored in a dictionary
    bool_dict={}
    for unique in uniques_at_level:
        bool_dict[unique]=taxtable[taxlevel]==unique
    #Sum the OTUs from the OTU table by passing the boolean filter to it
    #Store this information in new dataframe, with the uniques as the index
    df=pd.DataFrame(index=uniques_at_level, columns=otutable.columns)
    for unique in bool_dict:
        df.loc[unique]=np.sum(otutable[bool_dict[unique]])
    return df

def tax_normalise_filter(indf, filter_pct):
    dff=copy.deepcopy(indf)
    #filter by percentage
    taxastr='Taxa < ' + str(filter_pct)+"%"
    otherlist=[]
    for x in indf.columns:
        dff[x]=indf[x][((indf[x]/indf[x].sum())*100)>=filter_pct]    
        otherlist.append(indf[x][~(((indf[x]/indf[x].sum())*100)>=filter_pct)].sum(axis=0))
    row=pd.DataFrame([otherlist], columns=list(dff.columns), index=[taxastr])
    dff=dff.append(row)
    #remove empty index by sum
    dff=dff[dff.sum(axis=1)>0]
    #reorder columns as they originate from unordered dict
    #dff=dff[column_names]
    #normalise and keep otu df
    dfn=pd.DataFrame(index=dff.index, columns=dff.columns)
    for x in dff.columns:
        dfn[x]=((dff[x]/dff[x].sum())*100)
    return dff, dfn
 
def color_jenga(df):
    newpal=sns.xkcd_palette(["clear blue", "coral", "light aqua", "pale yellow", "light blue", "pinkish red"])
    length=int(len(df.index)/2)+1
    #plus 1 so as to not get caught out when int round 0.45 and 0.49 to 0...
    lightpal=sns.hls_palette(length, l=.5, s=.8, h=.5)
    darkpal=sns.hls_palette(length, l=.75, s=.9, h=0)
    for x in list(range(0, length)):
        newpal.append(lightpal[x])
        newpal.append(darkpal[x])
    random.shuffle(newpal)
    return newpal

def plot_my_df(indf, filterx):
    colors=color_jenga(indf)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(10,8))
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor='black', color=colors[x], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xticklabels(list(indf.columns))
    sns.despine(top=True, right=True)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Relative abundance (%)", fontsize=24)
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=1,title="Taxonomy (≥"+str(filterx)+"%)\n(Ascending abundance)", fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend

tax_table="Com_prof/tax_table_16S_2018.csv"
otu_table="Com_prof/otu_table_16S_2018.csv"

###From csv files
tax_table=pd.DataFrame.from_csv(tax_table)
otu_table=pd.DataFrame.from_csv(otu_table)

###Merge samples
column_names=['0', '1','2','3','4','5','6','8','10','16']
delimiter='-16S'
#delimiter='-ITS'
#delimiter='-18S'
letters=['A','B','C','D','E']
sample_dict={}
#Construct_sample_dict
for x in column_names:
    temp=[]
    if x=='0':
        sample_dict[x]=['D0'+delimiter]
    else:
        for y in letters:
            temp.append(x+y+delimiter)
        sample_dict[x]=temp
        
def merge_samples(table, sample_formats):
    newdf=pd.DataFrame()
    for x in sample_formats:
        newdf[x]=np.mean(table[sample_formats[x]], axis=1)
    return newdf

merged_otu_table=merge_samples(otu_table, sample_dict)
#re_order columns 
merged_otu_table=merged_otu_table[column_names]
#Fill 'nans' with 'NA'
tax_table=tax_table.fillna('NA')

filterx=2
df=tax_extract(merged_otu_table, tax_table, 'Class')
df_filtered_otu, df_filtered_normalised_class=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised_class, filterx)

##########################################REFERENCE KO TO BLAST HITS###########################
blastp_tophit='Proteomics/blastP/BlastP_result_files/ISDE_cat_blastp_out.xlsx'
f=xlrd.open_workbook(blastp_tophit)
sh=f.sheet_by_index(0)
blastp_map_dict={}
for rownum in range(0,sh.nrows):
        #ki=((''.join(reversed((''.join(reversed(k))).split(' ', 1)[0])).strip('()')))
        ki=sh.cell_value(rownum,0)
       # graphname='OUT_target_images/'+ki+'_molpct_graph_sterr.png'
       # vectorname='OUT_target_images/'+ki+'_file.png'
        blastp_map_dict[ki]=sh.cell_value(rownum,3)
        #target_ORF_dict[k]['BlastP_tophit']=sh.cell_value(rownum,3)
        #target_ORF_dict[k]['BlastP_e_value']=sh.cell_value(rownum,4)
        #target_ORF_dict[k]['BlastP_ID']=sh.cell_value(rownum,5)
        #target_ORF_dict[k]['BlastP_coverage']=sh.cell_value(rownum,6)
        #target_ORF_dict[k]['lineage']=taxa_master_dict[ki]['Taxonomy']['lineage']
       # target_ORF_dict[k]['totalrank']=df_master_mean.loc[ki]['rank']
       # plot_barplot(k)
        #target_ORF_dict[k]['graph']=graphname
        #target_ORF_dict[k]['vector']=vectorname



######################SULFUR STUFF####################
#Look at specific taxa 
dsub=df_taxa_master[df_taxa_master['taxa']=='Deltaproteobacteria']
#at a specific level
dsub2=dsub[dsub['pathway']=='Sulfur metabolism']
dsub2['blastp']=dsub['index'].map(blastp_map_dict)

########################PHOSP STUFF#######################
phosplist=[]
for x in blastp_map_dict:
    if bool(re.search('phosp', blastp_map_dict[x]))==True:
         print(x, blastp_map_dict[x])
         phosplist.append(x)
df_phosp=df_taxa_master[df_taxa_master['index'].isin(phosplist)]
df_phosp['blastp']=df_phosp['index'].map(blastp_map_dict)