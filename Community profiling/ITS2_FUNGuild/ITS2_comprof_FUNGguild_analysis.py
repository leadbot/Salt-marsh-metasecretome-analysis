# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 19:40:04 2020

@author: Dan
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import copy
import random
import pickle
####Importing and working with otu and tax tables from R

#tax_table="tax_table_16S_2018.csv"
#otu_table="otu_table_16S_2018.csv"

#tax_table="tax_table_18S_2018.csv"
#otu_table="otu_table_18S_2018.csv"

tax_table="tax_table_ITS2_2018.csv"
otu_table="otu_table_ITS2_2018.csv"
ITS_analysis='yes'
#NB silva 123 database has unusual taxonomy...
#Rank 1= Domain, Rank 2 =?, Rank 3=Superphyla?, Rank4=kingdom/phylum, Rank5='Class', Rank6='Order', Rank7='genus/family/species'



###From csv files
tax_table=pd.DataFrame.from_csv(tax_table)
otu_table=pd.DataFrame.from_csv(otu_table)

###Merge samples
column_names=['0', '1','2','3','4','5','6','8','10','16']
#delimiter='-16S'
delimiter='-ITS'
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

###format 18S index
def format_18S_index(df):
    newind=[]
    for x in df.index:
         if "_" in x:
             newind.append(x.split('_',3)[3])
         elif "-" not in x:
             newind.append(x)
    df=df.reset_index(drop=True)
    df['index']=newind
    df=df.set_index('index')
    return df

####################### ITS2 ############################################
####################### ITS2 ############################################
with open("Blast_all_unknown_hits_610.pickle", "rb") as input_file:
       all_blast_unknowns_dict=pickle.load(input_file)
   
with open("Blast_fungi_hits_only_351.pickle", "rb") as input_file:
       ITS2_fungi_unknowns_dict=pickle.load(input_file)
       
with open("Blast_offtarget_hits_only_42.pickle", "rb") as input_file:
       ITS2_offtarget_unknowns_dict=pickle.load(input_file)
test=[]            
if ITS_analysis.lower()=='yes':
    ##remove all unknowns from otu and tax table
    for otu in all_blast_unknowns_dict:
        if otu not in list(ITS2_fungi_unknowns_dict.keys()):
            if otu in list(tax_table.index):
                tax_table=tax_table.drop(otu, axis=0)
            if otu in list(otu_table.index):
                otu_table=otu_table.drop(otu, axis=0)
#drop additional otus          
tax_table=tax_table.drop(['OTU_1001','OTU_347'], axis=0)    
otu_table=otu_table.drop(['OTU_1001','OTU_347'], axis=0)                
###reassign column values with new taxa data
if ITS_analysis.lower()=='yes':
    ###drop unessesary columns from output
    droplist=['Rank1','Rank2','Rank3','Rank4','Rank5','Rank6','Rank7']
    for rnk in droplist:
         if rnk in tax_table.columns:
             tax_table=tax_table.drop(rnk, axis=1)
    ##change columns to lower
    cols=tax_table.columns
    colslower=[x.lower() for x in cols]
    tax_table.columns=colslower
    ##create replacement dict as {column: {index:value}}
    for otu in ITS2_fungi_unknowns_dict:
         for rank in ITS2_fungi_unknowns_dict[otu]:
             if not rank =='superkingdom':
                   tax_table.ix[otu,rank]=ITS2_fungi_unknowns_dict[otu][rank]
    #return columns to normal
    tax_table.columns=cols
tax_table=tax_table.replace(to_replace='unidentified', value='NA')
####################### ITS2 ############################################
####################### ITS2 ############################################
#Ensure data is float not str
otu_table=otu_table.astype(float)
merged_otu_table=merge_samples(otu_table, sample_dict)
#re_order columns 
merged_otu_table=merged_otu_table[column_names]
#Fill 'nans' with 'NA'
tax_table=tax_table.fillna('NA')
######Make output for funguild#########
FUNGuild_df=copy.deepcopy(otu_table)
FUNGuild_df=copy.deepcopy(merged_otu_table)
mapdict={}
for otu in list(merged_otu_table.index):
    linstr=''
    for tr in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']:
        prefix=tr[0].lower()+'_'
        if not tax_table.ix[otu][tr] =='NA':
            linstr=linstr+prefix+tax_table.ix[otu][tr]
        elif tax_table.ix[otu][tr] =='NA':
            linstr=linstr+prefix+'unidentified'
        if not tr =='Species':
            linstr=linstr+';'
    mapdict[otu]=linstr        
FUNGuild_df['taxonomy']=pd.Series(mapdict)
FUNGuild_df.to_csv("FUNGuild_input_file.tsv", header=True, sep='\t', index=True, index_label='OTU ID')
####################### ITS2 Read FUNGuild output ############################################
####################### ITS2 Read FUNGuild output ############################################
f=open('FUNGuild_output_file.guilds_matched.txt', 'r')
#NB the output helpfully has no linebreaks so need to write a function manually to parse the output txt
##There are 22 headers ([0:21]) and the 22nd is joined to the first column of the next line..
line=f.readline().strip('\n')
f.close()



splitdata=line.split('\t')
columnsf=[]
globalposition=0
for row in list(range(0,420)):
    if row == 0:
         for position in list(range(0,21)):
              if position==20:
                  columnsf.append(splitdata[position].split('OTU_')[0])
              if not position ==20:
                  columnsf.append(splitdata[position])
         fungdf=pd.DataFrame(columns=columnsf)
         globalposition+=20
    if not row==0 and not globalposition>=8400:
        rowf=[]
        for position in list(range(0,21)):
            if (globalposition+position) % 20 ==0 and position ==0:
                print('entereerereredloop')
                startstr='OTU_'+splitdata[globalposition+position].split('OTU_')[1]
                rowf.append('OTU_'+splitdata[globalposition+position].split('OTU_')[1])
            if (globalposition+position) % 20 !=0:
                print(globalposition, position)
                rowf.append(splitdata[globalposition+position])
            if (globalposition+position) % 20 ==0 and position ==20:
                rowf.append('OTU_'+splitdata[globalposition+position].split('OTU_')[0])
        print(rowf)
        print(startstr)
        globalposition+=20
        fungdf=fungdf.append(pd.DataFrame([rowf],columns=columnsf), ignore_index=True)   
##########################FORMAT DF FOR ANALYSIS###################################
fungdfformat=fungdf[['0', '1', '2', '3', '4', '5', '6', '8', '10', '16','Trophic Mode', 'Guild']]
fungdfRAW=fungdf[['0', '1', '2', '3', '4', '5', '6', '8', '10', '16','Trophic Mode', 'Guild', 'taxonomy', 'OTU ID']]
fungdfRAW=fungdfRAW.astype({'0':float, '1':float, '2':float, '3':float, '4':float, '5':float,\
                '6':float, '8':float, '10':float, '16':float})
its2_filter_percentage=0.05
fungdfformatorig=copy.deepcopy(fungdfformat)
otherlist=[]
for timepoint in ['0', '1', '2', '3', '4', '5', '6', '8', '10', '16']:
    fungdfformat[timepoint] = pd.to_numeric(fungdfformat[timepoint], downcast="float")
    fungdfformatorig[timepoint] = pd.to_numeric(fungdfformatorig[timepoint], downcast="float")
    #normalise and filter otus before groupby
    ###########
sum_column=fungdfformatorig.sum(axis=1)
sum_column_mask=(sum_column/sum_column.sum())*100>its2_filter_percentage
fungdfformat=fungdfformat[sum_column_mask]
otherlist_colfilt=fungdfformatorig[~sum_column_mask][['0','1','2','3','4','5','6','8','10','16']].sum(axis=0)
row_colfilt=pd.DataFrame([list(otherlist_colfilt)+['Modes < '+str(its2_filter_percentage)+'%']*2], columns=list(fungdfformat.columns))
fungdfformat=fungdfformat.append(row_colfilt)  
    #############
    #fungdfformat[timepoint]=fungdfformatorig[timepoint][((fungdfformatorig[timepoint]/fungdfformatorig[timepoint].sum())*100)>=its2_filter_percentage]    
    #otherlist.append(fungdfformatorig[timepoint][~(((fungdfformatorig[timepoint]/fungdfformatorig[timepoint].sum())*100)>=its2_filter_percentage)].sum(axis=0))
#row=pd.DataFrame([otherlist+['OTU > '+str(its2_filter_percentage)+'%']*2], columns=list(fungdfformat.columns))
#fungdfformat=fungdfformat.append(row)  
fungdfformat_grp=fungdfformat.groupby(['Trophic Mode', 'Guild']).sum()
#remove empty index by sum
fungdfformat_grp=fungdfformat_grp[fungdfformat_grp.sum(axis=1)>0]
#normalise to percent
fungdfformat_grp_nrm=copy.deepcopy(fungdfformat_grp)
for x in fungdfformat_grp.columns:
    fungdfformat_grp_nrm[x]=((fungdfformat_grp[x]/fungdfformat_grp[x].sum())*100)

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
    newpal=sns.xkcd_palette(["light grey", "white","clear blue", "coral", "light aqua", "pale yellow", "light blue", "pinkish red"])
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
    labels=[str(h).replace('_fam_Incertae_sedis','').replace('_cls_Incertae_sedis','') for h in labels]
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=1,title="Taxonomy (≥"+str(filterx)+"%)\n(Ascending abundance)", fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend

def plot_my_df_its2(indf, filterx):
    colors=color_jenga(indf)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(8,8))
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    '-', '+', 'x', '\\', '*', 'o', 'O', '.'
    hatchdict={'Saprotroph':'-', 'Pathotroph':'+', 'Pathotroph-Symbiotroph':'x', 'Pathotroph-Saprotroph-Symbiotroph':'\\','Pathotroph-Saprotroph':'*','Modes < 0.05%':''}
    for x in list(range(0, len(indf.index))):
        if indf.iloc[x].name[0] in hatchdict:
             ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor='black', hatch=hatchdict[indf.iloc[x].name[0]],color=colors[x], label=indf.iloc[x].name)
             yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
        elif not indf.iloc[x].name[0] in hatchdict:
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
    ###order the legend handle and labels:
    ha=[]
    la=[]
    alphalab=sorted(labels)
    for al in alphalab:
         for h, l in zip(handles, labels):
             if al==l:
                 ha.append(h)
                 la.append(l)
    newla=[]
    for l in la:
        newla.append(l.strip("()'").replace("'",'').replace(",",";").split('; ')[1])
        #if l=='Modes < 0.05%':
         #   newla.append('Modes < 0.05%')
        #elif not l =='Modes < 0.05%':
         #   newla.append(l.split('; ')[1])
    replace_dict={'Modes < 0.05%':'Modes < 0.05%',\
 'Fungal Parasite':'Fungal Parasite',\
 'Plant Pathogen':'Plant pathogen',\
 'Animal Pathogen-Fungal Parasite-Undefined Saprotroph':'Animal Pathogen-Fungal Parasite\nUndefined Saprotroph',\
 'Endophyte-Lichen Parasite-Plant Pathogen-Undefined Saprotroph':'Endophyte-Lichen Parasite\nPlant Pathogen-Undefined Saprotroph',\
 'Fungal Parasite-Litter Saprotroph':'Fungal Parasite-Litter Saprotroph',\
 'Fungal Parasite-Plant Pathogen0Plant Saprotroph':'Fungal Parasite-Plant Pathogen\nPlant Saprotroph',\
 'Animal Pathogen-Endophyte-Fungal Parasite-Plant Pathogen-Wood Saprotroph':'Animal Pathogen-Endophyte-Fungal Parasite\nPlant Pathogen-Wood Saprotroph',\
 'Animal Pathogen-Endophyte-Lichen Parasite-Plant Pathogen-Soil Saprotroph-Wood Saprotroph':'Animal Pathogen-Endophyte-Lichen Parasite\nPlant Pathogen-Soil Saprotroph-Wood Saprotroph',\
 'Animal Pathogen-Endophyte-Plant Pathogen-Wood Saprotroph':'Animal Pathogen-Endophyte-Plant Pathogen\nWood Saprotroph',\
 'Fungal Parasite-Undefined Saprotroph':'Fungal Parasite-Undefined Saprotroph',\
 'Endophyte-Plant Pathogen':'Endophyte-Plant Pathogen',\
 'Dung Saprotroph-Plant Saprotroph-Wood Saprotroph':'Dung Saprotroph-Plant Saprotroph\nWood Saprotroph',\
 'Leaf Saprotroph':'Leaf Saprotroph',\
 'Undefined Saprotroph':'Undefined Saprotroph',\
 'Wood Saprotroph':'Wood Saprotroph',\
 'Fungal Parasite-Plant Pathogen-Plant Saprotroph':'Fungal Parasite-Plant Pathogen\nPlant Saprotroph'}        
    newl=[]
    for l2 in newla:
        if l2 in replace_dict.keys():
            newl.append(replace_dict[l2])
    legend=ax.legend(ha[::-1], newl[::-1], loc='upper left', handlelength=2,handleheight=2,\
                     bbox_to_anchor=(1.1,1.025), ncol=1,fontsize=14)
    legend.get_title().set_fontsize(20)
    #plt.annotate("figure fraction",
     #       xy=(0.9,0.5), xycoords='figure fraction',
      #      xytext=(0.9, 0.8), textcoords='figure fraction',
       #     arrowprops=dict(arrowstyle="-",
        #              connectionstyle="arc3, rad=0"),
         #   )
    xvx=0.6075 
    xvv=xvx-0.05
    ax.annotate('', xy=(xvx, 0.93), xycoords='figure fraction',
             xytext=(xvx, 0.73), textcoords='figure fraction', 
             ha="center", va="center",
             arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=4))
    ax.annotate('Saprotroph', xy=(xvv+0.03, 0.88), rotation=90, size=14,xycoords='figure fraction')
    ax.annotate('', xy=(xvx, 0.66), xycoords='figure fraction',
             xytext=(xvx, 0.43), textcoords='figure fraction',
             ha="center", va="center",
             arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=4))
    ax.annotate('Pathotroph\nSaprotroph\nSymbiotroph', xy=(xvv, 0.54), rotation=90, size=14,xycoords='figure fraction')
    ax.annotate('', xy=(xvx, 0.4), xycoords='figure fraction',
             xytext=(xvx, 0.18), textcoords='figure fraction',
             ha="center", va="center",
             arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=4)) 
    ax.annotate('Pathotroph\nSaprotroph', xy=(xvv+0.015, 0.3), rotation=90, size=14,xycoords='figure fraction')

#    ax.annotate('', xy=(0.9, 0.175), xycoords='figure fraction',
 #            xytext=(0.9, 0.1), textcoords='figure fraction',
  #           ha="center", va="center",
   #          arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=4)) 
    #ax.annotate('Pathotroph', xy=(0.9, 0.1375), rotation=270, size=14,xycoords='figure fraction')
    return fig, legend, newl, newla, handles, labels,ha, la

figf, legf, newl, newla, handles, labels, h, k=plot_my_df_its2(fungdfformat_grp_nrm, its2_filter_percentage)
#figf.savefig("Trophic_modes_ITS2_final_fig.png", dpi=1200, bbox_inches='tight')
#
####assess changes in guild and trophic mode over time
fungdfformat_grp_tm=fungdfformat.groupby(['Trophic Mode']).sum()
fungdfformat_grp_guild=fungdfformat[['0', '1', '2', '3', '4', '5', '6', '8', '10', '16','Guild']].groupby(['Guild']).sum()
fungdfformat_grp_tm=fungdfformat_grp_tm[fungdfformat_grp_tm.sum(axis=1)>0]
fungdfformat_grp_guild=fungdfformat_grp_guild[fungdfformat_grp_guild.sum(axis=1)>0]
#normalise to percent
fungdfformat_grp_tm_nrm=copy.deepcopy(fungdfformat_grp_tm)
fungdfformat_grp_guild_nrm=copy.deepcopy(fungdfformat_grp_guild)
for x in fungdfformat_grp_tm.columns:
    fungdfformat_grp_tm_nrm[x]=((fungdfformat_grp_tm[x]/fungdfformat_grp_tm[x].sum())*100)
for x in fungdfformat_grp_guild.columns:
    fungdfformat_grp_guild_nrm[x]=((fungdfformat_grp_guild[x]/fungdfformat_grp_guild[x].sum())*100)
    
figtm, legtm=plot_my_df(fungdfformat_grp_tm_nrm, its2_filter_percentage)
figg, legg=plot_my_df(fungdfformat_grp_guild_nrm, its2_filter_percentage)


def give_me_trophic_mode_from_taxa(df, taxa):
    ###best to use fungdfRAW as df
    alltaxlist=set(df[df['taxonomy'].str.contains(str(taxa))]['taxonomy'])
    g=set(df[df['taxonomy'].str.contains(str(taxa))]['Guild'])
    tm=set(df[df['taxonomy'].str.contains(str(taxa))]['Trophic Mode'])
    gdf=df[df['taxonomy'].str.contains(str(taxa))][['0', '1', '2', '3', '4', '5', '6', '8', '10', '16', 'Guild']].groupby('Guild').sum()
    tmdf=df[df['taxonomy'].str.contains(str(taxa))][['0', '1', '2', '3', '4', '5', '6', '8', '10', '16', 'Trophic Mode']].groupby('Trophic Mode').sum()
    gdff, gdfn=tax_normalise_filter(gdf, 0)
    tmdff, tmdfn=tax_normalise_filter(tmdf,0)
    plot_my_df(gdfn, 0)
    plot_my_df(tmdfn,0)
    print(alltaxlist)
    print('Guilds: '+str(g))
    print('Trophic modes: ' + str(tm))
    return gdfn, tmdfn
give_me_trophic_mode_from_taxa(fungdfRAW, 'Pleosporaceae')  
give_me_trophic_mode_from_taxa(fungdfRAW, 'Saccharomycetales')    
give_me_trophic_mode_from_taxa(fungdfRAW, 'Hypocreales')  
give_me_trophic_mode_from_taxa(fungdfRAW, 'Capnodiales')  
give_me_trophic_mode_from_taxa(fungdfRAW, 'Tremellales')  
###############################OTU RICHNESS########################################
samples=[1,2,3,4,5,6,8,10,16]
otucountdict={}

def otu_richness(indf, samples, otucountdict):
    letters=['A','B','C','D','E']
    counts=otu_table[otu_table>0].count(axis=0)
    for x in samples:
        if x not in otucountdict:
            otucountdict[x]=[]
        for y in counts.index:
            for z in letters:
                 if str(y).startswith(str(x)+z):
                     otucountdict[x].append(counts[y])
    otucountdict[0]=[counts['D0'+str(delimiter)]]             
    return otucountdict
 
otu_richness(otu_table, samples, otucountdict)

data_to_plot=[]
sample_order=[0,1,2,3,4,5,6,8,10,16]
for x in sample_order:
    data_to_plot.append(otucountdict[x])

colopal=sns.hls_palette(10)
fig=plt.figure(1, figsize=(12, 8))

ax=fig.add_subplot(111)
bp=ax.boxplot(data_to_plot, patch_artist=True)
sns.set_style('white')
ax.set_ylabel('OTU richness (n)', fontsize=22)
ax.set_xlabel('Week', fontsize=22)
ax.spines['top'].set_visible(False)
ax.set_xticklabels(['0', '1', '2','3','4', '5','6','8','10','16'], fontsize=20)
ax.spines['right'].set_visible(False)
ax.set_ylim([0, 2500])
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=20)
for box in bp['boxes']:
    box.set(color='#000000', linewidth=2)
#set first line a bit wider to make more visible#
bp['boxes'][0].set(color='#000000', linewidth=5)
for patch, color in zip(bp['boxes'], colopal):
    patch.set(facecolor=color)
for whisker in bp['whiskers']:
    whisker.set(color='#000000', linewidth=2, linestyle='solid')
for cap in bp['caps']:
    cap.set(color='#000000', linewidth=2)
for median in bp['medians']:
    median.set(color='#000000', linewidth=2)
for flier in bp['fliers']:
    flier.set(marker='o', alpha=0.5)
    flier.set(color='#6aeae9')
    flier.set(markeredgecolor='#6aeae9')
    
#fig.savefig('OTU_richness_'+str(delimiter)+'.png', bbox_inches='tight', dpi=1200)
###########################OTU RICHNESS ##################################################

    
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

def within_taxa_filter_then_normalise(taxa_to_inspect, taxlevel, otutable, taxtable):
    taxdf=tax_extract(otutable, taxtable, taxlevel)
    dff=taxdf.copy(deep=True)
    for x in taxdf.columns:
        dff[x]=taxdf[x][((taxdf[x]/taxdf[x].sum())*100)==taxa_to_inspect]
    return dff

###
filterx=0.0001
df=tax_extract(merged_otu_table, tax_table, 'Kingdom')
#df=format_18S_index(df)
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)
#fig.savefig("Barplots/Kingdom_barplot_high_level2.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
###

###
filterx=0
df=tax_extract(merged_otu_table, tax_table, 'Phylum')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)
fig.savefig("Barplots/Phylum_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
###

###
filterx=0.25
df=tax_extract(merged_otu_table, tax_table, 'Class')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)
fig.savefig("Barplots/Class_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)

###

###
filterx=0.5
df=tax_extract(merged_otu_table, tax_table, 'Order')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)
fig.savefig("Barplots/Order_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
###

###
filterx=1
df=tax_extract(merged_otu_table, tax_table, 'Family')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)
fig.savefig("Barplots/Family_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
###

def plot_my_dfx2(indf, filterx,prelim, otudata):

    colors=color_jenga(indf)
    #reorder the dataframe by abundance (make it easier to plot)

    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)

    indf=indf.fillna(value=0)

    fig, ax = plt.subplots(figsize=(6,6))

    width=0.8 

    #Loop over index and plot the values

    yvalue_cumulative=np.array([0]*len(list(indf.columns)))

    position=list(range(0, len(indf.columns)))

    ##Remove 'NA' column##

    #indf=indf.loc[~(indf.index=='NA')]

    ### Plot the amount of data plotted

    ####Plot bars

    for x in list(range(0, len(indf.index))):

        #if indf.iloc[x].name in coldict:

             ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=colors[x],label=indf.iloc[x].name)

             yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])#color=coldict[indf.iloc[x].name],

        #else:

           # ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=resids[x], label=indf.iloc[x].name)

          #  yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])

    ax.set_ylim([0,100.5])

    ax.set_xticks(position)

    ax.set_xlim([0-width/2, max(position)+width/2+0.05])

    ax.set_xticklabels(list(indf.columns))

    sns.set(style="white")

    ax.spines['top'].set_visible(False)

    #ax2.spines['top'].set_visible(False)

    ax.tick_params(direction='out', length=6, width=4)  

    ax.tick_params(axis='both', labelsize=20)

    ax.set_ylabel("Relative abundance (%)", fontsize=24)

    ####boxplot#####

    colopal=sns.xkcd_palette(['clear blue'])*10

    ax2=ax.twinx()

    bp=ax2.boxplot(otudata, patch_artist=True, positions=position)

    ax2.set_ylabel('OTU richness (n)', fontsize=22)

    ax2.tick_params(direction='out', length=6, width=4)  

    ax2.tick_params(axis='both', labelsize=20)

    ax2.set_ylim([0,400])

    ax2.spines['top'].set_visible(False)

    ax2.set_xticks(position)

    ax2.set_xticklabels(list(indf.columns))

    #####box visuals###

    for box in bp['boxes']:

         box.set(color='#000000', linewidth=2)

    #set first line a bit wider to make more visible#

         bp['boxes'][0].set(color='#000000', linewidth=5)

    for patch, color in zip(bp['boxes'], colopal):

         patch.set(facecolor=color)

    for whisker in bp['whiskers']:

         whisker.set(color='#000000', linewidth=2, linestyle='solid')

    for cap in bp['caps']:

         cap.set(color='#000000', linewidth=2)

    for median in bp['medians']:

         median.set(color='#000000', linewidth=2)

    #for flier in bp['fliers']:

      #  flier.set(marker='o', alpha=0.5)

      #  flier.set(color='#6aeae9')

      #  flier.set(markeredgecolor='#6aeae9')

    #ax2.set_ylim([0,1.3])

    ax.set_xlabel('Week', fontsize=24)

    #call handles and labels, so the legend is inverted matching the vertical bar structure

    handles, labels = ax.get_legend_handles_labels()

    h2, l2 = ax2.get_legend_handles_labels()

    handles=h2+handles

    labels=l2+labels
    labels=[str(h).replace('_fam_Incertae_sedis','').replace('_cls_Incertae_sedis','').replace('_ord_Incertae_sedis','') for h in labels]


    legend=ax.legend(handles[::-1], labels[::-1], loc=(1.3, 0.7), ncol=1,title="Fungi "+str(prelim), fontsize=10)

    #legend.get_title().set_fontsize(20)

    plt.tight_layout()

    return fig#,legend
#Fam = yval = 0.25

filterx=0.3
df=tax_extract(merged_otu_table, tax_table, 'Order')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)

fig =plot_my_dfx2(df_filtered_normalised, filterx, 'Order', data_to_plot)

#fig.savefig('Fungi_FAMILIES_and_OTU_richness.png', dpi=750, bbox_inches='tight')

#fig.savefig("Fungi_ORDERS_and_OTU_richness.png", dpi=750, bbox_inches='tight')

df_u=tax_extract(otu_table, tax_table, 'Order')
df_fo,df_um=tax_normalise_filter(df_u, filterx)
tlist=[]
targettt='Microascales'
for x in [6,7,8,9]:
    for y in ['A','B','C','D','E']:
         tlist.append(str(x)+str(y)+'-ITS')
print(df_um.loc[targettt][tlist].mean())
print(df_um.loc[targettt][tlist].sem())

give_me_trophic_mode_from_taxa(fungdfRAW, 'Microascales')  
give_me_trophic_mode_from_taxa(fungdfRAW, 'Pleosporales')   
give_me_trophic_mode_from_taxa(fungdfRAW, 'Saccharomycetales')    
give_me_trophic_mode_from_taxa(fungdfRAW, 'Hypocreales')  
give_me_trophic_mode_from_taxa(fungdfRAW, 'Capnodiales')   
give_me_trophic_mode_from_taxa(fungdfRAW, 'Tremellales')  

###How much, in terms of abundance, did FUNGuild classify?
#normalise whole otu table
norm_otutabledf=pd.DataFrame(index=merged_otu_table.index, columns=merged_otu_table.columns)
for x in merged_otu_table.columns:
        norm_otutabledf[x]=((merged_otu_table[x]/merged_otu_table[x].sum())*100)

supdf=norm_otutabledf.loc[list(fungdfRAW['OTU ID'])]

fungdfformat_grp_nrm.index[8]

fungdfformat_grp_nrm.loc[('Pathotroph-Saprotroph-Symbiotroph',
 'Animal Pathogen-Endophyte-Lichen Parasite-Plant Pathogen-Soil Saprotroph-Wood Saprotroph')][[1,2,3]].mean()

fungdfformat_grp_nrm.loc[('Pathotroph-Saprotroph',
 'Endophyte-Lichen Parasite-Plant Pathogen-Undefined Saprotroph')][[2,3, 4, 5, 6, 7, 8, 9]].mean()

fungdfreps=otu_table.loc[fungdf['OTU ID']]
fungdfreps_nrm=pd.DataFrame(index=fungdfreps.index, columns=fungdfreps.columns)
for x in fungdfreps.columns:
        fungdfreps_nrm[x]=((fungdfreps[x]/fungdfreps[x].sum())*100)
fungdfreps_nrm['Guild']=fungdfRAW.set_index('OTU ID')['Guild']
fungdfreps_nrm_guildgrp=fungdfreps_nrm.groupby('Guild').sum()
fungdfreps_nrm_guildgrp[['16A-ITS','16B-ITS','16C-ITS','16D-ITS','16E-ITS']]['Undefined Saprotroph'].sem()