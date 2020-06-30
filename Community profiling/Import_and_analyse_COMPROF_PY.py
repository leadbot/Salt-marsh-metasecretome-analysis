# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 11:24:44 2018

@author: Dan
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import copy
import random
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec


otu_table_file="otu_table_16S_2018.csv"
silva_tax_table="formatted_16s_OTUs_2018_tax_assignments_SILVA_132.txt"

#tax_table="tax_table_18S_2018.csv"
#otu_table="otu_table_18S_2018.csv"
#
tax_table="tax_table_ITS2_2018.csv"
otu_table="otu_table_ITS2_2018.csv"


def create_dict_of_otus_at_tax_res_for_picrust(indf, levels):
    outdf={}
    for x in levels:
        if x not in outdf:
            outdf[x]={}
        uniques=set(list(tax_table[x]))
        for tax in uniques:
            if tax not in outdf[x]:
                outdf[x][tax]=list(tax_table[tax_table[x]==tax].index)
    return outdf
            
#tax_otu_master_dict=create_dict_of_otus_at_tax_res_for_picrust(tax_table, ['Phylum', 'Class', 'Family'])


###From csv files
#tax_table=pd.read_csv(tax_table_file, index_col=0)

################SILVA ANALYSIS######################
################SILVA ANALYSIS######################
################SILVA ANALYSIS######################
tax_template=pd.read_csv('formatted_16s_OTUs_2018_tax_assignments_SILVA_132.txt', sep='\t', index_col=0, header=None)

tax_table=pd.DataFrame(index=tax_template.index)
tax_table['Kingdom']=tax_template.iloc[:,0].str.split(';',1).str[0].str.split('__').str[1]
tax_table['Phylum']=tax_template.iloc[:,0].str.split(';',2).str[1].str.split('__').str[1]
tax_table['Class']=tax_template.iloc[:,0].str.split(';',3).str[2].str.split('__').str[1]
tax_table['Order']=tax_template.iloc[:,0].str.split(';',4).str[3].str.split('__').str[1]
tax_table['Family']=tax_template.iloc[:,0].str.split(';',5).str[4].str.split('__').str[1]
tax_table['Genus']=tax_template.iloc[:,0].str.split(';',6).str[5].str.split('__').str[1]
tax_table['Species']=tax_template.iloc[:,0].str.split(';',7).str[6].str.split('__').str[1]
################SILVA ANALYSIS######################
################SILVA ANALYSIS######################
################SILVA ANALYSIS######################
otu_table=pd.read_csv(otu_table_file, index_col=0)

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


###############################OTU RICHNESS########################################
samples=[1,2,3,4,5,6,8,10,16]
otucountdict={}

def otu_richness(indf, samples, otucountdict):
    otutuples=[]
    letters=['A','B','C','D','E']
    counts=indf[indf > 0].count(axis=0)
    for x in samples:
        if x not in otucountdict:
            otucountdict[x]=[]
        for y in counts.index:
            for z in letters:
                 if str(y).startswith(str(x)+z):
                     otucountdict[x].append(counts[y])
                     otutuples.append((x,counts[y]))
    otucountdict[0]=[counts['D0-16S']]             
    return otucountdict, otutuples

otucountdict={} 
otucountdict, otutuples=otu_richness(otu_table, samples, otucountdict)


data_to_plot=[]
sample_order=[0,1,2,3,4,5,6,8,10,16]
for x in sample_order:
    data_to_plot.append(otucountdict[x])

#colopal=sns.hls_palette(10)
colopal=sns.xkcd_palette(['clear blue'])*10
fig=plt.figure(1, figsize=(12, 8))

ax=fig.add_subplot(111)
bp=ax.boxplot(data_to_plot, patch_artist=True)
sns.set_style='white'
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
    
fig.savefig('OTU_richness_16S.png', bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all') 
###########################OTU RICHNESS ##################################################
def plot_my_df_cazy_ncol2(indf, filterx, heatmap_df, taxa, level, supplement, invaried, intotal, lAvel, prelim, varcel, varcelna):
    cazycol=[]
    noncazycol=[]
    for x in list(heatmap_df.index):
         if x in list(indf.index):
              cazycol.append(x)
         if x not in list(indf.index):
              noncazycol.append(x.strip('\n'))
    cazycol=cazycol+supplement
    indf=indf.loc[cazycol]
    #remove 'NA'
    indf=indf.loc[~(indf.index=='NA')]
    heatmap_df=heatmap_df.loc[~(heatmap_df.index=='NA')]
    molpct_before_filter=heatmap_df.sum(axis=0)
    #filter heatmap to contain only taxa found it amplicon sequencing
    print(len(heatmap_df.index))
    heatmap_df=heatmap_df.loc[cazycol]
    print(len(heatmap_df.index))
    #
    colors=color_jenga(indf)
    molpct_to_plot=(heatmap_df.loc[cazycol]).sum(axis=0)
    unaccounted_mol_pct=np.array(molpct_before_filter)-np.array(molpct_to_plot)
    molpct_err=(heatmap_df.loc[cazycol]).sem(axis=0)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(10,8))
    ax2=ax.twinx()
    ax2.errorbar(np.array([1,3,5,8]), np.array(intotal), fmt='-o', color='blue', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='CAZy total')
    ax2.errorbar(np.array([1,3,5,8]), np.array(varcelna), fmt='-o', color='red', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='Visible + Cellvibrionacea + NA')
    ax2.errorbar(np.array([1,3,5,8]), np.array(varcel), fmt='-o', color='grey', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='Visible + Cellvibrionacea')
    ax2.errorbar(np.array([1,3,5,8]), np.array(invaried), fmt='-o', color='black', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='CAZy contributions (visible taxa)')
    
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ##Remove 'NA' column##
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=0.75, edgecolor='black', color=colors[x], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_xticklabels(list(indf.columns))
    sns.set(style="white")
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Relative abundance (%)", fontsize=24)
    ax2.tick_params(direction='out', length=6, width=4)  
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_ylabel(u"\u2211" r'$\bar{x}$'" Molar percentage", fontsize=24)
    ax2.set_ylim([0,1.3])
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    handles=h2+handles
    labels=l2+labels
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.125,1.025), ncol=2,title="CAZy producing "+str(prelim)+"\n      At "+str(lAvel)+" level\n(Ascending abundance)", fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend, indf, noncazycol

def plot_my_df_cazy_ncol_fam(indf, filterx, heatmap_df, taxa, level, supplement):
    cazycol=[]
    noncazycol=[]
    for x in list(heatmap_df.index):
         if x in list(indf.index):
              cazycol.append(x)
         if x not in list(indf.index):
              noncazycol.append(x.strip('\n'))
    cazycol=cazycol+supplement
    indf=indf.loc[cazycol]
    #remove 'NA'
    indf=indf.loc[~(indf.index=='NA')]
    heatmap_df=heatmap_df.loc[~(heatmap_df.index=='NA')]
    molpct_before_filter=heatmap_df.sum(axis=0)
    #filter heatmap to contain only taxa found it amplicon sequencing
    print(len(heatmap_df.index))
    heatmap_df=heatmap_df.loc[cazycol]
    print(len(heatmap_df.index))
    #
    colors=color_jenga(indf)
    molpct_to_plot=(heatmap_df.loc[cazycol]).sum(axis=0)
    unaccounted_mol_pct=np.array(molpct_before_filter)-np.array(molpct_to_plot)
    molpct_err=(heatmap_df.loc[cazycol]).sem(axis=0)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(10,8))
    ax2=ax.twinx()
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ##Remove 'NA' column##
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=0.75, edgecolor='black', color=colors[x], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_xticklabels(list(indf.columns))
    sns.set(style="white")
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Relative abundance (%)", fontsize=24)
    ax2.tick_params(direction='out', length=6, width=4)  
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_ylabel(u"\u2211" r'$\bar{x}$'" Molar percentage", fontsize=24)
    ax2.set_ylim([0,1.3])
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    handles=h2+handles
    labels=l2+labels
    legend.get_title().set_fontsize(20)
    return fig, legend, indf

def plot_my_df_cazy_ncol_class(indf, filterx, heatmap_df, taxa, level):
    cazycol=[]
    noncazycol=[]
    for x in list(heatmap_df.index):
         if x in list(indf.index):
              cazycol.append(x)
         if x not in list(indf.index):
              noncazycol.append(x.strip('\n'))
    cazycol=cazycol
    indf=indf.loc[cazycol]
    #remove 'NA'
    indf=indf.loc[~(indf.index=='NA')]
    heatmap_df=heatmap_df.loc[~(heatmap_df.index=='NA')]
    molpct_before_filter=heatmap_df.sum(axis=0)
    #filter heatmap to contain only taxa found it amplicon sequencing
    print(len(heatmap_df.index))
    heatmap_df=heatmap_df.loc[cazycol]
    print(len(heatmap_df.index))
    #
    colors=color_jenga(indf)
    molpct_to_plot=(heatmap_df.loc[cazycol]).sum(axis=0)
    unaccounted_mol_pct=np.array(molpct_before_filter)-np.array(molpct_to_plot)
    molpct_err=(heatmap_df.loc[cazycol]).sem(axis=0)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(10,8))
    ax2=ax.twinx()
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ##Remove 'NA' column##
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=0.75, edgecolor='black', color=colors[x], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_xticklabels(list(indf.columns))
    sns.set(style="white")
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Relative abundance (%)", fontsize=24)
    ax2.tick_params(direction='out', length=6, width=4)  
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_ylabel(u"\u2211" r'$\bar{x}$'" Molar percentage", fontsize=24)
    ax2.set_ylim([0,1.3])
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    handles=h2+handles
    labels=l2+labels
    legend.get_title().set_fontsize(20)
    return fig, legend, indf


#####This barplot takes the heatmap from ###########################TAXA clustermap on ISDe main script
###To filter to only CAZy producing taxa and display the molar percent they produce
def plot_my_df_cazy(indf, filterx, heatmap_df, taxa, level, supplement):
    cazycol=[]
    for x in list(heatmap_df.index):
         if x in list(indf.index):
              cazycol.append(x)
    cazycol=cazycol+supplement
    indf=indf.loc[cazycol]
    colors=color_jenga(indf)
    molpct_to_plot=(heatmap_df.loc[cazycol]).sum(axis=0)
    molpct_err=(heatmap_df.loc[cazycol]).sem(axis=0)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(10,8))
    ax2=ax.twinx()
    ax2.errorbar(np.array([1,3,5,8]), np.array(molpct_to_plot), yerr=np.array(molpct_err), fmt='-o', color='black', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='CAZy (molar percentage)')

    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ##Remove 'NA' column##
    indf=indf.loc[~(indf.index=='NA')]
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor='black', color=colors[x], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_xticklabels(list(indf.columns))
    sns.set(style="white")
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Relative abundance (%)", fontsize=24)
    ax2.tick_params(direction='out', length=6, width=4)  
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_ylabel(u"\u2211" r'$\bar{x}$'" Molar percentage", fontsize=24)
    ax2.set_ylim([0,1.2])
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.125,1.025), ncol=1,title="CAZy producing\n"+str(taxa) +" (≥"+str(filterx)+"%)\nAt "+str(level)+" level\n(Ascending abundance)", fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend, molpct_to_plot


    
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

def within_taxa_filter_then_normalise(taxa_to_inspect, taxlevel, otutable, taxtable):
    taxdf=tax_extract(otutable, taxtable, taxlevel)
    dff=taxdf.copy(deep=True)
    for x in taxdf.columns:
        dff[x]=taxdf[x][((taxdf[x]/taxdf[x].sum())*100)==taxa_to_inspect]
    return dff
  
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
    
    
#Built around 16, 12 size, resizing will require new formatting
#def plot_my_df_with_line(indf, filterx, toggle, reltot, reltoterr, linelevel, focaltax, highlevelplot, highlevelx, highlevelxerr):
def plot_my_df_with_line(indf, filterx, toggle, linelevel, focaltax, indict):
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
    linecolors=['black', 'blue', 'green', 'red']
    if toggle =='yes':
        colocount=0
        for x in indict:
            ax.errorbar(position, indict[x]['x'], yerr=indict[x]['err'], fmt='-o', color=linecolors[colocount], ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='Total '+str(indict[x]['name']))
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
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=1,title=str(linelevel)+' (≥'+str(filterx)+"%) \n at " + str(focaltax) +' level' +'\n(Ascending abundance)', fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend

def plot_my_df_with_line_no_ylim(indf, filterx, toggle, linelevel, focaltax, indict):
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
    linecolors=['black', 'blue', 'green', 'red']
    if toggle =='yes':
        colocount=0
        for x in indict:
            ax.errorbar(position, indict[x]['x'], yerr=indict[x]['err'], fmt='-o', color=linecolors[colocount], ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='Total '+str(indict[x]['name']))
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=3, edgecolor='black', color=colors[x], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_xticks(position)
    ax.set_xticklabels(list(indf.columns))
    sns.despine(top=True, right=True)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Relative abundance (%)", fontsize=24)
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=1,title=str(linelevel)+' (≥'+str(filterx)+"%) \n at " + str(focaltax) +' level' +'\n(Ascending abundance)', fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend

def get_me_relative_stats(indf, mergedotutable, otutable):
    relative=indf.sum(axis=0)/mergedotutable.sum(axis=0)*100
    relative_err=[]
    cages=['A-','B-','C-','D-','E-']
    for number in list(indf.columns): 
        if not number == '0':
            temperr=[]
            for cage in list(otutable.columns):
                for letter in cages:
                    if re.search(str(number)+letter, cage):
                        temperr.append(indf[number].sum(axis=0)/otutable[cage].sum(axis=0)*100)
            relative_err.append((np.mean(np.array(temperr)))**(1/5))
    #append a 0 to error as d0 has no replicates
    relative_err=[0]+relative_err
    return relative, relative_err

def get_mean_get_err(df, tax):
    mean=df.loc[tax].mean()
    err=(df.loc[tax]).sem()
    print('Mean = ' + str(mean)+'\nStErr = ' +str(err))


def plot_my_df_ncol2(indf, filterx):
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
    legend=ax.legend(handles[::-1], labels[::-1], loc="upper left", bbox_to_anchor=(1.15,1.05), ncol=2,title="Taxonomy (≥"+str(filterx)+"%)\n(Ascending abundance)", fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend


#############Barplot_outputs#######################
taxonomic_index_list=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

###
filterx=0.0001
df=tax_extract(merged_otu_table, tax_table, 'Kingdom')
#df=format_18S_index(df)
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)
fig.savefig("Barplots/Kingdom_barplot_high_level2.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')
###

###
filterx=0.6
df=tax_extract(merged_otu_table, tax_table, 'Phylum')
df_filtered_otu_phylum, df_filtered_normalised_phylum=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised_phylum, filterx)
fig.savefig("Barplots/Phylum_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')
###

###
filterx=2
df=tax_extract(merged_otu_table, tax_table, 'Class')
df_filtered_otu_class, df_filtered_normalised_class=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised_class, filterx)
fig.savefig("Barplots/Class_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')
###

###
filterx=3.25
df=tax_extract(merged_otu_table, tax_table, 'Order')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised, filterx)
fig.savefig("Barplots/Order_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')
###

###
filterx=5.25
df=tax_extract(merged_otu_table, tax_table, 'Family')
df_filtered_otu_family, df_filtered_normalised_family=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df(df_filtered_normalised_family, filterx)
fig.savefig("Barplots/Family_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')
###
filterx=0
df=tax_extract(otu_table, tax_table, 'Family')
df_filtered_otu_family_unmerged_unfiltered, df_filtered_normalised_family_unmerged_unfiltered=tax_normalise_filter(df, filterx)
#fig, lgd=plot_my_df(df_filtered_normalised_family_unmerged_unfiltered, filterx)
#fig.savefig("Barplots/Family_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')
###
filterx=1.75
df=tax_extract(merged_otu_table, tax_table, 'Family')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df_ncol2(df_filtered_normalised, filterx)
fig.savefig("Barplots/Family_barplot_high_level_ncol2.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')

taxonomic_index_list=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

#############################Look at Archaea#################################
filterx=0.05
harvest_tax='Archaea'
harvest_tax_level='Kingdom' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Phylum'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Phylum')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}#'a':{'name':harvest_tax, 'x':reltot, 'err':reltoterr}}

fig, lgd=plot_my_df_with_line(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
fig.savefig("Barplots/Archaea_phylum_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')  
#########################################################################################

#############################Look at Archaeal phyla#################################
filterx=0.05
harvest_tax='Archaea'
harvest_tax_level='Kingdom' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Order'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Phylum')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   


linedict={}#'a':{'name':harvest_tax, 'x':reltot, 'err':reltoterr}}

fig, lgd=plot_my_df_with_line(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
fig.savefig("Barplots/Archaea_orders_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all') 
#########################################################################################
#PROTEOMIC CROSS REFERENCING
CAZy_family_list=['NA',
 'Vibrionaceae',
 'Flavobacteriaceae',
 'Cellvibrionaceae',
 'Alteromonadaceae',
 'Cytophagaceae',
 'Saccharospirillaceae',
 'Prolixibacteraceae',
 'Sphingobacteriaceae',
 'Marinilabiliaceae',
 'Bacteroidaceae',
 'Paludibacteraceae',
 'Chitinophagaceae',
 'Polyangiaceae',
 'Lentimicrobiaceae',
 'Flammeovirgaceae',
 'Melioribacteraceae',
 'Desulfovibrionaceae',
 'Peptococcaceae',
 'Desulfobulbaceae',
 'Planococcaceae',
 'Cyclobacteriaceae',
 'Odoribacteraceae',
 'Demequinaceae',
 'Chromatiaceae',
 'Paenibacillaceae',
 'Colwelliaceae',
 'Rubritaleaceae',
 'Pseudomonadaceae',
 'Streptosporangiaceae',
 'Lachnospiraceae',
 'Coccomyxaceae',
 'Oleiphilaceae',
 'Desulfuromonadaceae',
 'Glossiphoniidae',
 'Salinivirgaceae',
 'Oceanospirillaceae',
 'Ruminococcaceae',
 'Halomonadaceae',
 'Acholeplasmataceae',
 'Marinifilaceae',
 'Opitutaceae',
 'Lewinellaceae']

CAZy_order_list=['Vibrionales',
 'Alteromonadales',
 'NA',
 'Flavobacteriales',
 'Cellvibrionales',
 'Cytophagales',
 'Marinilabiliales',
 'Desulfuromonadales',
 'Bacteroidales',
 'Oceanospirillales',
 'Sphingobacteriales',
 'Chitinophagales',
 'Myxococcales',
 'Clostridiales',
 'Bacillales',
 'Ignavibacteriales',
 'Desulfovibrionales',
 'Desulfobacterales',
 'Micrococcales',
 'Chromatiales',
 'Verrucomicrobiales',
 'Pseudomonadales',
 'Streptosporangiales',
 'Rhynchobdellida',
 'Acholeplasmatales',
 'Opitutales',
 'Saprospirales']

###The input Gdff dataframe must be created to the same level, ORDER or FAMILY as the input heatmap dataframe ISDE script lines 2164 onwards

fig, lgd, molpct=plot_my_df_cazy(Gdff_filtered_normalised, 0, df_class_taxa_heatmap, 'Gammaproteobacteria', 'Order', [])
fig.savefig("CAZy_producing_gammaproteobacteria_order_level", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')  

######Filtering microbiome by CAZy#######
#NB/filter family, order level microbiome and use the df from above#
filterx=0
df=tax_extract(merged_otu_table, tax_table, 'Family')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd, dfcazy=plot_my_df_cazy_ncol_fam(df_filtered_normalised, filterx, df_class_taxa_heatmap, 'Family', 'Family', [])
fig.savefig("Microbiome_CAZY_FAMILY_LEVEL.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')  
#
#filterx=0
#df=tax_extract(merged_otu_table, tax_table, 'Order')
#df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
#fig, lgd, dfcazy=plot_my_df_cazy_ncol_class(df_filtered_normalised, filterx, df_class_taxa_heatmap, 'Order', 'Order')
#fig.savefig("Microbiome_CAZY_ORDER_LEVEL.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all')  
##
#filterx=0
#df=tax_extract(merged_otu_table, tax_table, 'Class')
#df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
#fig, lgd, dfcazy=plot_my_df_cazy_ncol_class(df_filtered_normalised, filterx, df_class_taxa_heatmap, 'Class', 'Class')
#fig.savefig("Microbiome_CAZY_CLASS_LEVEL.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all') 

#####Microbiome heatmap family level######
flist=[]
for x in list(df_filtered_normalised.index):
    if x in list(df_class_taxa_heatmap.index):
        flist.append(x)
df_molpct=df_class_taxa_heatmap.loc[flist]
df_molpct.columns=['1','3','5','10']
df_abundance=df_filtered_normalised.loc[flist][['1','3','5','10']]
#df_abundance.columns=df_abundance.columns.astype(int)
df_abundance.columns=['1','3','5','10']
#df_abundance.index.name='Taxonomy'


filterlist=[]#['Peptococcaceae','Salinivirgaceae']#'Cytophagaceae'
df_molpct=df_molpct.drop(filterlist)
df_abundance=df_abundance.drop(filterlist)
df_molpct=df_molpct.replace(0, np.nan)
df_abundance=df_abundance.replace(0, np.nan)
hm=(df_molpct/df_abundance).astype(float)

#hm.columns=['One','Three','Five','Ten']

fig, ax = plt.subplots(figsize=(7,12))
sns.set(font_scale=1.6)
hmp=sns.heatmap(hm.apply(np.log10), annot=True, annot_kws={"size": 16}, cbar_kws={'label': 'Producivity index'}, cmap='RdBu', linewidths=.5)
#fig.savefig("PNAS_plots/Productvity_index_heatmap_SILVA.png",bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')

#############################Look at Gammaproteobacteria#################################
filterx=1
harvest_tax='Gammaproteobacteria'
harvest_tax_level='Class' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Order'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Phylum')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

highlevelplot='Proteobacteria'
highlevelx=tax_extract(merged_otu_table, tax_table, 'Phylum').loc[highlevelplot]/tax_extract(merged_otu_table, tax_table, 'Phylum').sum(axis=0)*100
a, highlevelxerr=get_me_relative_stats(tax_extract(merged_otu_table, tax_table, 'Phylum'), merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={'a':{'name':harvest_tax, 'x':reltot, 'err':reltoterr}}

fig, lgd=plot_my_df_with_line(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
fig.savefig("Barplots/Gammaproteobacteria_family_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all') 

##########################################################################################
#filterx=5
#harvest_tax='Gammaproteobacteria'
#harvest_tax_level='Class' # tax level of harvest tax
##focal_taxonomy is the taxonomic level within harvest tax you want to look at
#focal_taxonomy='Family'
#Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
#GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Phylum')
#reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)
#
#highlevelplot='Proteobacteria'
#highlevelx=tax_extract(merged_otu_table, tax_table, 'Phylum').loc[highlevelplot]/tax_extract(merged_otu_table, tax_table, 'Phylum').sum(axis=0)*100
#a, highlevelxerr=get_me_relative_stats(tax_extract(merged_otu_table, tax_table, 'Phylum'), merged_otu_table, otu_table)
#
#Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
#a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   
#
#linedict={'a':{'name':harvest_tax, 'x':reltot, 'err':reltoterr}}
#
#fig, lgd=plot_my_df_with_line(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
#fig.savefig("Barplots/Gammaproteobacteria_family_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all') 

#############################################################################################


##Look at Bacteroidetes
#filterx=0
#harvest_tax='Bacteroidetes'
#harvest_tax_level='Phylum' # tax level of harvest tax
##focal_taxonomy is the taxonomic level within harvest tax you want to look at
#focal_taxonomy='Order'
#Bdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
#BDFtotal=tax_extract(merged_otu_table.loc[Bdf.index], Bdf, 'Phylum')
#reltot, reltoterr = get_me_relative_stats(BDFtotal, merged_otu_table, otu_table)
#
#highlevelplot='Bacteroidetes'
#highlevelx=tax_extract(merged_otu_table, tax_table, 'Phylum').loc[highlevelplot]/tax_extract(merged_otu_table, tax_table, 'Phylum').sum(axis=0)*100
#a, highlevelxerr=get_me_relative_stats(tax_extract(merged_otu_table, tax_table, 'Phylum'), merged_otu_table, otu_table)
#
#Bdff=tax_extract(merged_otu_table.loc[Bdf.index], Bdf, focal_taxonomy)   
#a, Bdff_filtered_normalised=tax_normalise_filter(Bdff, filterx)
#
#linedict={'a':{'name':harvest_tax, 'x':reltot, 'err':reltoterr}}
#
#fig, lgd=plot_my_df_with_line(Bdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
#fig.savefig("Barplots/Bacteroidetes_order_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all')
#
#fig, lgd, molpct=plot_my_df_cazy(Bdff_filtered_normalised, 0, df_class_taxa_heatmap, 'Bacteroidetes', 'Family', ['Marinilabiaceae'])
#fig.savefig("CAZy_producing_bacteroidetes_Family_level", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all')  

#############################################################################
#
#filterx=1
#harvest_tax='Bacteroidetes'
#harvest_tax_level='Phylum' # tax level of harvest tax
##focal_taxonomy is the taxonomic level within harvest tax you want to look at
#focal_taxonomy='Family'
#Bdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
#BDFtotal=tax_extract(merged_otu_table.loc[Bdf.index], Bdf, 'Phylum')
#reltot, reltoterr = get_me_relative_stats(BDFtotal, merged_otu_table, otu_table)
#
#highlevelplot='Bacteroidetes'
#highlevelx=tax_extract(merged_otu_table, tax_table, 'Phylum').loc[highlevelplot]/tax_extract(merged_otu_table, tax_table, 'Phylum').sum(axis=0)*100
#a, highlevelxerr=get_me_relative_stats(tax_extract(merged_otu_table, tax_table, 'Phylum'), merged_otu_table, otu_table)
#
#Bdff=tax_extract(merged_otu_table.loc[Bdf.index], Bdf, focal_taxonomy)   
#a, Bdff_filtered_normalised=tax_normalise_filter(Bdff, filterx)
#
#linedict={'a':{'name':harvest_tax, 'x':reltot, 'err':reltoterr}}

#fig, lgd=plot_my_df_with_line(Bdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
#fig.savefig("Barplots/Bacteroidetes_family_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all')  

#########################DELTYAPROTEOBACTERIA#####################################
#filterx=0.5
#harvest_tax='Deltaproteobacteria'
#harvest_tax_level='Class' # tax level of harvest tax
##focal_taxonomy is the taxonomic level within harvest tax you want to look at
#focal_taxonomy='Order'
#Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
#GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Phylum')
#reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)
#
#highlevelplot='Proteobacteria'
#highlevelx=tax_extract(merged_otu_table, tax_table, 'Phylum').loc[highlevelplot]/tax_extract(merged_otu_table, tax_table, 'Phylum').sum(axis=0)*100
#a, highlevelxerr=get_me_relative_stats(tax_extract(merged_otu_table, tax_table, 'Phylum'), merged_otu_table, otu_table)
#
#Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
#a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   
#
#linedict={'a':{'name':harvest_tax, 'x':reltot, 'err':reltoterr}}
#
#fig, lgd=plot_my_df_with_line(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
#fig.savefig("Barplots/Deltaproteobacteria_order_barplot_high_level.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all') 

######plot microbiome families, at class level...###############################
###############################################################################
###ensure df_class_taxa_heatmap has been resolved to the correct level #####
    
#familyindex=list(df_class_taxa_heatmap.index)
#familyindex.remove('NA')
#familyindex.append('Marinilabiaceae')
#filteredfamilyindexes=list(tax_table[tax_table['Family'].isin(familyindex)].index)
#
##normalise_otu_table_before_extracting
#a, otu_table_norm=tax_normalise_filter(merged_otu_table,0)
#
#otudf=otu_table_norm.loc[filteredfamilyindexes]
#
#newtt=tax_table.loc[filteredfamilyindexes]
#
#taxdfnew=tax_extract(otudf, newtt, 'Class')
#
totalmolpct=np.array([1.176125,1.248065,1.013474,1.209536])
varymolpct=np.array([0.617890,0.542846,0.600238,0.830602])
NApct=np.array([0.227886,  0.252550,  0.167165,  0.186258])
Cellvibrionaceaepct=np.array([0.062594,  0.158285,  0.068753,  0.079686])
var_cel=np.array(varymolpct)+np.array(Cellvibrionaceaepct)
var_cel_Na=var_cel+np.array(NApct)
#
#filterx=0
#fig, lgd=plot_my_dfx(taxdfnew, filterx, varymolpct, totalmolpct, 'Class', 'Families', var_cel, var_cel_Na)
#fig.savefig("MICROBIOME_CAZY_families_displayed_at_CLASS_level_dual_axis.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)
#plt.show()
#plt.clf()
#plt.cla()
#plt.close(fig)  
#plt.close('all') 

###############################################################################################
###############################################################################################
###############################################################################################




######mICROBIOME LEVEL FAMILY 4 PLOTS#####

filterx=0
df=tax_extract(merged_otu_table, tax_table, 'Family')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)
fig, lgd, dfcazy, noncazycol=plot_my_df_cazy_ncol2(df_filtered_normalised, filterx, df_class_taxa_heatmap, 'Family', 'Family', ['Marinilabiaceae'], varymolpct, totalmolpct, 'Class', 'Families', var_cel, var_cel_Na)
fig.savefig("Microbiome_CAZY_FAMILY_LEVEL.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=1200)

def plot_my_dfx(indf, filterx, invaried, intotal, level, prelim, varcel, varcelna):
    colors=color_jenga(indf)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(10,8))
    ax2=ax.twinx()
    ax2.errorbar(np.array([1,3,5,8]), np.array(intotal), fmt='-o', color='blue', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='CAZy total')
    ax2.errorbar(np.array([1,3,5,8]), np.array(varcelna), fmt='-o', color='red', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='Visible + Cellvibrionacea + NA')
    ax2.errorbar(np.array([1,3,5,8]), np.array(varcel), fmt='-o', color='grey', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='Visible + Cellvibrionacea')
    ax2.errorbar(np.array([1,3,5,8]), np.array(invaried), fmt='-o', color='black', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, label='CAZy contributions (visible taxa)')
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ##Remove 'NA' column##
    indf=indf.loc[~(indf.index=='NA')]
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=colors[x], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_xticklabels(list(indf.columns))
    sns.set(style="white")
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)
    ax.set_ylabel("Relative abundance (%)", fontsize=24)
    ax2.tick_params(direction='out', length=6, width=4)  
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_ylabel(u"\u2211" r'$\bar{x}$'" Molar percentage", fontsize=24)
    ax2.set_ylim([0,1.3])
    ax.set_xlabel('Week', fontsize=24)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    handles=h2+handles
    labels=l2+labels
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.125,1.025), ncol=1,title="CAZy producing "+str(prelim)+"\n      At "+str(level)+" level\n(Ascending abundance)", fontsize=20)
    legend.get_title().set_fontsize(20)
    return fig, legend



###########################     DOTPLOT     #############################################

def dotplot(formattedf, molarpctdf):
    molarpctdf.columns=[1,3,5,10]
    grey=(0,0,0,0.2)
    xtix=list(formattedf.columns)
    ytix=list(formattedf.index)
    xvals=list(range(0,len(formattedf.columns)))
    #reverse index to plot largest to smallest
    yvals=list(range(0,len(formattedf.index)))
    fig, ax = plt.subplots(figsize=(10,20))
    yaxiscount=max(yvals)
    multiplier=75
    molarmultiplier=multiplier*90
    alpha=0.4
    cols=sns.hls_palette(max(yvals)+1, l=.5, s=1, h=.5)[::-1]
    #add alpha
    cols=[list(cols[x])+[alpha] for x in list(range(0,len(cols)))]
    proteomelist=[(1,1),(3,3),(5,5),(10,8)]
    #Can also be done using plot and markersize=n
    for index in formattedf.index[::-1]:
        for x in xvals:
            ax.scatter(x, yaxiscount, s=(formattedf.loc[index][x])*multiplier, color=cols[yaxiscount])
        if ytix[yaxiscount] in list(molarpctdf.index):
            for x in proteomelist:
                 ax.scatter(x[1], yaxiscount, s=molarpctdf.loc[ytix[yaxiscount]][x[0]]*molarmultiplier, color=(1,1,1,0), edgecolor='black', linewidths=1)
        yaxiscount+=-1
    ax.set_xticks(xvals)
    ax.set_yticks(yvals)
    ax.set_yticklabels(ytix)
    ax.set_xlabel("Week", fontsize=24)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=20)   
    ax.set_xticklabels(xtix, fontsize=24)
    ax.yaxis.set_label_position("right")
    legendlabels=[0.1,'', 0.25, '',0.5,'', 1,'', 5, '',10, '',15, '', 20, '', 25]
    collabels=sns.hls_palette(len(legendlabels), l=.5, s=.9, h=.5)[::-1]
    collabels=[list(collabels[x])+[alpha] for x in list(range(0,len(collabels)))]

    #if use .plot and ms=n, index the plot e.g. [plt.plot(x, y, z)[0] for x in y]
    #patches=[plt.scatter([],[], s=(x*multiplier), color='blue', label=str(x)+'%') for x in legendlabels]
    patches=[plt.scatter([],[], s=1, color='white', label='white') if legendlabels[x] == '' else plt.scatter([],[], s=legendlabels[x]*multiplier, color=collabels[x], label=str(legendlabels[x])+'%') for x in list(range(0,len(legendlabels)))]
    legend=plt.legend(handles=patches, ncol=1,loc='upper left', bbox_to_anchor=(1.04,1), borderaxespad=0, fontsize=20, title="    Relative\n  abundance")
    legend.get_title().set_fontsize(20)
    whitelabels=[1,3,5,7,9, 11, 13, 15]
    for x in whitelabels:
        plt.setp(legend.get_texts()[x], color='w')
    molarlegendlabels=[0.01, '', 0.05, '', 0.1, '', 0.2, '', 0.3, '', 0.4, '', 0.5]
    patches2=[plt.scatter([],[], s=1, color='white', label='white') if x == '' else plt.scatter([],[], s=(x*molarmultiplier), color=(1,1,1,0), edgecolor='black', linewidths=1, label=str(x)+'%') for x in molarlegendlabels]
    legend2=plt.legend(handles=patches2, ncol=1, bbox_to_anchor=(1.04,0), loc="lower left", borderaxespad=0, fontsize=20, title=u"\u2211" r'$\bar{x}$'" CAZy molar\n    percentage")
    legend2.get_title().set_fontsize(20)
    for x in whitelabels[:-2]:
        plt.setp(legend2.get_texts()[x], color='w')
    ax.add_artist(legend)
    ax.add_artist(legend2)
    ###Third invisible legend...
    leg = plt.legend(handles=[], bbox_to_anchor=(1.4, 0))
    leg.get_frame().set_alpha(0) # hide legend
    return fig

###########################################   PAPER  DOTPLOT     #############################################
###########################################   PAPER  DOTPLOT     #############################################
###########################################   PAPER  DOTPLOT     #############################################
###########################################   PAPER  DOTPLOT     #############################################
##############                                                                              ##################
##############                           START SCRIPT FROM HERE                             ##################
##############                                                                              ##################
##############                                                                              ##################
   
#filter the df to cazy producers only

filterx=0
df=tax_extract(merged_otu_table, tax_table, 'Family')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, filterx)

filterx=0
df_genus=tax_extract(merged_otu_table, tax_table, 'Genus')
df_filtered_otu_genus, df_filtered_normalised_genus=tax_normalise_filter(df_genus, filterx)

##group uncultured/unknown into NA
df_filtered_normalised.loc['NA']=df_filtered_normalised.loc['NA']+df_filtered_normalised.loc['uncultured']+df_filtered_normalised.loc['uncultured bacterium']+df_filtered_normalised.loc['Unknown Family']
df_filtered_normalised=df_filtered_normalised.drop(['uncultured','uncultured bacterium', 'Unknown Family'], axis='rows')
df_filtered_normalised.loc['Clostridiales vadinBB60']=df_filtered_normalised.loc['Clostridiales vadinBB60 group']
df_filtered_normalised=df_filtered_normalised.drop(['Clostridiales vadinBB60 group'])

fig, lgd, dfcazy, noncazycol=plot_my_df_cazy_ncol2(df_filtered_normalised, filterx, df_class_taxa_heatmap, 'Family', 'Family', [], varymolpct, totalmolpct, 'Class', 'Families', var_cel, var_cel_Na)

fig_genus, lgd_genus, dfcazy_genus, noncazycol_genus=plot_my_df_cazy_ncol2(df_filtered_normalised_genus, filterx, df_genus_taxa_heatmap, 'Family', 'Family', [], varymolpct, totalmolpct, 'Class', 'Families', var_cel, var_cel_Na)

plt.show()
plt.clf()
plt.cla()
plt.close(fig)  
plt.close('all')

##swamp marini labels as they differ between dataframes, to Marinilabiaceae
#if 'Marinilabiliaceae' in df_class_taxa_heatmap.index:
#     print('changing Marinilabiliaceae to Marinilabiaceae')
#     df_class_taxa_heatmap.loc['Marinilabiaceae'] = df_class_taxa_heatmap.loc['Marinilabiliaceae']
#     df_class_taxa_heatmap=df_class_taxa_heatmap.drop('Marinilabiliaceae', axis=0)


#plot dotplot using dfcazy
#NB ensure they are same taxonomic resoltuion
#[::-1] to reverse order of plots
fig =dotplot(dfcazy.iloc[::-1], df_class_taxa_heatmap)
fig.savefig("Microbiome_DOTPLOT_family_level_resolution_1200DPI.png", bbox_inches='tight', dpi=1200)

noncazyproducers_families=[]
for x in df_class_taxa_heatmap.index:
    if x not in df_molpct.index:
        noncazyproducers_families.append(x)

noncazyproducers_genus=[]
for x in df_genus_taxa_heatmap.index:
    if x not in list(tax_table['Genus']):
        noncazyproducers_genus.append(x)

noncazydf_rowcol={'Gammaproteobacteria':['Oceanospirillaceae'],\
        'Flavobacteriia':['Cytophagaceae'],\
        'Bacteroidia':['Odoribacteraceae'],\
 'Deltaproteobacteria':['Polyangiaceae'],\
 'Actinobacteria':['Streptosporangiaceae'],\
 'Chlorophyceae':['Coccomyxaceae'],\
 'Clitellata':['Glossiphoniidae'],\
 'Saprospirae':['Lewinellaceae']}

noncazydf_rowcol_genus={}
for x in noncazyproducers_genus:
    if not x == 'NA':
         rowcollevel='class'
         for orf in taxa_master_dict:
             if taxa_master_dict[orf]['Taxonomy']['genus']==x:
                 result=taxa_master_dict[orf]['Taxonomy'][rowcollevel]
                 if result in noncazydf_rowcol_genus:
                     if x not in noncazydf_rowcol_genus[result]:
                            noncazydf_rowcol_genus[result].append(x)
                 if result not in noncazydf_rowcol_genus:
                       noncazydf_rowcol_genus[result]=[]
                       noncazydf_rowcol_genus[result].append(x)
                 continue
    
new_SILVA_fams=['Cellvibrionaceae','Prolixibacteraceae','Paludibacteraceae','Lentimicrobiaceae',\
                'Marinifilaceae','Melioribacteraceae','Demequinaceae','Rubritaleaceae']

def extract_x_level_for_col_color(indf, taxdf, inlevel, outlevel):
    taxlist=[]
    for x in list(indf.index):
        if not x == 'NA':
             tax=str(taxdf.loc[taxdf[taxdf[inlevel]==x].index[0]][outlevel])
             if tax not in taxlist:
                 taxlist.append(tax)
    #make color dict
    cols=sns.hls_palette(int((len(taxlist)+1)/2), l=.45, s=.7, h=1)
    alpha=0.4
    cols=[list(cols[x])+[alpha] for x in list(range(0,len(cols)))]
    cols2=sns.hls_palette(int((len(taxlist)+1)/2), l=.65, s=1, h=.5)
    cols=cols2+cols
    cols=cols[:len(taxlist)]
    coldict=dict(zip(taxlist, cols))
    #re-loop to construct dictionary with colors
    taxdict={}
    for x in list(indf.index):
        if not x == 'NA':
             tax=str(taxdf.loc[taxdf[taxdf[inlevel]==x].index[0]][outlevel])
             if tax not in taxdict:
                 taxdict[tax]={}
             if x not in taxdict[tax]:
                 taxdict[tax][x]={}
                 taxdict[tax][x]['colour']=coldict[tax]
    #taxdict['NA']={}
    #taxdict['NA']['NA']={}
    #taxdict['NA']['NA']['colour']=[1,1,1]
    return taxdict

def extract_x_level_for_col_color_global(indf, taxdf, inlevel, outlevel, noncazy):
    taxlist=[]
    for x in list(indf.index):
        if not x == 'NA':
             tax=str(taxdf.loc[taxdf[taxdf[inlevel]==x].index[0]][outlevel])
             if tax not in taxlist:
                 taxlist.append(tax)
    #make color dict
    for x in noncazy.keys():
        if x not in taxlist:
            taxlist.append(x)
    len_total_classes=len(set(list(indf.index)+list(noncazy.keys())))
    cols=sns.hls_palette(int((len_total_classes+1)/2), l=.45, s=.7, h=1)
    alpha=0.4
    cols=[list(cols[x])+[alpha] for x in list(range(0,len(cols)))]
    cols2=sns.hls_palette(int((len(taxlist)+1)/2), l=.65, s=1, h=.5)
    cols=cols2+cols
    cols=cols[:len(taxlist)]
    coldict=dict(zip(taxlist, cols))
    #re-loop to construct dictionary with colors
    taxdict={}
    for x in list(indf.index):
        if not x == 'NA':
             tax=str(taxdf.loc[taxdf[taxdf[inlevel]==x].index[0]][outlevel])
             if tax not in taxdict:
                 taxdict[tax]={}
             if x not in taxdict[tax]:
                 taxdict[tax][x]={}
                 taxdict[tax][x]['colour']=coldict[tax]
    taxdict['NA']={}
    taxdict['NA']['NA']={}
    taxdict['NA']['NA']['colour']=[0,0,0]
    for x in noncazy:
        if x in taxdict:
            for y in noncazy[x]:
                taxdict[x][y]={}
                taxdict[x][y]['colour']=coldict[x]
        if x not in taxdict:
            taxdict[x]={}
            for y in noncazy[x]:
                taxdict[x][y]={}
                taxdict[x][y]['colour']=coldict[x] 
    return taxdict

row_color_dict=extract_x_level_for_col_color(dfcazy, tax_table, 'Family', 'Class')
row_color_dict_genus=extract_x_level_for_col_color(dfcazy_genus, tax_table, 'Genus', 'Family')


row_color_dict_global=extract_x_level_for_col_color_global(dfcazy, tax_table, 'Family', 'Class', noncazydf_rowcol)
row_color_dict_global_genus=extract_x_level_for_col_color_global(dfcazy_genus, tax_table, 'Genus', 'Class', noncazydf_rowcol_genus)
row_color_dict_global_genus_family=extract_x_level_for_col_color_global(dfcazy_genus, tax_table, 'Genus', 'Family', noncazydf_rowcol_genus)

def dotplot_small(formattedf, molarpctdf, rowcoldict1):
    fontsizevar=14
    molarpctdf.columns=[1,3,5,10]
    grey=(0,0,0,0.2)
    xtix=list(formattedf.columns)
    ytix=list(formattedf.index)
    xvals=list(range(0,len(formattedf.columns)))
    #reverse index to plot largest to smallest
    yvals=list(range(0,len(formattedf.index)))
    fig, ax = plt.subplots(figsize=(3,10))
    yaxiscount=max(yvals)
    multiplier=17
    molarmultiplier=multiplier*50
    alpha=0.5
    #cols=sns.hls_palette(max(yvals)+1, l=.5, s=1, h=.5)[::-1]
    cols=sns.xkcd_palette(['clear blue'])*(max(yvals)+1)
    #add alpha
    cols=[list(cols[x])+[alpha] for x in list(range(0,len(cols)))]
    proteomelist=[(1,1),(3,3),(5,5),(10,8)]
    #Can also be done using plot and markersize=n
    for index in formattedf.index[::-1]:
        for x in xvals:
            ax.scatter(x, yaxiscount, s=(formattedf.loc[index][x])*multiplier, color=cols[yaxiscount])
        if ytix[yaxiscount] in list(molarpctdf.index):
            for x in proteomelist:
                 ax.scatter(x[1], yaxiscount, s=molarpctdf.loc[ytix[yaxiscount]][x[0]]*molarmultiplier, color=(1,1,1,0), edgecolor='black', linewidths=1)
                 #assign y axis value to row col dictionary, which can be toggled on off at the bottom of the script
                 for x in rowcoldict1:
                      if index in rowcoldict1[x]:
                          if 'yposition' not in rowcoldict1[x][index]:
                              rowcoldict1[x][index]['yposition']=yaxiscount
        yaxiscount+=-1
   
    ax.set_xticks(xvals)
    ax.set_yticks(yvals)
    ax.set_yticklabels(ytix)
    ax.set_xlabel("Week", fontsize=18)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=fontsizevar)   
    ax.set_xticklabels(xtix, fontsize=14)
    ax.yaxis.set_label_position("right")
    #ax.set_ylim([0.5, yaxiscount+0.5])
    legendlabels=[0.25,0.5, 1, 5, 10, 25]
    #collabels=sns.hls_palette(len(legendlabels), l=.5, s=.9, h=.5)[::-1]
    collabels=sns.xkcd_palette(['clear blue'])*(len(legendlabels))
    collabels=[list(collabels[x])+[alpha] for x in list(range(0,len(collabels)))]
    #if use .plot and ms=n, index the plot e.g. [plt.plot(x, y, z)[0] for x in y]
    #patches=[plt.scatter([],[], s=(x*multiplier), color='blue', label=str(x)+'%') for x in legendlabels]
    patches=[plt.scatter([],[], s=1, color='white', label='white') if legendlabels[x] == '' else plt.scatter([],[], s=legendlabels[x]*multiplier, color=collabels[x], label=str(legendlabels[x])+'%') for x in list(range(0,len(legendlabels)))]
    legend=plt.legend(handles=patches, ncol=1,loc='upper left', bbox_to_anchor=(1.0,1), borderaxespad=0, fontsize=fontsizevar, title="Abundance", frameon=False)
    legend.get_title().set_fontsize(fontsizevar)
    molarlegendlabels=[0.05, 0.1, 0.2, 0.3,  0.4,0.5]
    patches2=[plt.scatter([],[], s=1, color='white', label='white') if x == '' else plt.scatter([],[], s=(x*molarmultiplier), color=(1,1,1,0), edgecolor='black', linewidths=1, label=str(x)+'%') for x in molarlegendlabels]
    legend2=plt.legend(handles=patches2, ncol=1, bbox_to_anchor=(1.0,0.4), loc="lower left", borderaxespad=0, fontsize=fontsizevar, title=u"\u2211" r'$\bar{x}$'" CAZy mol%", frameon=False)
    legend2.get_title().set_fontsize(fontsizevar)
    ax.add_artist(legend)
    ax.add_artist(legend2)
    ###Third invisible legend...
    leg = plt.legend(handles=[], bbox_to_anchor=(1.4, 0))
    leg.get_frame().set_alpha(0) # hide legend
    ##############COLUMN COLOUR HEATPLOT STUFF#######
    offset=0.5
    patches3=[]
    patchlist=[]
    changedict={'Gammaproteobacteria':'γ-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria'}
    #Add row_colour patches
    for x in rowcoldict1:
        for y in rowcoldict1[x]:
            ax.add_patch(plt.Rectangle((-1.175,rowcoldict1[x][y]['yposition']-offset),0.25, 0.85,facecolor=rowcoldict1[x][y]['colour'],clip_on=False,linewidth = 0))
            #if statement to stop duplicate legend handles
            if x not in patchlist:
                if x in changedict:
                     patchlist.append(x)
                     patches3.append(mpatches.Patch(color=rowcoldict1[x][y]['colour'], label=str(changedict[x]).strip('[]')))
                else:
                     patchlist.append(x)
                     patches3.append(mpatches.Patch(color=rowcoldict1[x][y]['colour'], label=str(x).strip('[]')))
    legend3=plt.legend(handles=patches3, ncol=1,bbox_to_anchor=(1.0,-0.075), loc="lower left", borderaxespad=0, fontsize=fontsizevar-2, title='Class', frameon=False, handlelength=0.35, handleheight=1.5)
    legend3.get_title().set_fontsize(fontsizevar)
    return fig, rowcoldict1

fig, rcd =dotplot_small(dfcazy.iloc[::-1], df_class_taxa_heatmap, row_color_dict)
fig_genus, rcd_genus =dotplot_small(dfcazy_genus.iloc[::-1], df_genus_taxa_heatmap, row_color_dict_genus)
#fig.savefig("DOTPLOT_updated_5_10.png", dpi=1200, bbox_inches='tight')

def give_me_Xleveltaxa_that_are_cazy_producing_from_X_level(taxa, inlevel, tax_table, Xleveltaxa, cazydf):
    result=set(tax_table[tax_table[inlevel]==taxa][Xleveltaxa])
    print('All taxa recovered = ' + str(result))
    cazyresult=[x for x in result if x in cazydf.index]
    print('CAZy producing only = ' + str(cazyresult))
give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Flammeovirgaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)



#####################################LINKED LEGEND PLOTS#######################################
###########PAPER VERSION####################

######plot microbiome families, at class level...###############################
###############################################################################
###ensure df_class_taxa_heatmap has been resolved to the correct level #####
    

###extract cazy producing taxa to separate dataframe###

familyindex=list(df_class_taxa_heatmap.index)
familyindex.remove('NA')
familyindex.append('Marinilabiaceae')
filteredfamilyindexes=list(tax_table[tax_table['Family'].isin(familyindex)].index)

#normalise_otu_table_before_extracting
a, otu_table_norm=tax_normalise_filter(merged_otu_table,0)

otudf=otu_table_norm.loc[filteredfamilyindexes]

newtt=tax_table.loc[filteredfamilyindexes]

#Change to class for class...
taxdfnew=tax_extract(otudf, newtt, 'Family')
####UNMERGED OTU TABLE FOR STATS FOR PAPER#####
a1, otu_table_norm_UNMERGED=tax_normalise_filter(otu_table,0)
otudf_UNMERGED=otu_table_norm_UNMERGED.loc[filteredfamilyindexes]
taxdfnew_UNMERGED=tax_extract(otudf_UNMERGED, newtt, 'Family')
#taxdfnew_UNMERGED[['1A-16S','1B-16S','1C-16S','1D-16S','1E-16S']].sum(axis=0).std()
#(taxdfnew_UNMERGED[['1A-16S','1B-16S','1C-16S','1D-16S','1E-16S']].sum(axis=0).reset_index(drop=True)/float(taxdfnew_UNMERGED[['D0-16S']].sum(axis=0))).std()
###Raw plot at X level resolution ###

abundance_filterx=1.3
df=tax_extract(merged_otu_table, tax_table, 'Family')
df_filtered_otu, df_filtered_normalised=tax_normalise_filter(df, abundance_filterx)
df_filtered_normalised=df_filtered_normalised.fillna(0).astype(float)
df_filtered_normalised.loc['NA']=df_filtered_normalised.loc['NA']+df_filtered_normalised.loc['uncultured']+df_filtered_normalised.loc['uncultured bacterium']+df_filtered_normalised.loc['Unknown Family']
df_filtered_normalised=df_filtered_normalised.drop(['uncultured','uncultured bacterium', 'Unknown Family'], axis='rows')
df_filtered_normalised.loc['Clostridiales vadinBB60']=df_filtered_normalised.loc['Clostridiales vadinBB60 group']
df_filtered_normalised=df_filtered_normalised.drop(['Clostridiales vadinBB60 group'])

##############################GENUS###################################
familyindex_genus=list(df_genus_taxa_heatmap.index)
familyindex_genus.remove('NA')
filteredfamilyindexes_genus=list(tax_table[tax_table['Genus'].isin(familyindex_genus)].index)
#normalise_otu_table_before_extracting
a_genus, otu_table_norm_genus=tax_normalise_filter(merged_otu_table,0)
otudf_genus=otu_table_norm.loc[filteredfamilyindexes_genus]
newtt_genus=tax_table.loc[filteredfamilyindexes_genus]
#Change to class for class...
taxdfnew_genus=tax_extract(otudf_genus, newtt_genus, 'Genus')
####UNMERGED OTU TABLE FOR STATS FOR PAPER#####
a1_genus, otu_table_norm_UNMERGED_genus=tax_normalise_filter(otu_table,0)
otudf_UNMERGED_genus=otu_table_norm_UNMERGED.loc[filteredfamilyindexes_genus]
taxdfnew_UNMERGED_genus=tax_extract(otudf_UNMERGED_genus, newtt_genus, 'Genus')
abundance_filterx_genus=0.5
df_genus=tax_extract(merged_otu_table, tax_table, 'Genus')
df_filtered_otu_genus, df_filtered_normalised_genus=tax_normalise_filter(df_genus, abundance_filterx_genus)
df_filtered_normalised_genus=df_filtered_normalised_genus.fillna(0).astype(float)
##############################GENUS###################################

##############################CLASS###################################
familyindex_class=list(df_classx_taxa_heatmap.index)
familyindex_class.remove('NA')
filteredfamilyindexes_class=list(tax_table[tax_table['Class'].isin(familyindex_class)].index)
#normalise_otu_table_before_extracting
a_class, otu_table_norm_class=tax_normalise_filter(merged_otu_table,0)
otudf_class=otu_table_norm.loc[filteredfamilyindexes_class]
newtt_class=tax_table.loc[filteredfamilyindexes_class]
#Change to class for class...
taxdfnew_class=tax_extract(otudf_class, newtt_class, 'Class')
####UNMERGED OTU TABLE FOR STATS FOR PAPER#####
a1_class, otu_table_norm_UNMERGED__class=tax_normalise_filter(otu_table,0)
otudf_UNMERGED_class=otu_table_norm_UNMERGED.loc[filteredfamilyindexes_class]
taxdfnew_UNMERGED_class=tax_extract(otudf_UNMERGED_class, newtt_class, 'Class')
abundance_filterx_class=0.1
df_class=tax_extract(merged_otu_table, tax_table, 'Class')
df_filtered_otu_class, df_filtered_normalised_class=tax_normalise_filter(df_class, abundance_filterx_class)
df_filtered_normalised_class=df_filtered_normalised_class.fillna(0).astype(float)
##################################CLASS###################################

#Global colour dict should be performed on the UNFILTERED taxa dictionary to link colors to taxa permanently across graphs
def global_color_dict(indf, cazyproducingtaxadf, modified_df):
    constant_palette=sns.xkcd_palette(["sky blue", "strawberry", "cerise", "light pink", "light yellow", "light green", 'orchid', 'black', 'clear blue','white', 'light brown', 'light grey', 'bright orange', 'pink/purple', 'off white'])
    coldict={'Flavobacteriaceae': constant_palette[8], 'Desulfobulbaceae': constant_palette[0],'Marinilabiliaceae': constant_palette[2], 'Saccharospirillaceae': constant_palette[3], 'Alteromonadaceae': constant_palette[4], 'Vibrionaceae': constant_palette[5], 'Cyclobacteriaceae':constant_palette[13], 'NA':constant_palette[7], 'Prolixibacteraceae':constant_palette[6], 'Cellvibrionaceae':constant_palette[9], 'Rhodobacteraceae':constant_palette[10], 'Spirochaetaceae':constant_palette[11], 'Demequinaceae':constant_palette[12], 'Taxa < 1.3%':constant_palette[14]}
    resids=sns.xkcd_palette(['grey', 'light grey', 'slate', 'bright orange', 'fluorescent green', 'bright pink', 'neon blue', 'bright red', 'desert'])
    uniques=set(list(indf.index)+list(cazyproducingtaxadf.index)+list(modified_df.index))
    length=int((len(uniques)-len(coldict)-len(resids))/2)+2
    #plus 1 so as to not get caught out when int round 0.45 and 0.49 to 0...
    lightpal=sns.hls_palette(length, l=.75, s=.9, h=0.1)[::-1]
    darkpal=sns.hls_palette(length, l=.5, s=.7, h=.7)
    newpal=list(resids)
    for x in list(range(0, length)):
        newpal.append(lightpal[x])
        newpal.append(darkpal[x])
    indexcount=0
    for x in uniques:
        if x not in coldict:
            coldict[x]=newpal[indexcount]
            indexcount+=1
    return coldict

#Global colour dict should be performed on the UNFILTERED taxa dictionary to link colors to taxa permanently across graphs
def global_color_dict_genus(indf, cazyproducingtaxadf, modified_df):
    constant_palette=sns.xkcd_palette(["sky blue", "strawberry", "cerise", "light pink", "light yellow", "light green", 'orchid', 'black', 'clear blue','white', 'light brown', 'light grey', 'bright orange', 'pink/purple', 'off white'])
    coldict={}
    resids=sns.xkcd_palette(['grey', 'light grey', 'slate', 'bright orange', 'fluorescent green', 'bright pink', 'neon blue', 'bright red', 'desert'])
    uniques=set(list(indf.index)+list(cazyproducingtaxadf.index)+list(modified_df.index))
    length=int((len(uniques)-len(resids))/2)+2
    #plus 1 so as to not get caught out when int round 0.45 and 0.49 to 0...
    lightpal=sns.hls_palette(length, l=.75, s=.9, h=0.1)[::-1]
    darkpal=sns.hls_palette(length, l=.5, s=.7, h=.7)
    newpal=list(resids)
    for x in list(range(0, length)):
        newpal.append(lightpal[x])
        newpal.append(darkpal[x])
    indexcount=0
    for x in uniques:
        if x not in coldict:
            coldict[x]=newpal[indexcount]
            indexcount+=1
    return coldict
#########Use raw data to create linked colors to taxa so that graphs are consistent and can share a legend ##########

x_level_tax_col_dict=global_color_dict(df_filtered_otu, taxdfnew, df_filtered_normalised)
x_level_tax_col_dict_class=global_color_dict_genus(df_filtered_otu_class, taxdfnew_class, df_filtered_normalised_class)
x_level_tax_col_dict_genus=global_color_dict_genus(df_filtered_otu_genus, taxdfnew_genus, df_filtered_normalised_genus)
###update definition to include newly constructed dictionary and select colours based on this

def plot_my_df_ncol3(indf, filterx, coldict):
    colors=color_jenga(indf)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(5,6))
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=coldict[indf.iloc[x].name], label=indf.iloc[x].name)
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
   # legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=2,title="Taxonomy (≥"+str(filterx)+"%)\n(Ascending abundance)", fontsize=20)
  #  legend.get_title().set_fontsize(20)
    plt.tight_layout()
    return fig#, legend

fig =plot_my_df_ncol3(df_filtered_normalised, filterx, x_level_tax_col_dict)
fig_genus=plot_my_df_ncol3(df_filtered_normalised_genus, filterx, x_level_tax_col_dict_genus)
#fig.savefig("LINKED_family_microbiome_plot.png", dpi=1200, bbox_to_inches="tight")


def plot_my_dfx2(indf, filterx, level, prelim, otudata, coldict):
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
    indf=indf.loc[~(indf.index=='NA')]
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        #if indf.iloc[x].name in coldict:
             ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=coldict[indf.iloc[x].name], label=indf.iloc[x].name)
             yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
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
    ax2.set_ylim([0,2500])
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
    legend=ax.legend(handles[::-1], labels[::-1], loc=(1.2, 0.5), ncol=1,title="CAZy producing "+str(prelim)+"\n      At "+str(level)+" level\n(Ascending abundance)", fontsize=10)
    #legend.get_title().set_fontsize(20)
    plt.tight_layout()
    return fig#,legend

fig =plot_my_dfx2(taxdfnew, filterx, 'Class', 'Families', data_to_plot,  x_level_tax_col_dict)
fig_class =plot_my_dfx2(taxdfnew_class, filterx, 'Class', 'Families', data_to_plot,  x_level_tax_col_dict_class)
fig_genus =plot_my_dfx2(taxdfnew_genus, filterx, 'doesntmatter', 'doesntmatter', data_to_plot,  x_level_tax_col_dict_genus)

#fig.savefig("CAZy_producing_families_OTU_richness.png", dpi=1200)

######################################PAPER##ALL PLOTS IN ONE ####################################################
######################################PAPER##ALL PLOTS IN ONE ####################################################
######################################PAPER##ALL PLOTS IN ONE ####################################################
######################################PAPER##ALL PLOTS IN ONE ####################################################
def all_in_one(fstuple, noncazycoldf, dotplotfilter, formattedf, molarpctdf, rowcoldict1, indf, filterx, coldict,indf2, otudata, global_legend_patches):
    fontsizevar=12
    #split the figure into two gridspec boxes, top and bottom
    fig = plt.figure(figsize=fstuple)#Gridspec format = row (y), col (x)
    gs_left = GridSpec(5, 7, wspace=0, hspace=0.05)
    
    #(gridsize, row, col), (row pos, col pos)
    ax = fig.add_subplot(gs_left[0:5, 0:2])
    
    gs_right = GridSpec(5, 7, wspace=0.5, hspace=0.35)
    ax2 = fig.add_subplot(gs_right[0:3, 3:5])
    ax3 = fig.add_subplot(gs_right[0:3, 5:7])
    ax4 = fig.add_subplot(gs_right[3:5, 3:7])
    #########DOTPLOT#######
    #########DOTPLOT#######
    #########DOTPLOT#######
    ##filter
    tax_less_than_molpct_filter=molarpctdf[~(molarpctdf.sum(axis=1)>dotplotfilter)].index
    molarpctdf.loc['']=[0,0,0,0]
    molarpctdf.columns=[1,3,5,10]
    grey=(0,0,0,0.2)
    multiplier=6
    molarmultiplier=multiplier*30
    alpha=0.5
    #Can also be done using plot and markersize=n
    ###add data with no abundance''''
    formattedf=noncazycoldf.append(formattedf)
    
    #for x in tax_less_than_molpct_filter:
    #    if not x == '': 
     #        if x in formattedf:
      #           formattedf.drop(x)
      
    ###formattedf=formattedf[formattedf.sum(axis=1)>dotplotfilter]
    #####
    yvals=list(range(0,len(formattedf.index)))
    yaxiscount=max(yvals)
    xtix=list(formattedf.columns)
    ytix=list(formattedf.index)
    xvals=list(range(0,len(formattedf.columns)))
    #cols=sns.hls_palette(max(yvals)+1, l=.5, s=1, h=.5)[::-1]
    cols=sns.xkcd_palette(['clear blue'])*(max(yvals)+1)
    colo=sns.xkcd_palette(['clear blue'])
    #add alpha
    cols=[list(cols[x])+[alpha] for x in list(range(0,len(cols)))]
    colo=[list(colo[x])+[alpha] for x in list(range(0,len(colo)))]
    proteomelist=[(1,1),(3,3),(5,5),(10,8)]
    abundance_marker='o'
    cazy_marker='s'
    for index in formattedf.index[::-1]:
        for x in xvals:
            ax.scatter(x, yaxiscount, s=(formattedf.loc[index][x])*multiplier, color=colo, marker=abundance_marker)
        if ytix[yaxiscount] in list(molarpctdf.index):
            for x in proteomelist:
                 ax.scatter(x[1], yaxiscount, s=molarpctdf.loc[ytix[yaxiscount]][x[0]]*molarmultiplier, color=(1,1,1,0), edgecolor='black', linewidths=1, marker=cazy_marker)
                 #assign y axis value to row col dictionary, which can be toggled on off at the bottom of the script
                 for x in rowcoldict1:
                      if index in rowcoldict1[x]:
                          if 'yposition' not in rowcoldict1[x][index]:
                              rowcoldict1[x][index]['yposition']=yaxiscount
        if index == '':
            ax.plot([0, max(xvals)],[yaxiscount, yaxiscount], '--', color='grey')
        yaxiscount+=-1
    
    ax.set_xticks(xvals)
    ax.set_yticks(yvals)
    ax.set_yticklabels(ytix)
    ax.set_xticklabels(['0','1','2','3','4','5','6','8','10','16'])
    ax.set_xlabel("Week", fontsize=fontsizevar, labelpad=0) 
    #ax.set_xticklabels(xtix, fontsize=10)
    ax.tick_params(axis='x',which='major', labelsize=10, pad=-3)
    ax.tick_params(axis='y',which='major', labelsize=9, pad=-6)
    legendlabels=[0.25,0.5, 1, 5, 10, 25]
    #collabels=sns.hls_palette(len(legendlabels), l=.5, s=.9, h=.5)[::-1]
    collabels=sns.xkcd_palette(['clear blue'])*(len(legendlabels))
    collabels=[list(collabels[x])+[alpha] for x in list(range(0,len(collabels)))]
    #if use .plot and ms=n, index the plot e.g. [plt.plot(x, y, z)[0] for x in y]
    #patches=[plt.scatter([],[], s=(x*multiplier), color='blue', label=str(x)+'%') for x in legendlabels]
    patches=[plt.scatter([],[], s=1, color='white', marker=abundance_marker, label='white') if legendlabels[x] == '' else plt.scatter([],[], marker=abundance_marker, s=legendlabels[x]*multiplier, color=collabels[x], label=str(legendlabels[x])+'%') for x in list(range(0,len(legendlabels)))]
    legend=ax.legend(handles=patches, ncol=1,loc='upper left', bbox_to_anchor=(1.0,1), borderaxespad=0, fontsize=fontsizevar-2, title="Abundance", frameon=False)
    legend.get_title().set_fontsize(fontsizevar-2)
    molarlegendlabels=[0.05, 0.1, 0.2, 0.3,  0.4,0.5]
    patches2=[plt.scatter([],[], s=1, marker=cazy_marker, color='white', label='white') if x == '' else plt.scatter([],[], marker=cazy_marker, s=(x*molarmultiplier), color=(1,1,1,0), edgecolor='black', linewidths=1, label=str(x)+'%') for x in molarlegendlabels]
    legend2=ax.legend(handles=patches2, ncol=1, bbox_to_anchor=(1,0.325), loc="lower left", borderaxespad=0, fontsize=fontsizevar-2, title=u"\u2211" r'$\bar{x}$'" CAZy\n     mol%", frameon=False)
    legend2.get_title().set_fontsize(fontsizevar-2)
    ax.add_artist(legend)
    ax.add_artist(legend2)
    ###Third invisible legend...
    #leg = plt.legend(handles=[], bbox_to_anchor=(1.4, 0))
    #leg.get_frame().set_alpha(0) # hide legend
    ##############COLUMN COLOUR HEATPLOT STUFF#######
    offset=0.6
    patches3=[]
    patchlist=[]
    changedict={'Gammaproteobacteria':'γ-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria'}
    #Add row_colour patches
    for x in rowcoldict1:
        for y in rowcoldict1[x]:
            ax.add_patch(plt.Rectangle((-0.55,rowcoldict1[x][y]['yposition']-offset),0.25, 0.85,facecolor=rowcoldict1[x][y]['colour'],clip_on=False,linewidth = 0))
            #if statement to stop duplicate legend handles
            if x not in patchlist:
                if x in changedict:
                     patchlist.append(x)
                     patches3.append(mpatches.Patch(color=rowcoldict1[x][y]['colour'], label=str(changedict[x]).strip('[]')))
                else:
                     patchlist.append(x)
                     patches3.append(mpatches.Patch(color=rowcoldict1[x][y]['colour'], label=str(x).strip('[]')))
    legend3=ax.legend(handles=patches3, ncol=1,bbox_to_anchor=(1.0,-0.075), loc="lower left", borderaxespad=0, fontsize=fontsizevar-4.2, title='Class', frameon=False, handlelength=0.25, handleheight=.75, columnspacing=0.2, handletextpad=0.1,labelspacing=0.15)
    legend3.get_title().set_fontsize(fontsizevar-2)
    #########BAR1#######
    #########BAR1#######
    #########BAR1#######
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    sns.set(style="white")
    width=0.8
    #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax2.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=coldict[indf.iloc[x].name], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax2.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax2.set_ylim([0,100.5])
    ax2.set_xticks(position)
    ax2.set_xticklabels(list(indf.columns))
    #sns.despine(top=True, right=True)
    ax2.tick_params(direction='out', length=0, width=0)
    ax2.tick_params(axis='y',which='major', labelsize=10, pad=2)
    ax2.tick_params(axis='x',which='major', labelsize=fontsizevar, pad=2)
    ax2.yaxis.labelpad=0
    ax2.set_ylabel("Relative abundance (%)", fontsize=fontsizevar)
    ax2.set_xlabel('Week', fontsize=12, labelpad=0)
    ax2.set_yticks([20,40,60,80])
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax2.get_legend_handles_labels()
   # legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=2,title="Taxonomy (≥"+str(filterx)+"%)\n(Ascending abundance)", fontsize=20)
  #  legend.get_title().set_fontsize(20)
    #########BAR_2#######
    #########BAR_2#######
    #########BAR_2#######
    #reorder the dataframe by abundance (make it easier to plot)
    indf2=indf2.reindex(index=indf2.sum(axis=1).rank(ascending=0).sort_values().index)
    indf2=indf2.fillna(value=0)
    width=0.8 
    #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf2.columns)))
    position=list(range(0, len(indf2.columns)))
    ##Remove 'NA' column##
    indf2=indf2.loc[~(indf2.index=='NA')]
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf2.index))):
             ax3.bar(position, indf2.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=coldict[indf2.iloc[x].name], label=indf2.iloc[x].name)
             yvalue_cumulative=np.array(yvalue_cumulative+indf2.iloc[x])
    
    ax3.set_ylim([0,100.5])
    ax3.set_xticks(position)
    ax3.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax3.set_xticklabels(list(indf2.columns))
    sns.set(style="white")
    ax3.spines['top'].set_visible(False)
    #ax2.spines['top'].set_visible(False)
    ax3.tick_params(direction='out', length=0, width=0)  
    ax3.tick_params(axis='both', labelsize=20)
    ax3.set_ylabel("Relative abundance (%)", fontsize=fontsizevar)
    ####boxplot#####
    colopal=sns.xkcd_palette(['clear blue'])*10
    ax3_1=ax3.twinx()
    bp=ax3_1.boxplot(otudata, patch_artist=True, positions=position)
    ax3_1.set_ylabel('OTU richness (n)', fontsize=fontsizevar)
    ax3_1.tick_params(direction='out', length=0, width=0)  
    ax3_1.tick_params(axis='both', labelsize=fontsizevar-2)
    ax3_1.set_ylim([0,2500])
    ax3_1.spines['top'].set_visible(False)
    ax3_1.set_xticks(position)
    ax3_1.set_xticklabels(list(indf2.columns))
    ax3.set_yticks([20,40,60,80])
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
    ax3.tick_params(axis='y',which='major', labelsize=10, pad=2)
    ax3.tick_params(axis='x',which='major', labelsize=fontsizevar, pad=2)
    ax3.yaxis.labelpad=0
    ax3_1.tick_params(axis='y',which='major', labelsize=10, pad=2)
    ax3_1.tick_params(axis='x',which='major', labelsize=fontsizevar, pad=2)
    ax3_1.yaxis.labelpad=0
    ax3.set_xlabel('Week', fontsize=12,labelpad=0)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax3.get_legend_handles_labels()
    h2, l2 = ax3_1.get_legend_handles_labels()
    handles=h2+handles
    labels=l2+labels
    #legend.get_title().set_fontsize(20)
    ###Turn off ax4 to import global legend
    ax4.set_axis_off()
    legend4=ax4.legend(handles=patches4[::-1], loc=(-0.065,-0.215),ncol=3,fontsize=fontsizevar-2, frameon=False, edgecolor='black',handlelength=0.75, handleheight=0.75, columnspacing=0.2, handletextpad=0.1,labelspacing=0.25)#
     #add text
    font_dict= {'family': 'arial','color':  'black','weight': 'bold','size': 16}
    y=.8625
    fig.text(0.1075,y, 'a', fontdict=font_dict)
    fig.text(0.455,y, 'b', fontdict=font_dict)
    fig.text(0.685,y, 'c', fontdict=font_dict)

    ####add colored backgrounds
    alphaz=0.1
    #(xy)z,m = (lowerleft), width, height
    fancybox1 = mpatches.FancyBboxPatch((0.0175,0.09), 0.4, 0.78,zorder=0, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor='b',edgecolor='none',alpha=alphaz, transform=fig.transFigure)
    fancybox2 = mpatches.FancyBboxPatch((0.46,0.09), 0.47, 0.78,zorder=0, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor='r',edgecolor='none',alpha=alphaz, transform=fig.transFigure)
    fig.patches.append(fancybox1)
    fig.patches.append(fancybox2)
    return fig, formattedf, molarpctdf







cazycol=[]
noncazycol=[]
for x in list(df_class_taxa_heatmap.index):
         if x in list(dfcazy.index):
              cazycol.append(x)
         if x not in list(dfcazy.index):
              noncazycol.append(x.strip('\n'))
noncazycol=noncazycol[:-1]
noncazycol=['']+noncazycol


cazycol_genus=[]
noncazycol_genus=[]
for x in list(df_genus_taxa_heatmap.index):
         if x in list(dfcazy_genus.index):
              cazycol_genus.append(x)
         if x not in list(dfcazy_genus.index):
              noncazycol_genus.append(x.strip('\n'))
noncazycol_genus=noncazycol_genus[:-1]
noncazycol_genus=['']+noncazycol_genus



#noncazycol_df=df_class_taxa_heatmap.loc[noncazycol[::-1]
taxa_not_in_otu_table_df=pd.DataFrame(columns=dfcazy.columns, index=noncazycol[::-1])
taxa_not_in_otu_table_df=taxa_not_in_otu_table_df.fillna(0)
dotplotfilter=0.025
filterx=1

#noncazycol_df=df_class_taxa_heatmap.loc[noncazycol[::-1]
taxa_not_in_otu_table_df_genus=pd.DataFrame(columns=dfcazy_genus.columns, index=noncazycol_genus[::-1])
taxa_not_in_otu_table_df_genus=taxa_not_in_otu_table_df_genus.fillna(0)
dotplotfilter_genus=0.025
filterx_genus=1

####legend_plot
#append dagger to cazy producing taxa and order by abunbance to make it visually clear
coldict_supplement={}
ordered_taxa=list(df_filtered_normalised.sum(axis=1).rank(ascending=False).sort_values().index)
new_ordered_taxa=[]
for x in ordered_taxa:
    if x in taxdfnew.index:
        #x2=str(x)+ r"$\bf{" + '†' + "}$".strip('[]')
        x2=str(x)+'†'.strip('[]')
        coldict_supplement[x2]=x_level_tax_col_dict[x]
        new_ordered_taxa.append(x2)
    else:
        coldict_supplement[x.strip('[]')]=x_level_tax_col_dict[x]
        new_ordered_taxa.append(x.strip('[]'))
patches4=[]        
for x in new_ordered_taxa:
         patches4.append(mpatches.Patch(facecolor=coldict_supplement[x], label=str(x), edgecolor='black'))

##########GENUS###############
#append dagger to cazy producing taxa and order by abunbance to make it visually clear
coldict_supplement_genus={}
ordered_taxa_genus=list(df_filtered_normalised_genus.sum(axis=1).rank(ascending=False).sort_values().index)
new_ordered_taxa_genus=[]
for x in ordered_taxa_genus:
    if x in taxdfnew_genus.index:
        #x2=str(x)+ r"$\bf{" + '†' + "}$".strip('[]')
        x2=str(x)+'†'.strip('[]')
        coldict_supplement_genus[x2]=x_level_tax_col_dict_genus[x]
        new_ordered_taxa_genus.append(x2)
    else:
        coldict_supplement_genus[x.strip('[]')]=x_level_tax_col_dict_genus[x]
        new_ordered_taxa_genus.append(x.strip('[]'))
patches4_genus=[]        
for x in new_ordered_taxa_genus:
         patches4_genus.append(mpatches.Patch(facecolor=coldict_supplement_genus[x], label=str(x), edgecolor='black'))
##########GENUS###############
         
         
fig, ax = plt.subplots(figsize=(8.5,3))
ax.set_axis_off()
legend4=ax.legend(handles=patches4, loc=(0,0),ncol=4,fontsize=11, frameon=False, edgecolor='black',handlelength=1, handleheight=1)
plt.show()

fig, t, m=all_in_one((10,7),taxa_not_in_otu_table_df, dotplotfilter,dfcazy.iloc[::-1], df_class_taxa_heatmap, row_color_dict_global, df_filtered_normalised, filterx, x_level_tax_col_dict,taxdfnew, data_to_plot, patches4)
#fig.savefig("PNAS_plots/ALL_IN_ONE_PAPER_10_7_1200dpi_SILVA.png", dpi=1200, bbox_inches='tight', pad_inches = 0)



def all_in_one_genus(fstuple, noncazycoldf, dotplotfilter, formattedf, molarpctdf, rowcoldict1, indf, filterx, coldict,indf2, otudata, global_legend_patches):
    fontsizevar=12
    #split the figure into two gridspec boxes, top and bottom
    fig = plt.figure(figsize=fstuple)#Gridspec format = row (y), col (x)
    gs_left = GridSpec(5, 7, wspace=0, hspace=0.05)
    
    #(gridsize, row, col), (row pos, col pos)
    ax = fig.add_subplot(gs_left[0:5, 0:2])
    
    gs_right = GridSpec(5, 7, wspace=0.5, hspace=0.35)
    ax2 = fig.add_subplot(gs_right[0:3, 3:5])
    ax3 = fig.add_subplot(gs_right[0:3, 5:7])
    ax4 = fig.add_subplot(gs_right[3:5, 3:7])
    #########DOTPLOT#######
    #########DOTPLOT#######
    #########DOTPLOT#######
    ##filter
    tax_less_than_molpct_filter=molarpctdf[~(molarpctdf.sum(axis=1)>dotplotfilter)].index
    molarpctdf.loc['']=[0,0,0,0]
    molarpctdf.columns=[1,3,5,10]
    grey=(0,0,0,0.2)
    multiplier=6
    molarmultiplier=multiplier*30
    alpha=0.5
    #Can also be done using plot and markersize=n
    ###add data with no abundance''''
    formattedf=noncazycoldf.append(formattedf)
    
    #for x in tax_less_than_molpct_filter:
    #    if not x == '': 
     #        if x in formattedf:
      #           formattedf.drop(x)
      
    ###formattedf=formattedf[formattedf.sum(axis=1)>dotplotfilter]
    #####
    yvals=list(range(0,len(formattedf.index)))
    yaxiscount=max(yvals)
    xtix=list(formattedf.columns)
    ytix=list(formattedf.index)
    xvals=list(range(0,len(formattedf.columns)))
    #cols=sns.hls_palette(max(yvals)+1, l=.5, s=1, h=.5)[::-1]
    cols=sns.xkcd_palette(['clear blue'])*(max(yvals)+1)
    colo=sns.xkcd_palette(['clear blue'])
    #add alpha
    cols=[list(cols[x])+[alpha] for x in list(range(0,len(cols)))]
    colo=[list(colo[x])+[alpha] for x in list(range(0,len(colo)))]
    proteomelist=[(1,1),(3,3),(5,5),(10,8)]
    abundance_marker='o'
    cazy_marker='s'
    for index in formattedf.index[::-1]:
        for x in xvals:
            ax.scatter(x, yaxiscount, s=(formattedf.loc[index][x])*multiplier, color=colo, marker=abundance_marker)
        if ytix[yaxiscount] in list(molarpctdf.index):
            for x in proteomelist:
                 ax.scatter(x[1], yaxiscount, s=molarpctdf.loc[ytix[yaxiscount]][x[0]]*molarmultiplier, color=(1,1,1,0), edgecolor='black', linewidths=1, marker=cazy_marker)
                 #assign y axis value to row col dictionary, which can be toggled on off at the bottom of the script
                 for x in rowcoldict1:
                      if index in rowcoldict1[x]:
                          if 'yposition' not in rowcoldict1[x][index]:
                              rowcoldict1[x][index]['yposition']=yaxiscount
        if index == '':
            ax.plot([0, max(xvals)],[yaxiscount, yaxiscount], '--', color='grey')
        yaxiscount+=-1
    
    ax.set_xticks(xvals)
    ax.set_yticks(yvals)
    ax.set_yticklabels(ytix)
    ax.set_xticklabels(['0','1','2','3','4','5','6','8','10','16'])
    ax.set_xlabel("Week", fontsize=fontsizevar, labelpad=0) 
    #ax.set_xticklabels(xtix, fontsize=10)
    ax.tick_params(axis='x',which='major', labelsize=10, pad=-3)
    ax.tick_params(axis='y',which='major', labelsize=9, pad=-6)
    legendlabels=[0.25,0.5, 1, 5, 10, 25]
    #collabels=sns.hls_palette(len(legendlabels), l=.5, s=.9, h=.5)[::-1]
    collabels=sns.xkcd_palette(['clear blue'])*(len(legendlabels))
    collabels=[list(collabels[x])+[alpha] for x in list(range(0,len(collabels)))]
    #if use .plot and ms=n, index the plot e.g. [plt.plot(x, y, z)[0] for x in y]
    #patches=[plt.scatter([],[], s=(x*multiplier), color='blue', label=str(x)+'%') for x in legendlabels]
    patches=[plt.scatter([],[], s=1, color='white', marker=abundance_marker, label='white') if legendlabels[x] == '' else plt.scatter([],[], marker=abundance_marker, s=legendlabels[x]*multiplier, color=collabels[x], label=str(legendlabels[x])+'%') for x in list(range(0,len(legendlabels)))]
    legend=ax.legend(handles=patches, ncol=1,loc='upper left', bbox_to_anchor=(1.0,1), borderaxespad=0, fontsize=fontsizevar-2, title="Abundance", frameon=False)
    legend.get_title().set_fontsize(fontsizevar-2)
    molarlegendlabels=[0.05, 0.1, 0.2, 0.3,  0.4,0.5]
    patches2=[plt.scatter([],[], s=1, marker=cazy_marker, color='white', label='white') if x == '' else plt.scatter([],[], marker=cazy_marker, s=(x*molarmultiplier), color=(1,1,1,0), edgecolor='black', linewidths=1, label=str(x)+'%') for x in molarlegendlabels]
    legend2=ax.legend(handles=patches2, ncol=1, bbox_to_anchor=(1,0.525), loc="lower left", borderaxespad=0, fontsize=fontsizevar-2, title=u"\u2211" r'$\bar{x}$'" CAZy\n     mol%", frameon=False)
    legend2.get_title().set_fontsize(fontsizevar-2)
    ax.add_artist(legend)
    ax.add_artist(legend2)
    ###Third invisible legend...
    #leg = plt.legend(handles=[], bbox_to_anchor=(1.4, 0))
    #leg.get_frame().set_alpha(0) # hide legend
    ##############COLUMN COLOUR HEATPLOT STUFF#######
    offset=0.6
    patches3=[]
    patchlist=[]
    changedict={'Gammaproteobacteria':'γ-proteobacteria', 'Deltaproteobacteria':'δ-proteobacteria'}
    #Add row_colour patches
    for x in rowcoldict1:
        for y in rowcoldict1[x]:
            ax.add_patch(plt.Rectangle((-0.55,rowcoldict1[x][y]['yposition']-offset),0.25, 0.85,facecolor=rowcoldict1[x][y]['colour'],clip_on=False,linewidth = 0))
            #if statement to stop duplicate legend handles
            if x not in patchlist:
                if x in changedict:
                     patchlist.append(x)
                     patches3.append(mpatches.Patch(color=rowcoldict1[x][y]['colour'], label=str(changedict[x]).strip('[]')))
                else:
                     patchlist.append(x)
                     patches3.append(mpatches.Patch(color=rowcoldict1[x][y]['colour'], label=str(x).strip('[]')))
    legend3=ax.legend(handles=patches3, ncol=1,bbox_to_anchor=(1.0,-0.05), loc="lower left", borderaxespad=0, fontsize=fontsizevar-4.7, title='Family/Class', frameon=False, handlelength=0.25, handleheight=.75, columnspacing=0.2, handletextpad=0.1,labelspacing=0.15)
    legend3.get_title().set_fontsize(fontsizevar-2)
    #########BAR1#######
    #########BAR1#######
    #########BAR1#######
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    sns.set(style="white")
    width=0.8
    #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax2.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=coldict[indf.iloc[x].name], label=indf.iloc[x].name)
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax2.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax2.set_ylim([0,100.5])
    ax2.set_xticks(position)
    ax2.set_xticklabels(list(indf.columns))
    #sns.despine(top=True, right=True)
    ax2.tick_params(direction='out', length=0, width=0)
    ax2.tick_params(axis='y',which='major', labelsize=10, pad=2)
    ax2.tick_params(axis='x',which='major', labelsize=fontsizevar, pad=2)
    ax2.yaxis.labelpad=0
    ax2.set_ylabel("Relative abundance (%)", fontsize=fontsizevar)
    ax2.set_xlabel('Week', fontsize=12, labelpad=0)
    ax2.set_yticks([20,40,60,80])
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax2.get_legend_handles_labels()
   # legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0.99,1.025), ncol=2,title="Taxonomy (≥"+str(filterx)+"%)\n(Ascending abundance)", fontsize=20)
  #  legend.get_title().set_fontsize(20)
    #########BAR_2#######
    #########BAR_2#######
    #########BAR_2#######
    #reorder the dataframe by abundance (make it easier to plot)
    indf2=indf2.reindex(index=indf2.sum(axis=1).rank(ascending=0).sort_values().index)
    indf2=indf2.fillna(value=0)
    width=0.8 
    #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf2.columns)))
    position=list(range(0, len(indf2.columns)))
    ##Remove 'NA' column##
    indf2=indf2.loc[~(indf2.index=='NA')]
    ### Plot the amount of data plotted
    ####Plot bars
    for x in list(range(0, len(indf2.index))):
             ax3.bar(position, indf2.iloc[x],width, bottom=yvalue_cumulative, linewidth=.5, edgecolor='black', color=coldict[indf2.iloc[x].name], label=indf2.iloc[x].name)
             yvalue_cumulative=np.array(yvalue_cumulative+indf2.iloc[x])
    
    ax3.set_ylim([0,100.5])
    ax3.set_xticks(position)
    ax3.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax3.set_xticklabels(list(indf2.columns))
    sns.set(style="white")
    ax3.spines['top'].set_visible(False)
    #ax2.spines['top'].set_visible(False)
    ax3.tick_params(direction='out', length=0, width=0)  
    ax3.tick_params(axis='both', labelsize=20)
    ax3.set_ylabel("Relative abundance (%)", fontsize=fontsizevar)
    ####boxplot#####
    colopal=sns.xkcd_palette(['clear blue'])*10
    ax3_1=ax3.twinx()
    bp=ax3_1.boxplot(otudata, patch_artist=True, positions=position)
    ax3_1.set_ylabel('OTU richness (n)', fontsize=fontsizevar)
    ax3_1.tick_params(direction='out', length=0, width=0)  
    ax3_1.tick_params(axis='both', labelsize=fontsizevar-2)
    ax3_1.set_ylim([0,2500])
    ax3_1.spines['top'].set_visible(False)
    ax3_1.set_xticks(position)
    ax3_1.set_xticklabels(list(indf2.columns))
    ax3.set_yticks([20,40,60,80])
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
    ax3.tick_params(axis='y',which='major', labelsize=10, pad=2)
    ax3.tick_params(axis='x',which='major', labelsize=fontsizevar, pad=2)
    ax3.yaxis.labelpad=0
    ax3_1.tick_params(axis='y',which='major', labelsize=10, pad=2)
    ax3_1.tick_params(axis='x',which='major', labelsize=fontsizevar, pad=2)
    ax3_1.yaxis.labelpad=0
    ax3.set_xlabel('Week', fontsize=12,labelpad=0)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax3.get_legend_handles_labels()
    h2, l2 = ax3_1.get_legend_handles_labels()
    handles=h2+handles
    labels=l2+labels
    #legend.get_title().set_fontsize(20)
    ###Turn off ax4 to import global legend
    ax4.set_axis_off()
    legend4=ax4.legend(handles=patches4_genus[::-1], loc=(-0.065,-0.175),ncol=3,fontsize=fontsizevar-3, frameon=False, edgecolor='black',handlelength=0.75, handleheight=0.75, columnspacing=0.05, handletextpad=0.1,labelspacing=0.25)#
     #add text
    font_dict= {'family': 'arial','color':  'black','weight': 'bold','size': 16}
    y=.8625
    fig.text(0.1075,y, 'a', fontdict=font_dict)
    fig.text(0.455,y, 'b', fontdict=font_dict)
    fig.text(0.685,y, 'c', fontdict=font_dict)

    ####add colored backgrounds
    alphaz=0.1
    #(xy)z,m = (lowerleft), width, height
    fancybox1 = mpatches.FancyBboxPatch((0.0175,0.09), 0.4, 0.78,zorder=0, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor='b',edgecolor='none',alpha=alphaz, transform=fig.transFigure)
    fancybox2 = mpatches.FancyBboxPatch((0.46,0.09), 0.47, 0.78,zorder=0, boxstyle=mpatches.BoxStyle("Round", pad=0.02), facecolor='r',edgecolor='none',alpha=alphaz, transform=fig.transFigure)
    #fig.patches.append(fancybox1)
    #fig.patches.append(fancybox2)
    return fig, formattedf, molarpctdf



fig_genus, t_genus, m_genus=all_in_one_genus((10,12),taxa_not_in_otu_table_df_genus, dotplotfilter_genus,dfcazy_genus.iloc[::-1], df_genus_taxa_heatmap, row_color_dict_global_genus, df_filtered_normalised_genus, filterx_genus, x_level_tax_col_dict_genus,taxdfnew_genus, data_to_plot, patches4_genus)
##fig_genus.savefig("PNAS_plots/ALL_IN_ONE_PAPER_10_7_1200dpi_SILVA_GENUS.png", dpi=1200, bbox_inches='tight', pad_inches = 0)

fig_genusF, t_genusF, m_genusF=all_in_one_genus((10,12),taxa_not_in_otu_table_df_genus, dotplotfilter_genus,dfcazy_genus.iloc[::-1], df_genus_taxa_heatmap, row_color_dict_global_genus_family, df_filtered_normalised_genus, filterx_genus, x_level_tax_col_dict_genus,taxdfnew_genus, data_to_plot, patches4_genus)
##fig_genusF.savefig("PNAS_plots/ALL_IN_ONE_PAPER_10_7_1200dpi_SILVA_GENUS_F.png", dpi=1200, bbox_inches='tight', pad_inches = 0)

















####legend_plot
#append dagger to cazy producing taxa and order by abunbance to make it visually clear
coldict_supplement={}
ordered_taxa=list(df_filtered_normalised.sum(axis=1).rank(ascending=False).sort_values().index)
new_ordered_taxa=[]
for x in ordered_taxa:
    if x in taxdfnew.index:
        #x2=str(x)+ r"$\bf{" + '†' + "}$".strip('[]')
        x2=str(x)+'†'.strip('[]')
        coldict_supplement[x2]=x_level_tax_col_dict[x]
        new_ordered_taxa.append(x2)
    else:
        coldict_supplement[x.strip('[]')]=x_level_tax_col_dict[x]
        new_ordered_taxa.append(x.strip('[]'))
patches4=[]        
for x in new_ordered_taxa:
         patches4.append(mpatches.Patch(facecolor=coldict_supplement[x], label=str(x), edgecolor='black'))

fig, ax = plt.subplots(figsize=(8.5,3))
ax.set_axis_off()
legend4=ax.legend(handles=patches4, loc=(0,0),ncol=4,fontsize=11, frameon=False, edgecolor='black',handlelength=1, handleheight=1)
plt.show()
#fig.savefig('Comprof_legend_only_1%.png', dpi=1200, bbox_inches='tight')
      




########################################## N M D S nmds N M D S #################################################    
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 
###########BRAY CRUTIS MATRIX#################
#filter otu_table
otu_table_xform=otu_table#.transform('sqrt')
otu_table_f, otu_table_n=tax_normalise_filter(otu_table_xform, 0)
#fill NaNs
otu_table_f=otu_table_f.fillna(0)

from scipy.spatial import distance

#distance.braycurtis(otu_table['10A-16S'], otu_table['10B-16S'])

sample_list=list(otu_table.columns)
iteration=0
bray_curtis_matrix=[]
#store bray curtis values as list of lists, by looping over all samples in a particular order
target_sample=0
for target in sample_list:
     templist=[]
     for x in list(range(0, len(sample_list))):
          templist.append(distance.braycurtis(otu_table_f[target], otu_table_f[sample_list[x]]))
     bray_curtis_matrix.append(templist)
    
from sklearn import manifold

stress=1e-12
seed = np.random.RandomState(seed=3)

#generate nmds scaffold
nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=stress,
                    dissimilarity="precomputed", random_state=seed, n_jobs=1,
                    n_init=1)
#positioning for nmds
npos = nmds.fit_transform(bray_curtis_matrix, init=None)

        
#plt.scatter(npos[:, 0], npos[:, 1], color='darkorange', s=50, lw=0, label='NMDS')
marker_list=('x','o', '^','H','o', '^','H','o', '^','H','o')
fig, ax = plt.subplots(figsize=(6,6))
colorlist=sns.xkcd_palette(["clear blue", 'red', 'orchid', 'green', 'coral', 'grey', 'yellow', 'shocking pink', 'salmon', 'black'])
colorlist=sns.hls_palette(10)
label_list=['W10','W16','W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'W8', 'D0']
reps=5
rep0=0 
ax.tick_params(axis='both', labelsize=14, direction='in', length=1, width=1, color='black')
ax.set_xlabel('MDS1', fontsize=16)    
ax.set_ylabel('MDS2', fontsize=16)  
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
for x in list(range(0, 10)):
    print(label_list[x])
    if reps == 5:
         ax.scatter(npos[:reps][:, 0], npos[:reps][:, 1], color=colorlist[x],marker=marker_list[x],edgecolors='black', s=50, lw=1.5, label=label_list[x])
    if reps >5 and reps <=45:
         ax.scatter(npos[rep0:reps][:, 0], npos[rep0:reps][:, 1], color=colorlist[x],marker=marker_list[x],edgecolors='black', s=50, lw=1.5, label=label_list[x])
    if reps == 50: 
         #D0 timepoint
         ax.scatter(npos[-1:][:, 0], npos[-1:][:, 1], color=colorlist[x], s=50, lw=1.5,marker=marker_list[x],edgecolors='black', label=label_list[x])
    reps+=5
    rep0+=5

plt.legend(loc=(.95,0.25), ncol=1, frameon=True)

########################################## N M D S nmds N M D S #################################################    
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 
########################################## N M D S nmds N M D S ################################################# 


###############################################################################
#############################Look at Gammaproteobacteria not normalised#################################
filterx=0
harvest_tax='Gammaproteobacteria'
harvest_tax_level='Class' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Family'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Family')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(a, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
#############################Look at VIBRIO#################################
filterx=0
harvest_tax='Vibrionaceae'
harvest_tax_level='Family' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Genus'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Genus')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)
#fig.savefig('nitrogen_fixers_vibrio_genus_norm', dpi=500)


#############################Look at Alteromonadaceae#################################
filterx=0
harvest_tax='Alteromonadaceae'
harvest_tax_level='Family' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Genus'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Genus')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)


#############################Look at Saccharospirillaceae#################################
filterx=0
harvest_tax='Saccharospirillaceae'
harvest_tax_level='Family' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Genus'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Genus')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)

#############################Look at Cellvibrionaceae#################################
filterx=0
harvest_tax='Cellvibrionaceae'
harvest_tax_level='Family' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Genus'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Genus')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)

#############################Look at Firmicutes#################################
filterx=0
harvest_tax='Firmicutes'
harvest_tax_level='Phylum' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Family'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Family')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)

#############################Look at Deltaproteobacteria#################################
filterx=0
harvest_tax='Deltaproteobacteria'
harvest_tax_level='Class' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Family'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Family')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)

#############################Look at Bacteroidetes#################################
filterx=0
harvest_tax='Bacteroidetes'
harvest_tax_level='Phylum' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Family'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Family')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)

#############################Look at Prolixibacteraceae #################################
filterx=0
harvest_tax='Prolixibacteraceae'
harvest_tax_level='Family' # tax level of harvest tax
#focal_taxonomy is the taxonomic level within harvest tax you want to look at
focal_taxonomy='Genus'
Gdf=tax_table[tax_table[harvest_tax_level]==harvest_tax]
GDFtotal=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, 'Genus')
reltot, reltoterr = get_me_relative_stats(GDFtotal, merged_otu_table, otu_table)

Gdff=tax_extract(merged_otu_table.loc[Gdf.index], Gdf, focal_taxonomy)
a, Gdff_filtered_normalised=tax_normalise_filter(Gdff, filterx)   

linedict={}

fig, lgd=plot_my_df_with_line_no_ylim(Gdff_filtered_normalised, filterx, 'yes', harvest_tax, focal_taxonomy, linedict)




################# ALL BIOLOGICAL REPLICATES ##############################
def color_jenga_all(df):
    newpal=sns.xkcd_palette(["clear blue", "coral", "light aqua", "pale yellow", "light blue", "pinkish red", "white", "light grey", "black"])
    length=int(len(df.index)/2)+1
    patterns=['-', '+', 'x', '\\', '*', 'o', 'O', '.']
    patterns=patterns*(int(length/len(patterns))+1)
    hatches=[]
    #plus 1 so as to not get caught out when int round 0.45 and 0.49 to 0...
    lightpal=sns.hls_palette(length, l=.3, s=.8, h=.5)
    darkpal=sns.hls_palette(length, l=.6, s=.9, h=0)
    temp_pal=[]
    for x in list(range(0, length)):
        temp_pal.append(lightpal[x])
        temp_pal.append(darkpal[x])
        hatches.append(patterns[x])
        hatches.append('')       
    random.shuffle(temp_pal)
    newpal=newpal+temp_pal
    return newpal, hatches



def plot_my_df_ALL(indf, filterx, cols, colorder):
    colors, hatch_list=color_jenga_all(indf)
    #reorder the dataframe by abundance (make it easier to plot)
    indf=indf.reindex(index=indf.sum(axis=1).rank(ascending=0).sort_values().index)
    indf=indf.fillna(value=0)
    fig, ax = plt.subplots(figsize=(12,12))
    sns.set(style="white")
    width=0.8
 #Loop over index and plot the values
    yvalue_cumulative=np.array([0]*len(list(indf.columns)))
    position=list(range(0, len(indf.columns)))
    ### Plot the amount of data plotted
    ###sort columns
    indf=indf[cols]
    ####Plot bars
    for x in list(range(0, len(indf.index))):
        ax.bar(position, indf.iloc[x],width, bottom=yvalue_cumulative, linewidth=2, edgecolor='black', color=colors[x], label=indf.iloc[x].name, hatch=hatch_list[x])
        yvalue_cumulative=np.array(yvalue_cumulative+indf.iloc[x])
    ax.set_xlim([0-width/2, max(position)+width/2+0.05])
    ax.set_ylim([0,100.5])
    ax.set_xticks(position)
    ax.set_xticklabels(colorder, rotation=90)
    sns.despine(top=True, right=True)
    ax.tick_params(direction='out', length=6, width=4)  
    ax.tick_params(axis='both', labelsize=16)
    ax.set_ylabel("Relative abundance (%)", fontsize=16)
    ax.set_xlabel('Week', fontsize=16)
    #call handles and labels, so the legend is inverted matching the vertical bar structure
    handles, labels = ax.get_legend_handles_labels()
    legend=ax.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(0,-0.08), ncol=3,title="Taxonomy (≥"+str(filterx)+"%)\n(Ascending abundance)", fontsize=14, columnspacing=0.1, frameon=True, handleheight=1.5)
    legend.get_title().set_fontsize(12)
    return fig, legend

filterx=1
column_order=['D0-16S','1A-16S', '1B-16S', '1C-16S', '1D-16S', '1E-16S','2A-16S', '2B-16S', '2C-16S', '2D-16S', '2E-16S','3A-16S', '3B-16S', '3C-16S', '3D-16S', '3E-16S','4A-16S', '4B-16S', '4C-16S', '4D-16S', '4E-16S','5A-16S', '5B-16S', '5C-16S', '5D-16S', '5E-16S','6A-16S', '6B-16S', '6C-16S', '6D-16S', '6E-16S','8A-16S', '8B-16S', '8C-16S', '8D-16S', '8E-16S','10A-16S', '10B-16S', '10C-16S', '10D-16S', '10E-16S','16A-16S', '16B-16S', '16C-16S', '16D-16S', '16E-16S']
column_shorts=[x.split('-')[0] for x in column_order]
df=tax_extract(otu_table, tax_table, 'Family')
df_filtered_otu_familyx, df_filtered_normalised_familyx=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df_ALL(df_filtered_normalised_familyx, filterx,column_order, column_shorts)
#fig.savefig("PNAS_plots/PNAS_supplementary_families_technical_replicates.png", bbox_inches='tight', dpi=1200)


filterx=0.5
column_order=['D0-16S','1A-16S', '1B-16S', '1C-16S', '1D-16S', '1E-16S','2A-16S', '2B-16S', '2C-16S', '2D-16S', '2E-16S','3A-16S', '3B-16S', '3C-16S', '3D-16S', '3E-16S','4A-16S', '4B-16S', '4C-16S', '4D-16S', '4E-16S','5A-16S', '5B-16S', '5C-16S', '5D-16S', '5E-16S','6A-16S', '6B-16S', '6C-16S', '6D-16S', '6E-16S','8A-16S', '8B-16S', '8C-16S', '8D-16S', '8E-16S','10A-16S', '10B-16S', '10C-16S', '10D-16S', '10E-16S','16A-16S', '16B-16S', '16C-16S', '16D-16S', '16E-16S']
column_shorts=[x.split('-')[0] for x in column_order]
df=tax_extract(otu_table, tax_table, 'Genus')
df_filtered_otu_familyx, df_filtered_normalised_familyx=tax_normalise_filter(df, filterx)
fig, lgd=plot_my_df_ALL(df_filtered_normalised_familyx, filterx,column_order, column_shorts)
#fig.savefig("PNAS_plots/PNAS_supplementary_genus_technical_replicates.png",bbox_inches='tight', dpi=1200)



give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Alteromonadaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Vibrionaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Flavobacteriaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Cellvibrionaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Prolixibacteraceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Sphingobacteriaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Bacteroidaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Rhodobacteraceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Marinilabiliaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

give_me_Xleveltaxa_that_are_cazy_producing_from_X_level('Saccharospirillaceae', 'Family', tax_table, 'Genus', df_genus_taxa_heatmap)

set(tax_table[tax_table['Genus']=='Teredinibacter']['Family'])
set(tax_table[tax_table['Genus']=='Sporocytophaga']['Family'])
set(tax_table[tax_table['Genus']=='Aquimarina']['Family'])
set(tax_table[tax_table['Genus']=='Hyunsooleela']['Family'])
set(tax_table[tax_table['Genus']=='Planococcus']['Family'])
set(tax_table[tax_table['Genus']=='Pseudosphingobacterium']['Family'])
set(tax_table[tax_table['Genus']=='Desulfosporosinus']['Family'])
set(tax_table[tax_table['Genus']=='Simiduia']['Family'])
set(tax_table[tax_table['Genus']=='Sorangium']['Family'])
set(tax_table[tax_table['Genus']=='Lentimicrobium ']['Family'])
set(tax_table[tax_table['Genus']=='Arcticbacter']['Family'])
set(tax_table[tax_table['Genus']=='Formosa']['Family'])
v

figg, axg = plt.subplots(figsize=(7,25))
sns.set(font_scale=0.8)
hmgen=sns.clustermap(df_genus_taxa_heatmap.drop(['NA', 'Vibrio']), col_cluster=False, annot=False, annot_kws={"size": 14}, figsize=(7,16), cbar_kws={'label': 'mol%'}, cmap='RdBu', linewidths=.5)

