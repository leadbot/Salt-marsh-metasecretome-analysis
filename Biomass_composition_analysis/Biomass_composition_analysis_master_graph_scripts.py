# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 21:20:08 2018

@author: Dan
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
##########BIOMASS COMPOSITION ANALYSIS BARPLOTS

cell=[155.9637819,155.8700673,156.0859155,174.077912,190.0748804,209.5865546,169.664331,130.7570381,163.1465918,141.1761482]
cellerr=[5.443618197,4.395747838,4.700194006,4.818834271,5.594223348,5.81779328,9.291900992,10.6773415,6.23139119,3.17359572]


lig=[269.971831,239.6056338,205.1906103,254.3887324,261.4985915,286.1671362,273.8666667,236.7549296,314.5765258,294.2046948]
ligerr=[1.854796428,4.648269128,9.727393753,10.17262299,7.02224129,7.257130872,13.22727854,13.64159479,12.94147343,6.077445456]

ash=[92.15770148,205.4031225,284.1049052,139.6713576,117.6029625,127.3158398,185.5361049,292.3638708,225.2337059,226.2591326]
asherr=[7.397486854,31.34257549,22.91783689,16.58743141,14.39844853,3.296678966,20.17558755,50.16256722,35.89376022,21.53962017]

weeks=['0','1','2','3','4','5','6','8','10','16']
X=list(range(0,len(weeks)))
colour_barplot=sns.xkcd_palette(['light violet'])
colour_barplot2=sns.xkcd_palette(['ivory'])
colour_barplot3=sns.xkcd_palette(['light grey'])
sns.set_style('white')

errorbar_format=dict(ecolor='black', lw=1, capsize=3, capthick=1)
#Use gridspec ratios to ensure bars are the same size
fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(24,6),dpi=700) 
C=ax.bar(X,cell,edgecolor='black',width=0.8, yerr=cellerr, error_kw=errorbar_format, tick_label=weeks, color=colour_barplot)
for bar in C:
    bar.set_edgecolor("black")
    bar.set_linewidth(2)
ax.set_ylabel('Crystalline cellulose\n(µg glucose.mg biomass$^-$$^1$)', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax.set_xlim([0-.4, len(weeks)-.4])
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=16)
###LIGNIN
L=ax2.bar(X,lig,width=0.8, yerr=ligerr, error_kw=errorbar_format, tick_label=weeks, color=colour_barplot2)
for bar in L:
    bar.set_edgecolor("black")
    bar.set_linewidth(2)
ax2.set_ylabel('Acetyl bromide soluble lignin\n(µg.mg biomass$^-$$^1$)', fontsize=20)
ax2.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax2.set_xlim([0-.4, len(weeks)-.4])
ax2.tick_params(direction='out', length=3, width=2)  
ax2.yaxis.tick_left()
ax2.tick_params(labelright='off')
ax2.tick_params(direction='out', length=3, width=2)
ax2.tick_params(axis='both', labelsize=16)
###ASH
A=ax3.bar(X,ash,width=0.8, yerr=asherr, error_kw=errorbar_format, tick_label=weeks, color=colour_barplot3)
for bar in A:
    bar.set_edgecolor("black")
    bar.set_linewidth(2)
ax3.set_ylabel('Ash\n(µg.mg biomass$^-$$^1$)', fontsize=20)
ax3.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax3.set_xlim([0-.4, len(weeks)-.4])
ax3.tick_params(direction='out', length=3, width=2)  
ax3.yaxis.tick_left()
ax3.tick_params(labelright='off')
ax3.tick_params(direction='out', length=3, width=2)
ax3.tick_params(axis='both', labelsize=16)
#plt.text(position[2]-0.2, -10, "Substrate", fontsize=16)
plt.tight_layout()
fig.savefig("ug per mg cel_lig_ash.png", dpi=1200)

