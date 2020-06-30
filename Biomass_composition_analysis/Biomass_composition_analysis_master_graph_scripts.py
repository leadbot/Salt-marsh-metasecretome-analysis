# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 21:20:08 2018

@author: Dan
"""
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches
#####Masss remaining
import statsmodels.api as sm
import math
from scipy.optimize import curve_fit
datafile='MASTER_Biomass_composition_analysis_CC_MBTH_DIONEX_ABSL_ASH_aug17.xls'
workbook=xlrd.open_workbook(datafile)
sheet=workbook.sheet_by_name("mass_remaining_python")

biomass_remaining_pct=np.array([x.value for x in sheet.row_slice(1, start_colx=1, end_colx=11)]).astype(np.float)
biomass_remaining_g=np.array([x.value for x in sheet.row_slice(2, start_colx=1, end_colx=11)]).astype(np.float)

biomass_remaining_pct_sterr=np.array([x.value for x in sheet.row_slice(21, start_colx=1, end_colx=11)]).astype(np.float)
biomass_remaining_g_sterr=np.array([x.value for x in sheet.row_slice(22, start_colx=1, end_colx=11)]).astype(np.float)
xpos=np.array([0,1,2,3,4,5,6,8,10,16]).astype(np.float)

errorbar_format=dict(ecolor='black', lw=1, capsize=3, capthick=1)
fillupper_coord=np.array([x+(y) for x,y in zip(biomass_remaining_g, biomass_remaining_g_sterr)]).astype(np.float)
filllower_coord=np.array([x-(y) for x,y in zip(biomass_remaining_g, biomass_remaining_g_sterr)]).astype(np.float)

#NB when writing the function X or Y MUST come first, e.g. x, a, b, c NOT a, b, c, x
def expo(x, a, b):
        return a * np.exp(-b*x)
         
popt, pcov = curve_fit(expo, xpos, biomass_remaining_g)
xcurve= np.linspace(0, 16, 50)

###get R squared value###
residuals = biomass_remaining_g- expo(xpos, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((biomass_remaining_g-np.mean(biomass_remaining_g))**2)
r_squared = 1 - (ss_res / ss_tot)

r2string=str(round(float(popt[0]),2))+'e$^-$'+'$^{'+str(round(float(popt[1]),2))+'}$'+'$^x$\nR$^2$='+str(round(float(r_squared),3))

fig, ax = plt.subplots(figsize=(8,6),dpi=500) 
sns.set_style('white')
ax.errorbar(xpos, biomass_remaining_g, yerr=biomass_remaining_g_sterr, label='Mass remaining', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, color=(0, 0.88, 0))
ax.plot(xcurve, expo(xcurve, *popt), color='black')
ax.fill_between(xpos, filllower_coord, fillupper_coord, color=(0, 0.88, 0, 0.2))
ax.set_ylim([0,51])
ax2=ax.twinx()
ax2.set_ylim([0,101])
ax.set_xlim([0, 16])
ax.set_xticks(xpos)
ax.set_ylabel("Mass (g)", fontsize=24)
ax.set_xlabel('Week', fontsize=24)
ax2.set_ylabel('Mass remaining (%)', fontsize=24)
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.tick_params(direction='out', length=3, width=2)  
ax2.tick_params(direction='out', length=3, width=2)  
ax.tick_params(axis='both', labelsize=20)
ax2.tick_params(axis='both', labelsize=20)
plt.text(7.5, 50, r2string)
plt.tight_layout()
fig.savefig("mass_remaining_only_matplot_8_6.png", dpi=700)

####degradation rate

datafile='MASTER_Biomass_composition_analysis_CC_MBTH_DIONEX_ABSL_ASH_aug17.xls'
workbook=xlrd.open_workbook(datafile)
sheet=workbook.sheet_by_name("mass_remaining_python")

biomass_remaining_g_sterr=[x.value for x in sheet.row_slice(22, start_colx=1, end_colx=11)]
biomass_remaining_g=[x.value for x in sheet.row_slice(2, start_colx=1, end_colx=11)]

biomass_rate_pct_loss_per_day=[4.053307836,	0.348205037	,1.552383179	,0.605395155	,0.431655023	,0.696245725	,0.222152686,	0.379605315	,0.175591362]
biomass_rate_pct_err=[0.22713557	,0.428860519	,0.305440359,	0.132542559,	0.063452483,	0.07779356	,0.08695321,	0.111912123,	0.029670451]
error=[1.333273044,	1.832038333,	1.546107669	,1.01848397,	0.704694425	,0.780276175	,0.824934246,	0.935869061,	0.481879662]

count=0
weekdecay=[]
weeks=['1','2','3','4','5','6','8','10','16']
weekdiv=[1,1,1,1,1,1,2,2,6]
position=np.array([1,2,3,4,5,6,8,10,16])
for x in list(range(0,len(biomass_remaining_g))):
    if count+1 <= len(biomass_remaining_g)-1:
         weekdecay.append(biomass_remaining_g[count+1]-biomass_remaining_g[count])
         count+=1
weekdecay_norm=np.array(weekdecay)/np.array(weekdiv)
sns.set_style("white")
f, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))

# plot the same data on both axes
ax.errorbar(position, biomass_rate_pct_loss_per_day, yerr=biomass_rate_pct_err,ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25)
ax2.errorbar(position, biomass_rate_pct_loss_per_day, yerr=biomass_rate_pct_err,ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25)

# zoom-in / limit the view to different portions of the data
ax.set_ylim(3.6, 4.5)  # outliers only
ax2.set_ylim(-0.25, 2)  # most of the data
ax.text(-1.5, 4.2, "Decay (mass loss (%).day$^-$$^1$)", ha='center', rotation='vertical', fontsize=20)
ax2.set_xlabel("Week", fontsize=24)
ax.set_xticks(position)
ax.set_xticklabels(weeks)
# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop='off')  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
ax.tick_params(axis='y',direction='out', length=3, width=2)
ax2.tick_params(direction='out', length=3, width=2)
ax.set_xticks(position)
ax2.set_xticks(position)
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
f.savefig("decomposition_rates_pct_per_day.png", dpi=1200)

##########BIOMASS COMPOSITION ANALYSIS BARPLOTS

cell=[155.9637819,155.8700673,156.0859155,174.077912,190.0748804,209.5865546,169.664331,130.7570381,163.1465918,141.1761482]
cellerr=[5.443618197,4.395747838,4.700194006,4.818834271,5.594223348,5.81779328,9.291900992,10.6773415,6.23139119,3.17359572]


lig=[269.971831,239.6056338,205.1906103,254.3887324,261.4985915,286.1671362,273.8666667,236.7549296,314.5765258,294.2046948]
ligerr=[1.854796428,4.648269128,9.727393753,10.17262299,7.02224129,7.257130872,13.22727854,13.64159479,12.94147343,6.077445456]

ash=[92.15770148,205.4031225,284.1049052,139.6713576,117.6029625,127.3158398,185.5361049,292.3638708,225.2337059,226.2591326]
asherr=[7.397486854,31.34257549,22.91783689,16.58743141,14.39844853,3.296678966,20.17558755,50.16256722,35.89376022,21.53962017]

totalhemi=[192.17,208.47,224.70,213.14,219.32,270.35,268.03,239.82,251.82,234.37]
totalhemierr=[0,4.61,1.41,3.88,3.14,16.21,5.44,8.09,4.08,3.08]


weeks=['0','1','2','3','4','5','6','8','10','16']
X=list(range(0,len(weeks)))
colour_barplot=sns.xkcd_palette(['light violet'])
colour_barplot2=sns.xkcd_palette(['ivory'])
colour_barplot3=sns.xkcd_palette(['light grey'])
sns.set_style('white')

errorbar_format=dict(ecolor='black', lw=1, capsize=3, capthick=1)
#Use gridspec ratios to ensure bars are the same size

hemicol=sns.xkcd_palette(['clear blue'])
errorbar_format=dict(ecolor='black', lw=1, capsize=3, capthick=1)
#Use gridspec ratios to ensure bars are the same size
fig, ax = plt.subplots(figsize=(8,6),dpi=700) 
C=ax.bar(X,totalhemi,edgecolor='black',width=0.8, yerr=totalhemierr, error_kw=errorbar_format, tick_label=weeks, color=hemicol)
for bar in C:
    bar.set_edgecolor("black")
    bar.set_linewidth(2)
ax.set_ylabel('Total matrix polysaccharides\n(µg monosaccharide.mg biomass$^-$$^1$)', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax.set_xlim([0-.4, len(weeks)-.4])
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=16)
plt.tight_layout()
fig.savefig("total_matrix_polysaccharides.png", dpi=1200)




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

#########################################
#Hemicellulose fractions
xylose=[45.11,	63.34,	69.58,	67.89,	66.53,	71.19,	79.64,	73.17,	81.35,	80.11]
galactose=[33.35,	31.77,	31.50,	28.73,	30.65,	40.69,	38.56,	34.75,	36.37,	32.72]
glucose=[28.49,	29.55,	30.18,	28.34,	28.15,	43.55,	40.06,	33.86,	34.23,	29.96]
rhamnose=[23.99,	29.85,	33.85,	31.34,	31.52,	28.66,	34.88,	32.13,	32.40,	27.71]
arabinose=[24.95,	21.78,	20.14,	21.00,	22.69,	45.42,	37.57,	31.73,	30.91,	28.39]
gala=[14.72,	19.90,	22.04,	18.27,	23.47,	23.42,	23.60,	22.59,	23.16,	22.85]
manose=[14.60,	12.32,	11.57,	9.51,	10.89,	10.00,	9.10,	8.06,	9.28,	9.04]
glua=[5.04,	5.22,	4.41,	6.67,	3.57,	5.31,	2.40,	1.18,	1.71,	1.21]
fucose=[1.91,	2.16,	1.85,	1.75,	2.34,	2.69,	2.85,	2.90,	3.04,	3.03]

xylose_err=[0,	1.84,	0.60,	1.50,	0.93,	5.52,	2.36,	3.74,	2.36,	2.29]
galactose_err=[0,	0.51,	0.27,	0.63,	0.53,	2.60,	1.00,	0.94,	0.86,	1.20]
glucose_err=[0,	0.60,	0.18,	0.45,	0.55,	3.85,	0.76,	0.84,	1.40,	1.35]
rhamnose_err=[0,	1.68,	0.64,	0.64,	1.26,	3.15,	1.07,	1.50,	1.24,	0.79]
arabinose_err=[0,	0.61,	0.70,	0.98,	0.59,	4.75,	1.28,	0.91,	1.62,	1.54]
gala_err=[0,	0.45,	0.70,	2.20,	0.51,	1.22,	0.60,	2.14,	0.90,	1.32]
manose_err=[0,	0.34,	0.24,	0.35,	0.73,	0.41,	0.33,	0.54,	0.31,	0.21]
glua_err=[0,	0.46,	0.52,	0.68,	0.63,	0.18,	0.34,	0.58,	0.30,	0.23]
fucose_err=[0,	0.44,	0.07,	0.05,	0.09,	0.06,	0.15,	0.14,	0.15,	0.12]

total_hemi=

position=list(range(0, len(xylose)))
weeks=['0','1','2','3','4','5','6','8','10', '16']

e2=np.array(xylose)+np.array(galactose)
e3=e2+np.array(glucose)
e4=e3+np.array(rhamnose)
e5=e4+np.array(arabinose)
e6=e5+np.array(gala)
e7=e6+np.array(manose)
e8=e7+np.array(glua)
e9=e8+np.array(fucose)

allcol=[[182, 225, 255],[153,204,255],[153,195,255],[51, 153, 255],[0, 128, 255],[0,102,204],[0, 76, 153],[0, 51, 102],[0, 25, 51]]
def convert(x):
    ret=[]
    for y in x:
        ret.append(y/255)
    return ret

hemicolorz=[]
for x in allcol:
    hemicolorz.append(convert(x))
names=['Xylose', 'Galactose', 'Glucose', 'Rhamnose', 'Arabinose', 'Galacturonic acid', 'Mannose', 'Glucuronic acid', 'Fucose']

patches=[]
count=0
for x in names:
    patches.append(mpatches.Patch(color=hemicolorz[count], label=x))
    count+=1

fig, ax = plt.subplots(figsize=(8,6),dpi=700) 
ax.stackplot(position, xylose, galactose, glucose,rhamnose,arabinose,gala,manose,glua,fucose, colors=hemicolorz)
ax.errorbar(position, xylose, xylose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e2, galactose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e3, glucose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e4, rhamnose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e5, arabinose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e6, gala_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e7, manose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e8, glua_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, e9, fucose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
fig.subplots_adjust(bottom=0.3) 
ax.legend(handles=patches, ncol=3, fontsize=16, loc=(-0.15,-.5))
#ax.bar(position, xylose, yerr=xyloseerr, color=colors[0], error_kw=errorbar_format, width=barWidth, align='center', edgecolor='white', label=labelz[0])
ax.set_ylabel('µg.mg biomass$^-$$^1$', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax.set_xlim([0, len(xylose)-1])
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.set_xticklabels(weeks)
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=16)
fig.savefig("matrix_polysaccharides.png", dpi=1200)

#############################################################################################
#ug per mg for all fractions


fraccolors=colour_barplot=sns.xkcd_palette(['light violet'])+hemicolorz+sns.xkcd_palette(['ivory', 'light grey'])
position=np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 16]).astype(float)
weeks=['0','1','2','3','4','5','6','8','10', '16']
ticks=[0, 1, 2, 3, 4, 5, 6, 8, 10, 16]
namez=['Cellulose', 'Xylose', 'Galactose', 'Glucose', 'Rhamnose', 'Arabinose', 'Galacturonic acid', 'Mannose', 'Glucuronic acid', 'Fucose', 'Lignin', 'Ash']

patches=[]
count=0
for x in namez:
    patches.append(mpatches.Patch(color=fraccolors[count], label=x))
    count+=1

a=np.array(cell)
a1=a+np.array(xylose)
a2=a1+np.array(galactose)
a3=a2+np.array(glucose)
a4=a3+np.array(rhamnose)
a5=a4+np.array(arabinose)
a6=a5+np.array(gala)
a7=a6+np.array(manose)
a8=a7+np.array(glua)
a9=a8+np.array(fucose)
a10=a9+np.array(lig)
a11=a10+np.array(ash)

fig, ax = plt.subplots(figsize=(16,8),dpi=700) 
ax.stackplot(position, cell, xylose, galactose, glucose,rhamnose,arabinose,gala,manose,glua,fucose,lig, ash, colors=fraccolors)

ax.errorbar(position, a, cellerr, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a1, xylose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a2, galactose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a3, glucose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a4, rhamnose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a5, arabinose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a6, gala_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a7, manose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a8, glua_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a9, fucose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a10, ligerr, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a11, asherr, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)

ax.set_ylabel('µg.mg biomass$^-$$^1$', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax.set_xlim([0, 16])
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.set_xticks(ticks)
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=20)
ax.legend(handles=patches, ncol=6, fontsize=17, loc=(-0.1,-.25))
fig, ax = plt.subplots(figsize=(16,8),dpi=700) 
ax.stackplot(position, cell, xylose, galactose, glucose,rhamnose,arabinose,gala,manose,glua,fucose,lig, ash, colors=fraccolors)

ax.errorbar(position, a, cellerr, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a1, xylose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a2, galactose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a3, glucose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a4, rhamnose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a5, arabinose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a6, gala_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a7, manose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a8, glua_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a9, fucose_err, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a10, ligerr, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, a11, asherr, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)

ax.set_ylabel('µg.mg biomass$^-$$^1$', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax.set_xlim([0, 16])
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.set_xticks(ticks)
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=20)
ax.legend(handles=patches, ncol=6, fontsize=17, loc=(-0.1,-.25))
fig.subplots_adjust(bottom=0.175) 
fig.savefig("ug_mg_all_fractions_ash.png", dpi=1200)

############################Proportion fractions########################################################
######################################
#####         ASH EXCLUDED        ####
######################################
import xlrd
datafile='MASTER_Biomass_composition_analysis_CC_MBTH_DIONEX_ABSL_ASH_aug17.xls'
workbook=xlrd.open_workbook(datafile)
sheet=workbook.sheet_by_name("normalised_compositions_python")

###DATA_VALUES####
cellp=[x.value for x in sheet.row_slice(97, start_colx=1, end_colx=11)]
xylosep=[x.value for x in sheet.row_slice(98, start_colx=1, end_colx=11)]
galactosep=[x.value for x in sheet.row_slice(99, start_colx=1, end_colx=11)]
glucosep=[x.value for x in sheet.row_slice(100, start_colx=1, end_colx=11)]
rhamnosep=[x.value for x in sheet.row_slice(101, start_colx=1, end_colx=11)]
arabinosep=[x.value for x in sheet.row_slice(102, start_colx=1, end_colx=11)]
galap=[x.value for x in sheet.row_slice(103, start_colx=1, end_colx=11)]
manosep=[x.value for x in sheet.row_slice(104, start_colx=1, end_colx=11)]
gluap=[x.value for x in sheet.row_slice(105, start_colx=1, end_colx=11)]
fucosep=[x.value for x in sheet.row_slice(106, start_colx=1, end_colx=11)]
ligp=[x.value for x in sheet.row_slice(107, start_colx=1, end_colx=11)]
#ash=[x.value for x in sheet.row_slice(14, start_colx=1, end_colx=None)]
###############################################################################################
#Standard error values
cell_errp=[x.value for x in sheet.row_slice(57, start_colx=1, end_colx=11)]
xylose_errp=[x.value for x in sheet.row_slice(58, start_colx=1, end_colx=11)]
galactose_errp=[x.value for x in sheet.row_slice(59, start_colx=1, end_colx=11)]
glucose_errp=[x.value for x in sheet.row_slice(60, start_colx=1, end_colx=11)]
rhamnose_errp=[x.value for x in sheet.row_slice(61, start_colx=1, end_colx=11)]
arabinose_errp=[x.value for x in sheet.row_slice(62, start_colx=1, end_colx=11)]
gala_errp=[x.value for x in sheet.row_slice(63, start_colx=1, end_colx=11)]
manose_errp=[x.value for x in sheet.row_slice(64, start_colx=1, end_colx=11)]
glua_errp=[x.value for x in sheet.row_slice(65, start_colx=1, end_colx=11)]
fucose_errp=[x.value for x in sheet.row_slice(66, start_colx=1, end_colx=11)]
lignin_errp=[x.value for x in sheet.row_slice(67, start_colx=1, end_colx=11)]
#ash_err=[x.value for x in sheet.row_slice(14, start_colx=1, end_colx=None)]

b=np.array(cellp)
b1=b+np.array(xylosep)
b2=b1+np.array(galactosep)
b3=b2+np.array(glucosep)
b4=b3+np.array(rhamnosep)
b5=b4+np.array(arabinosep)
b6=b5+np.array(galap)
b7=b6+np.array(manosep)
b8=b7+np.array(gluap)
b9=b8+np.array(fucosep)
b10=b9+np.array(ligp)


fig, ax = plt.subplots(figsize=(16,8),dpi=700) 
ax.stackplot(position, cellp, xylosep, galactosep, glucosep,rhamnosep,arabinosep,galap,manosep,gluap,fucosep,ligp, colors=fraccolors)

ax.errorbar(position, b, cell_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b1, xylose_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b2, galactose_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b3, glucose_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b4, rhamnose_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b5, arabinose_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b6, gala_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b7, manose_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b8, glua_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b9, fucose_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, b10, lignin_errp, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)

ax.set_ylabel('Proportion (%)', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax.set_xlim([0, 16])
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.set_xticks(ticks)
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=20)
ax.legend(handles=patches, ncol=6, fontsize=17, loc=(-0.1,-.25))
plt.tight_layout
fig.subplots_adjust(bottom=0.175) 
fig.savefig("normalised_fraction_EXCLUDING_ASH_MATPLOT.png", dpi=1200)

######################################
#####         ASH I    NCLUDED   ####
######################################

cellg=[x.value for x in sheet.row_slice(112, start_colx=1, end_colx=11)]
xyloseg=[x.value for x in sheet.row_slice(113, start_colx=1, end_colx=11)]
galactoseg=[x.value for x in sheet.row_slice(114, start_colx=1, end_colx=11)]
glucoseg=[x.value for x in sheet.row_slice(115, start_colx=1, end_colx=11)]
rhamnoseg=[x.value for x in sheet.row_slice(116, start_colx=1, end_colx=11)]
arabinoseg=[x.value for x in sheet.row_slice(117, start_colx=1, end_colx=11)]
galag=[x.value for x in sheet.row_slice(118, start_colx=1, end_colx=11)]
manoseg=[x.value for x in sheet.row_slice(119, start_colx=1, end_colx=11)]
gluag=[x.value for x in sheet.row_slice(120, start_colx=1, end_colx=11)]
fucoseg=[x.value for x in sheet.row_slice(121, start_colx=1, end_colx=11)]
ligning=[x.value for x in sheet.row_slice(122, start_colx=1, end_colx=11)]
ashg=[x.value for x in sheet.row_slice(123, start_colx=1, end_colx=11)]
###############################################################################################
#Standard error values
cell_errg=[x.value for x in sheet.row_slice(70, start_colx=1, end_colx=11)]
xylose_errg=[x.value for x in sheet.row_slice(71, start_colx=1, end_colx=11)]
galactose_errg=[x.value for x in sheet.row_slice(72, start_colx=1, end_colx=11)]
glucose_errg=[x.value for x in sheet.row_slice(73, start_colx=1, end_colx=11)]
rhamnose_errg=[x.value for x in sheet.row_slice(74, start_colx=1, end_colx=11)]
arabinose_errg=[x.value for x in sheet.row_slice(75, start_colx=1, end_colx=11)]
gala_errg=[x.value for x in sheet.row_slice(76, start_colx=1, end_colx=11)]
manose_errg=[x.value for x in sheet.row_slice(77, start_colx=1, end_colx=11)]
glua_errg=[x.value for x in sheet.row_slice(78, start_colx=1, end_colx=11)]
fucose_errg=[x.value for x in sheet.row_slice(79, start_colx=1, end_colx=11)]
lignin_errg=[x.value for x in sheet.row_slice(80, start_colx=1, end_colx=11)]
ash_errg=[x.value for x in sheet.row_slice(81, start_colx=1, end_colx=11)]

k=np.array(cellg)
k1=k+np.array(xyloseg)
k2=k1+np.array(galactoseg)
k3=k2+np.array(glucoseg)
k4=k3+np.array(rhamnoseg)
k5=k4+np.array(arabinoseg)
k6=k5+np.array(galag)
k7=k6+np.array(manoseg)
k8=k7+np.array(gluag)
k9=k8+np.array(fucoseg)
k10=k9+np.array(ligning)
k11=k10+ashg


fig, ax = plt.subplots(figsize=(16,8),dpi=700) 
ax.stackplot(position, cellg, xyloseg, galactoseg, glucoseg,rhamnoseg,arabinoseg,galag,manoseg,gluag,fucoseg,ligning,ashg, colors=fraccolors)

ax.errorbar(position, k, cell_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k1, xylose_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k2, galactose_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k3, glucose_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k4, rhamnose_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k5, arabinose_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k6, gala_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k7, manose_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k8, glua_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k9, fucose_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k10, lignin_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, k11, ash_errg, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)

ax.set_ylabel('Proportion (%)', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
sns.despine(top=True, right=True)
ax.set_xlim([0, 16])
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.set_xticks(ticks)
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=20)
ax.legend(handles=patches, ncol=6, fontsize=17, loc=(-0.1,-.25))
fig.subplots_adjust(bottom=0.175)

fig.savefig("normalised_fraction_INCLUDING_ASH_MATPLOT.png", dpi=1200)

##########################################################
#####         COMPISITION OF MASS REMAINING           ####
#########################################################


datafile='MASTER_Biomass_composition_analysis_CC_MBTH_DIONEX_ABSL_ASH_aug17.xls'
workbook=xlrd.open_workbook(datafile)
sheet=workbook.sheet_by_name("mass_remaining_python")

###DATA_VALUES####
biomass_remaining_pct=[x.value for x in sheet.row_slice(1, start_colx=1, end_colx=None)]
biomass_remaining_g=[x.value for x in sheet.row_slice(2, start_colx=1, end_colx=None)]
cello=[x.value for x in sheet.row_slice(3, start_colx=1, end_colx=None)]
xyloseo=[x.value for x in sheet.row_slice(4, start_colx=1, end_colx=None)]
galactoseo=[x.value for x in sheet.row_slice(5, start_colx=1, end_colx=None)]
glucoseo=[x.value for x in sheet.row_slice(6, start_colx=1, end_colx=None)]
rhamnoseo=[x.value for x in sheet.row_slice(7, start_colx=1, end_colx=None)]
arabinoseo=[x.value for x in sheet.row_slice(8, start_colx=1, end_colx=None)]
galao=[x.value for x in sheet.row_slice(9, start_colx=1, end_colx=None)]
manoseo=[x.value for x in sheet.row_slice(10, start_colx=1, end_colx=None)]
gluao=[x.value for x in sheet.row_slice(11, start_colx=1, end_colx=None)]
fucoseo=[x.value for x in sheet.row_slice(12, start_colx=1, end_colx=None)]
lignino=[x.value for x in sheet.row_slice(13, start_colx=1, end_colx=None)]
asho=[x.value for x in sheet.row_slice(14, start_colx=1, end_colx=None)]
# Add original data
x=['0', '1', '2', '3', '4', '5', '6', '8','10', '16']
x_rev=x[::-1]
####STD_ERR_VALUES####
biomass_remaining_pct_sterr=[x.value for x in sheet.row_slice(21, start_colx=1, end_colx=None)]
biomass_remaining_g_sterr=[x.value for x in sheet.row_slice(22, start_colx=1, end_colx=None)]

cell_erro=[x.value for x in sheet.row_slice(23, start_colx=1, end_colx=None)]
xylose_erro=[x.value for x in sheet.row_slice(24, start_colx=1, end_colx=None)]
galactose_erro=[x.value for x in sheet.row_slice(25, start_colx=1, end_colx=None)]
glucose_erro=[x.value for x in sheet.row_slice(26, start_colx=1, end_colx=None)]
rhamnose_erro=[x.value for x in sheet.row_slice(27, start_colx=1, end_colx=None)]
arabinose_erro=[x.value for x in sheet.row_slice(28, start_colx=1, end_colx=None)]
gala_erro=[x.value for x in sheet.row_slice(29, start_colx=1, end_colx=None)]
manose_erro=[x.value for x in sheet.row_slice(30, start_colx=1, end_colx=None)]
glua_erro=[x.value for x in sheet.row_slice(31, start_colx=1, end_colx=None)]
fucose_erro=[x.value for x in sheet.row_slice(32, start_colx=1, end_colx=None)]
lignin_erro=[x.value for x in sheet.row_slice(33, start_colx=1, end_colx=None)]
ash_erro=[x.value for x in sheet.row_slice(34, start_colx=1, end_colx=None)]

o=np.array(cello)
o1=o+np.array(xyloseo)
o2=o1+np.array(galactoseo)
o3=o2+np.array(glucoseo)
o4=o3+np.array(rhamnoseo)
o5=o4+np.array(arabinoseo)
o6=o5+np.array(galao)
o7=o6+np.array(manoseo)
o8=o7+np.array(gluao)
o9=o8+np.array(fucoseo)
o10=o9+np.array(lignino)
o11=o10+asho

fillupper_coord=np.array([x+(y) for x,y in zip(biomass_remaining_g, biomass_remaining_g_sterr)]).astype(np.float)
filllower_coord=np.array([x-(y) for x,y in zip(biomass_remaining_g, biomass_remaining_g_sterr)]).astype(np.float)

fig, ax = plt.subplots(figsize=(16,8),dpi=700) 
ax.stackplot(position, cello, xyloseo, galactoseo, glucoseo,rhamnoseo,arabinoseo,galao,manoseo,gluao,fucoseo,lignino,asho, colors=fraccolors)

ax.errorbar(position, o, cell_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o1, xylose_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o2, galactose_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o3, glucose_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o4, rhamnose_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o5, arabinose_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o6, gala_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o7, manose_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o8, glua_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o9, fucose_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o10, lignin_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)
ax.errorbar(position, o11, ash_erro, ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt=None)

ax2=ax.twinx()

ax.errorbar(position, biomass_remaining_g, yerr=biomass_remaining_g_sterr, label='Total biomass remaining', ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, color=(0, 0.88, 0))
ax.fill_between(position, filllower_coord, fillupper_coord, color=(0, 0.88, 0, 0.2))

ax.set_ylabel('Mass remaining (g)', fontsize=20)
ax.set_xlabel('Week', fontsize=20)
#ax2.set_ylabel("Mass remaining (g)", fontsize=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax.set_xlim([0, 16])
ax.set_ylim([0,50])
ax2.set_ylim([0,50])
ax.tick_params(direction='out', length=3, width=2)  
#ax2.tick_params(direction='out', length=3, width=2)
ax2.set_yticks([])  
ax.yaxis.tick_left()
ax.set_xticks(ticks)
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=20)
#ax2.tick_params(axis='both', labelsize=20)

#get fillbetween legend
lines, labels = ax.get_legend_handles_labels()
ax.legend(handles=patches, ncol=6, fontsize=17, loc=(-0.1,-.25))
ax2.legend(lines, labels, fontsize=17, loc=(0.025, .95))
fig.subplots_adjust(bottom=0.175)

fig.savefig("mass remaining stackarea.png", dpi=1200)

#####Ash vs tides
import statsmodels.api as sm
seven_day_tide=[4.268965517,4.227586207,4.234482759,4.25862069,4.303448276]
mean_change_in_ash_w2_to_w6=[78.70178274,-144.4335476,-22.06839505,9.712877311,58.22026507]
mean_change_err=[38.36223643,29.49553009,19.59071211,13.63177026,21.07587802]


4_results=sm.OLS(mean_change_in_ash_w2_to_w6,sm.add_constant(seven_day_tide)).fit()

#print (regression_results.summary())
fig, ax = plt.subplots(figsize=(8,6),dpi=700) 
ax.errorbar(seven_day_tide,mean_change_in_ash_w2_to_w6, yerr=mean_change_err,  ecolor='black', capsize=3, elinewidth=1.5, markeredgewidth=1.25, fmt='o')
k= np.linspace(4.32,4.2,10)
ax.plot(k, k*regression_results.params[1] + regression_results.params[0], c='black')
ax.set_xlabel("Seven day mean ACD (m)", fontsize=20)
ax.set_ylabel('Inter week ash deviation\nµg Ash.mg biomass$^-$$^1$', fontsize=20)
#ax.set_xlim([1,4])
sns.despine(top=True, right=True)
ax.tick_params(direction='out', length=3, width=2)  
ax.yaxis.tick_left()
ax.tick_params(labelright='off')
ax.tick_params(direction='out', length=3, width=2)
ax.tick_params(axis='both', labelsize=20)
equation="y="+str(round(float(regression_results.params[1]),2))+"x+"+str(round(float(regression_results.params[0]),2))+"\nR$^2$="+str(round(float(regression_results.rsquared),2))
plt.text(4.2, -5, equation, fontsize=20, color='black')
#plt.text(position[2]-0.2, -10, "Substrate", fontsize=16)
plt.tight_layout()
fig.savefig("weekly_ash_change_against_seven_day_mean_tide_high", dpi=1200)
