import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as gridspec
import bbi
import bioframe
import matplotlib as mpl

#This function is used to assigns values 0-1 to P arm and values 1-0 or 1-2 to arm Q
# 0 represents telomeric regions and 1 is centromeric region
# 0 and 2 are telomeric regions and 1 is centromeric region
def create_Centromere_Telomere_vector(fl,cis,ct,gn,op,mp=False):
    long = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12']
    short = ['chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']
    #read EV file cis/trans
    f1 = pd.read_csv(cis+fl+"_"+gn+".cis.vecs.tsv", sep='\t', header=0)
    #get the columns required and add column CT to maindf
    f1 = f1.reindex(columns = ['chrom','start','end','E1','CT']) 
    #file having centromeric regions
    file_c = pd.read_csv(ct, sep='\t', header=0)
    for index,row in file_c.iterrows():
        #extract rows for single chr arm 1
        if mp:
            f2=f1[(f1['chrom']==row['chrom']) & (f1['end']<row['midpoint'])]
        else:
            f2=f1[(f1['chrom']==row['chrom']) & (f1['end']<=row['chromStart'])]
        
        if(len(f2)!=0):
            # assign values between 0-1 for all bins , 0-telomere and 1 centromere
            x=np.arange(0,1,1/len(f2))
            #get index value from main df to fill CT column
            s=f2.index[0]
            #set CT column 
            for i in range(len(f2)):
                f1.iloc[s+i,f1.columns.get_loc('CT')]=x[i]
        
        #extract rows for single chr arm 2
        if mp:
            f3=f1[(f1['chrom']==row['chrom']) & (f1['end']>row['midpoint'])]
        else:
            f3=f1[(f1['chrom']==row['chrom']) & (f1['start']>row['chromEnd'])]
        
        if(len(f3)!=0):
            #assign values 1-0 to arm q, 1 being centromere
            #y=np.flip(np.arange(0,1,1/len(f3)))
            #assign values 1-2, 1 being centromere and 2 is telomere
            y=np.arange(1,2,1/len(f3))
            s=f3.index[0]
            for i in range(len(f3)):
                f1.iloc[s+i,f1.columns.get_loc('CT')]=y[i]

    flong = f1[f1['chrom'].isin(long)]
    fshort = f1[f1['chrom'].isin(short)]
    f1.to_csv(op+fl+"_"+gn+".Cent.vecs.tsv",sep='\t',index=False)
    flong.to_csv(op+fl+"_"+gn+".longCent.vecs.tsv",sep='\t',index=False)
    fshort.to_csv(op+fl+"_"+gn+".shortCent.vecs.tsv",sep='\t',index=False)
