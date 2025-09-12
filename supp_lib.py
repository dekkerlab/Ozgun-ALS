import pandas as pd
import numpy as np
from numpy import savetxt
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cooler
import bioframe
import bbi
import matplotlib as mpl
from collections import OrderedDict
mpl.style.use('seaborn-v0_8-white')
import os

# make matplotlib pdf-s text recognizable by evil-Adobe
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def norm_shoulder(f):
    '''This function normalizes the flanks by subtracting 
    the avg of 3 flanked bin from every row'''
    _ddd=f
    shoulders = np.hstack( [_ddd[:,:3], _ddd[:,-3:]] )
    shoulders = np.nanmean(shoulders, axis=1, keepdims=True)
    _ddd -= shoulders
    return _ddd 

def plot_stackup(i,reg,binsize,flank,bw,im,title,ylim,norm,hlines,gsIns):
    AsyncBounds=reg
    mids = (AsyncBounds['start'] + AsyncBounds['end']) // 2
    #print(type(mids[0]))
    binsize=binsize
    flank = flank
    nbins = flank*2 // binsize
    bwfile = bw
    stacks = bbi.stackup(bwfile, AsyncBounds['chrom'], mids - flank, mids + flank, bins=nbins)
    
    if norm=='shoulders':
        stacks=norm_shoulder(stacks)
    elif norm=='entire mean':
        stacks = stacks-(np.nanmean(stacks,axis=1)[...,None])
    else:
        stacks=stacks
        
    x = np.linspace(-flank/1e6, flank/1e6, nbins)
    cmap = plt.cm.get_cmap('coolwarm')
    cmap.set_bad('#777777')
    im_opts = dict(
        vmin=im[0], 
        vmax=im[1],
        extent=[-flank/1e6, flank/1e6, len(AsyncBounds), 0],
        cmap=cmap )
    
    # heatmap
    ax = ax1 = plt.subplot(gsIns[0, i])
    if not norm:
        X = stacks-(np.nanmean(stacks,axis=1)[...,None])
    else:
        X=stacks
    img = ax.matshow(X, **im_opts, rasterized=True)
    ax.axvline(0, c='grey', lw=0.5)
    ax.grid(False)
    ax.set_aspect('auto')
    ax.tick_params()
    ax.set_title(title)
    
    # summary
    ax = plt.subplot(gsIns[1, i], sharex=ax1)    
    for h in hlines:
        ax.axhline(h, c='#777777', lw=1, ls='--')
    ax.plot(x, np.nanmean(X, axis=0), c='k', lw=2)
    ax.set_xlim(-flank/1e6, flank/1e6)
    ax.xaxis.set_visible(True)
    ax.tick_params(bottom=False,labelbottom=False)
    ax.set_ylim(ylim[0], ylim[1])
    
    # color bar
    cax = plt.subplot(gsIns[2, i])
    cb = plt.colorbar(img, cax=cax, orientation='horizontal')
    cb.locator = mpl.ticker.MaxNLocator(nbins=3)
    cb.update_ticks()
    
def get_EVsorted_regions(bin,flank,reg,bw_file,how,norm):

    '''This function sorts the regions based on EV : 
    There are three option
    1) subtracting means of flank
    2) sort by middle bin
    3) sort by compartment based on changes at the mid bin (a_b,b_a,a_a,b_b)
    4) Default sort is based on mean of entire row 
    
    Return regions, chrom start end
    '''
    ip = pd.read_csv(reg,sep='\t')
    reg=ip[~ip['chrom'].isin(['chrX','chrY','chrM'])].reset_index(drop=True)
    #print("Org:\n",reg)
    
    binsize=bin
    flank = flank
    nbins = flank*2 // binsize
    bwfile = bw_file
    
    AsyncBounds = reg
    mids = (AsyncBounds['start'] + AsyncBounds['end']) // 2
    f = bbi.stackup(bwfile, AsyncBounds['chrom'], mids - flank, mids + flank, bins=nbins)  
    print(np.argwhere(np.isnan(f[:,3])),np.argwhere(np.isnan(f[:,4])))
    if norm:
            f=norm_shoulder(f)
            
    if how =="flank_diff" :                 
        mid=int(f.shape[1]/2)
        f1=np.nanmean(f[:,mid:],axis=1)
        f2=np.nanmean(f[:,:mid],axis=1)
        diff=f1-f2
        idx=np.argsort(diff)
        regs_fin=reg.loc[idx]
    elif how == "mid":
        idx = np.argsort(f[:, f.shape[1]//2])
        regs_fin=reg.loc[idx]
    elif how == "compartment-changes": 
        #savetxt('data.csv', f, delimiter=',')
        idx=np.arange(f.shape[0])
        print(len(idx))
        a_b=idx[((f[:,3]>=0) & (f[:,4]<0))]
        b_a=idx[((f[:,3]<0) & (f[:,4]>=0))]
        a_a=idx[(f[:,3]>=0) & (f[:,4]>=0) ]
        b_b=idx[(f[:,3]<0) & (f[:,4]<0) ]
        print(len(a_b),len(b_a),len(a_a),len(b_b))
        lst=np.concatenate((a_b,b_a,a_a,b_b), axis=None)
        fin=reg.loc[lst].reset_index(drop=True)
        c1=len(a_b)
        c2=c1+len(b_a)
        c3=c2+len(a_a) 
        AsyncBounds = fin
        mids = (AsyncBounds['start'] + AsyncBounds['end']) // 2
        st = bbi.stackup(bwfile, AsyncBounds['chrom'], mids - flank, mids + flank, bins=nbins)
        X=st
        idx1 = np.argsort(X[:c1, 3])[::-1]
        #print(idx1)
        idx2 = np.argsort(X[c1:c2, 3])
        idx2=idx2+c1
        #print(idx2)
        idx3 = np.argsort(X[c2:c3, 3])[::-1]
        idx3=idx3+c2
        #print(idx3)
        idx4 = np.argsort(X[c3:, 3])
        idx4=idx4+c3
        #print(idx4)
        ids=np.concatenate((idx1,idx2,idx3,idx4), axis=None)
        print(len(ids))
        regs_fin=fin.loc[ids]
    else:
        avg=np.nanmean(f[:,:],axis=1)
        idx=np.argsort(avg)
        regs_fin=reg.loc[idx]
    return regs_fin


def plot_stackups_sets(
                  extra_plots,
                  hmss, # will become a dictionary now (or list) ...
                  titles,
                  limss,
                  cmps,
                  norms=None,
                  binsizes=None,
                  extra_order = None,
                  hmss_order = None,
                  fillmissing=False,
                  interpolation="nearest",
                  inch_per_1k_stack= 1.2, #inches per 1000 elements stackup ...
                  fig_fontsize = 50,
                 ):
    """
    plot a buch of stackups ...
    """
    # rewrite everyhting assuming hmss is a dict of stackup groups !
    # groups are plotted on top of each other ...
    
    if extra_plots is None:
        extra_plots = []
    # regardless - claculate number of axes for stackups ...
    num_stackup_groups = len(hmss)
    # pick in every stackup group and see how many are there
    num_stackups = max(len(hmss[k]) for k in hmss)
    num_rows = num_stackups + len(extra_plots)
    # let's figure out - how tall is this stackup
    # get heights of stackups from each groups
    stackup_group_heights = [ len(hmss[k][0]) for k in hmss ]
    stackup_height = sum(stackup_group_heights) * inch_per_1k_stack / 1_000
    service_height = 3.5
    width_per_col = 3.7
    figure_height = stackup_height + service_height
    fig = plt.figure(
        figsize=(width_per_col*num_rows, figure_height + 2*.8 + 2*.8),
        facecolor="white",
        tight_layout=False,
        constrained_layout=False,
    )
    gs = fig.add_gridspec(
        num_stackup_groups+2,
        num_rows,
        width_ratios=[1]*num_rows,
        height_ratios = \
            [0.95*service_height/figure_height] + \
            [(_h/sum(stackup_group_heights))*(figure_height-service_height)/figure_height for _h in stackup_group_heights] + \
            [0.05*service_height/figure_height],
        left=0,
        right=1,
        top=1 - .8/(figure_height + 2*.8 + 2*0.8),
        bottom=.8/(figure_height + 2*.8 + 2*0.8),
        hspace=((len(stackup_group_heights)+2)/figure_height)*.8,
        wspace=0.07,
    )

    ax_profile = {}
    ax_stackup = {}
    ax_xtra = OrderedDict()
    ax_cbar = {}
    # let's define order
    if extra_order is None:
        extra_order = list( range(len(extra_plots)) )
    if hmss_order is None:
        hmss_order = list( range(len(extra_plots), num_rows) )
    # replace following with the pre-defined column indexes ...
    for idx in hmss_order:
        ax_profile[idx] = fig.add_subplot(gs[0,idx])
        ax_stackup[idx] = [fig.add_subplot(gs[_i+1,idx]) for _i in range(num_stackup_groups)] # stackup groups ...
        ax_cbar[idx] = fig.add_subplot(gs[-1,idx])
    for idx in extra_order:
        ax_xtra[idx] = [fig.add_subplot(gs[_i+1,idx]) for _i in range(num_stackup_groups)] # stackup groups ...

    # turning some of the input parameters into "oredered" or labeled dicts ...
    if norms is None:
        norms = { _i:None for _i in hmss_order}
    else:
        norms = { _i:norms[i] for i,_i in enumerate(hmss_order)}
    vlims = { _i:limss[i] for i,_i in enumerate(hmss_order)}
    titles = { _i:titles[i] for i,_i in enumerate(hmss_order)}
    if binsizes is None:
        binsizes = { _i:1 for _i in hmss_order}
    else:
        binsizes = { _i:binsizes[i] for i,_i in enumerate(hmss_order)}

    hm_arr = {}
    profile_hm = {}
    # for each group of stackups (vertically set)
    for group_id, k in enumerate(hmss):
        hm_arr[group_id] = {}
        profile_hm[group_id] = {}
        # for every stackup in each group (horizontal set)
        for idx, hm in zip(hmss_order, hmss[k]):
            if fillmissing:
                X = hm[:]
                missing = ~np.isfinite(X)
                mu = np.nanmean(X, axis=0, keepdims=True) # axis 0 or 1 - rows or columns ?!
                hm_arr[group_id][idx] = np.where(missing, mu, X)
            else:
                hm_arr[group_id][idx] = hm[:]
            if norms[idx] is None:
                profile_hm[group_id][idx] = np.nanmean(hm_arr[group_id][idx],axis=0)
            else:
                profile_hm[group_id][idx] = np.exp(np.nanmean(np.log(hm_arr[group_id][idx]),axis=0))

    for idx, cmap in zip(hmss_order, cmps):
        # plot profiles from every group on a single common axis for profiles...
        for _i in range(num_stackup_groups):
            ax_profile[idx].plot(profile_hm[_i][idx],linewidth=4)
        ax_profile[idx].set_yscale("linear" if norms[idx] is None else "log")
        # stackups for every group ...
        for _i in range(num_stackup_groups):
            stack_hm = ax_stackup[idx][_i].imshow(
                              hm_arr[_i][idx],
                              norm=norms[idx],
                              aspect="auto",
                              vmin=vlims[idx][0] if norms[idx] is None else None,
                              vmax=vlims[idx][1] if norms[idx] is None else None,
                              cmap=cmap,
                              interpolation=interpolation,
            )
        # beautify ...
        group_id_beautify = 0
        first_bin = 0-.5
        center_bin = hm_arr[group_id_beautify][idx].shape[1]/2 - .5
        last_bin = hm_arr[group_id_beautify][idx].shape[1]-.5
        ax_profile[idx].set_xlim([first_bin, last_bin])
        ax_profile[idx].set_ylim(vlims[idx])
        ax_profile[idx].set_title(titles[idx])
        ax_profile[idx].tick_params(axis="y", length=0, direction="in", pad=-5)
        ax_profile[idx].tick_params(axis="x", length=6)
        ax_profile[idx].set_yticks(vlims[idx])
        ax_profile[idx].set_yticklabels(vlims[idx],fontsize=fig_fontsize)
        for _tidx, tick in enumerate(ax_profile[idx].yaxis.get_majorticklabels()):
            tick.set_horizontalalignment("left")
            if _tidx == 0:
                tick.set_verticalalignment("bottom")
            elif _tidx == 1:
                tick.set_verticalalignment("top")
        # human readable kb stuff:
        flank_in_kb = int((center_bin+.5)*binsizes[idx]/1000)
        flank_ticks = [first_bin, center_bin, last_bin]
        flank_ticklabels = [-flank_in_kb, "", flank_in_kb]
        ax_profile[idx].set_xticks(flank_ticks)
        ax_profile[idx].set_xticklabels(flank_ticklabels,fontsize=fig_fontsize)
        for _tidx, tick in enumerate(ax_profile[idx].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 2:
                tick.set_horizontalalignment("right")
            else:
                tick.set_horizontalalignment("center")
        for _i in range(num_stackup_groups-1):
            ax_stackup[idx][_i].set_xticks([])
            ax_stackup[idx][_i].set_xticklabels([])
            ax_stackup[idx][_i].set_yticks([])
            ax_stackup[idx][_i].set_yticklabels([])
        # bottom one - show ticks for now ...
        _i = num_stackup_groups-1
        ax_stackup[idx][_i].set_xticks(flank_ticks)
        ax_stackup[idx][_i].set_xticklabels(flank_ticklabels,fontsize=fig_fontsize)
        ax_stackup[idx][_i].tick_params(axis="x", length=6)        
        ax_stackup[idx][_i].set_yticks([])
        ax_stackup[idx][_i].set_yticklabels([])
        for _tidx, tick in enumerate(ax_stackup[idx][_i].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 2:
                tick.set_horizontalalignment("right")
            else:
                tick.set_horizontalalignment("center")
        plt.colorbar(stack_hm,
                    cax=ax_cbar[idx],
                    orientation="horizontal",
                    ticks=vlims[idx])
        ax_cbar[idx].tick_params(axis="x", length=6)
        ax_cbar[idx].set_xticklabels(vlims[idx],fontsize=fig_fontsize)
        for _tidx, tick in enumerate(ax_cbar[idx].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 1:
                tick.set_horizontalalignment("right")

    return ax_xtra

def drop_5above_below(f1,f2,org):
    #This function is requires 2 files as input, since it checks for weights of
    # 5 bins above and below (10 bins) of a given bin/region, and if any of the two files have
    # nans in 10bins , that region is dropped
    d1 = pd.read_csv("weights/"+f1+".10kbweights.tsv",sep='\t')
    d2 = pd.read_csv("weights/"+f2+".10kbweights.tsv",sep='\t')
    to_drop = []
    for idx, row in org.iterrows():
        t1 = d1[(d1['chrom']==row['chrom']) & (d1['start']==row['start']) & \
                    (d1['end']==row['end'])].index[0]
        s1 = d1.iloc[t1-4:t1+6]['weight'].isna().sum()
        t2 = d2[(d2['chrom']==row['chrom']) & (d2['start']==row['start']) & \
                    (d2['end']==row['end'])].index[0]
        s2 = d2.iloc[t2-4:t2+6]['weight'].isna().sum()
        #if either of d1 or d2 has nans then drop that entry in both
        if(s1!=0 or s2!=0):
            to_drop.append(idx)
        if(idx%1000==0):
            print(idx)
    return to_drop

def weight_filtering(f1,f2,s1,s2):
    # This function takes 2 insulations files and merges them.
    # if any of bins in either of these have NaN weights in 5 bins around the considered bin, that region is dropped from the merged df

    d1 = pd.read_csv("cooltools54/Insulation/10kb_bins_100kb_win/"+f1+"_10kb_bin_100kb_win_insulation",\
                    sep='\t').drop(['is_bad_bin', 'n_valid_pixels_100000'], axis=1)
    d2 = pd.read_csv("cooltools54/Insulation/10kb_bins_100kb_win/"+f2+"_10kb_bin_100kb_win_insulation",\
                    sep='\t').drop(['is_bad_bin', 'n_valid_pixels_100000'], axis=1)
    m=d1.merge(d2, how='inner', on=['chrom','start','end'],suffixes=('_'+s1, '_'+s2))


    wt1=pd.read_csv("weights/"+f1+".10kbweights.tsv",sep='\t')
    wt2=pd.read_csv("weights/"+f2+".10kbweights.tsv",sep='\t')
    wt=wt1.merge(wt2, how='inner', on=['chrom','start','end'],suffixes=('_'+s1, '_'+s2))

    fm = m.merge(wt, how='inner', on=['chrom','start','end'])
    fm = fm[['chrom','start','end','weight_'+s1,'weight_'+s2,'log2_insulation_score_100000_'+s1,\
            'log2_insulation_score_100000_'+s2,'boundary_strength_100000_'+s1,'boundary_strength_100000_'+s2]]
    fm['diff-insul']=fm['log2_insulation_score_100000_'+s2]-fm['log2_insulation_score_100000_'+s1]
    fm['diff-strength']=fm['boundary_strength_100000_'+s2]-fm['boundary_strength_100000_'+s1]
    print("Entire df length:",len(fm))
    #fm.head(500).tail()

    #First Filtering
    rs=fm.dropna(subset=['weight_'+s1,'weight_'+s2])
    print(len(rs))
    discarded = fm[~fm.index.isin(rs.index)]
    #second filtering
    rs=rs.reset_index(drop=True)
    discarded = discarded.reset_index(drop=True)
    print("Dropping starts")
    to_drop = drop_5above_below(f1,f2,rs)
    final = rs.drop(to_drop).sort_values('diff-insul',ascending=False,key=abs).reset_index(drop=True)
    dropped = pd.concat([discarded, rs.loc[to_drop, :]], ignore_index=True).sort_values('diff-insul', ascending=False, key=abs)
    final.to_csv("Filtered-insul/"+f1+"_"+f2+"_weight_filtered.tsv",sep='\t',index=False,header=True)
    dropped.to_csv("Filtered-insul/"+f1+"_"+f2+"_weight_filtered_dropped.tsv",sep='\t',index=False,header=True)
    assert len(fm) == len(dropped)+len(final)