import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statistics import mode
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from matplotlib import ticker
import rstr
import matplotlib.gridspec as gridspec

def plot_rank_relation (df, perdictor1, perdictor2, threshold, kinase, dir):
    """
    plot scatter plot and histogram of the same-kinase Jaccard's Index rank between the two compared predictors
    
    Parameters 
    ---------
    df: dataframe
        rank summary
    perdictor1, perdictor2: str
        the two compared prodictor, also as the column names in the given df
    threshold: str
        stringency
    kinase: str
        kinase type
    dir: str
        output directory
      
    """
    if perdictor1 == 'pp':
        x_name = 'PhosphoPICK'
    elif perdictor1 == 'nw':
        x_name = 'NetworKIN'
    elif perdictor1 == 'gps':
        x_name = 'GPS'
        
    if perdictor2 == 'pp':
        y_name = 'PhosphoPICK'
    elif perdictor2 == 'nw':
        y_name = 'NetworKIN'
    elif perdictor2 == 'gps':
        y_name = 'GPS'
    
    # plot the rank relationship of the compared predictors in rank percentile
    df['percentile in '+ x_name] =  (df['Rank in ' + x_name]/df['number of kinase in ' + x_name])*100
    df['percentile in '+ y_name] =  (df['Rank in ' + y_name]/df['number of kinase in ' + y_name])*100

    x = df['percentile in ' + y_name].tolist()
    y = df['percentile in ' + x_name].tolist()

    thresh_x =  (df.at[0,'number of overlapped kinase']/df.at[0,'number of kinase in ' + y_name])*100
    thresh_y =  (df.at[0,'number of overlapped kinase']/df.at[0,'number of kinase in ' + x_name])*100
    
    sns.set(context='paper', style='white', palette='muted', font_scale=1)
     # plot histogram 
    plt.hist([x, y], bins = 10, label=['ref kinase: ' + x_name + ', rank in ' + y_name, 'ref kinase: ' + y_name + ', rank in ' + x_name])
    plt.legend(loc='upper left')
    plt.suptitle('Overlapped Kinases Similarity Rank between ' + x_name + ' and ' + y_name)
    plt.title('Threshold: ' + threshold + ', ' + kinase, fontsize=10)
    plt.xlabel('Top Percentile')
    plt.ylabel('Count')
    plt.xlim(1,100)
    plt.axvline(x=thresh_x, color='C0', linestyle='dashed', linewidth=1)
    plt.axvline(x=thresh_y, color='C1', linestyle='dashed', linewidth=1)
    plt.savefig(dir + x_name + '_' + y_name  + '_' + threshold + '_' + kinase + '_hist.svg', format='svg', bbox_inches="tight")
    plt.show()
    
    # plot the scatter plot
    df.plot(kind='scatter', x = 'percentile in ' + y_name , y = 'percentile in ' + x_name , alpha=0.3)
    plt.suptitle('Overlapped Kinases Similarity Rank between ' + x_name + ' and ' + y_name)
    plt.title('Threshold: ' + threshold + ', ' + kinase, fontsize=10)
    plt.xlim(0,100)
    plt.ylim(0,100)
    plt.axhline(y=thresh_y, color='C1', linestyle='dashed', linewidth=1)
    plt.axvline(x=thresh_x, color='C0', linestyle='dashed', linewidth=1)
    plt.xlabel('Top Percentile \n (' + 'ref kinase: ' + x_name + ', rank in: ' + y_name + ')')
    plt.ylabel('Top Percentile \n (' + 'ref kinase: ' + y_name + ', rank in: ' + x_name + ')')
    plt.savefig(dir + x_name + '_' + y_name  + '_' + threshold + '_' + kinase + '_scatter.svg', format='svg', bbox_inches="tight")
    plt.show()

def plot_rank_density (dfs, dir):
    """
    plot multiple CDFs in as subplots of same-kinase rank in the two compared predictors
    
    Parameters 
    ---------
    dfs: list of dataframes of  rank summaries
        [dfs[Y kin: nw vs pp], dfs[ST kin: nw vs pp],
         dfs[Y kin: nw vs gps], dfs[ST kin: nw vs gps],
         dfs[Y kin: pp vs gps], dfs[ST kin: pp vs gps]]
        
    dir: str
        output directory
        
    """
    pp = 'PhosphoPICK'
    nw = 'NetworKIN'    
    gps = 'GPS'
    
    y_xlim = 10
    st_xlim = 20

    sns.set(context='paper', style='white', palette='muted', font_scale=2.5)
    sns.set_style("ticks")
    plt.figure(0, figsize=(20, 15))

    if len(dfs) == 3:
        output = dir + 'density_sameStringency.svg'
    else:
        output = dir + 'density_mixedStringency.svg'
    
    # NetworKIN vs. GPS
    # Y kinases
    # rank of a nw kinase in gps predictions
    cumulative_df, label = get_cumelative_count(dfs[0], 'Rank in ' + gps)
    set_up_subplot ((0,0), cumulative_df, label, nw, gps, y_xlim)
    # rank of a gps kinase in nw predictions
    cumulative_df, label = get_cumelative_count(dfs[0], 'Rank in ' + nw)
    set_up_subplot ((0,1), cumulative_df, label, gps, nw, y_xlim)
    # ST kinases
    # rank of a nw kinase in gps predictions
    cumulative_df, label = get_cumelative_count(dfs[1], 'Rank in ' + gps)
    set_up_subplot ((0,2), cumulative_df, label, nw, gps, st_xlim)
    # rank of a gps kinase in nw predictions
    cumulative_df, label = get_cumelative_count(dfs[1], 'Rank in ' + nw)
    set_up_subplot ((0,3), cumulative_df, label, gps, nw, st_xlim)

    # PhosphoPICK vs. NetworKIN
    # Y kinases
    # rank of a nw kinase in pp predictions
    cumulative_df, label = get_cumelative_count(dfs[2], 'Rank in ' + pp)
    set_up_subplot ((1,0), cumulative_df, label, nw, pp, y_xlim)
    # rank of a pp kinase in nw predictions
    cumulative_df, label = get_cumelative_count(dfs[2], 'Rank in ' + nw)
    set_up_subplot ((1,1), cumulative_df, label, pp, nw, y_xlim)
    # ST kinases
    # rank of a nw kinase in pp predictions
    cumulative_df, label = get_cumelative_count(dfs[3], 'Rank in ' + pp)
    set_up_subplot ((1,2), cumulative_df, label, nw, pp, st_xlim)
    # rank of a pp kinase in nw predictions
    cumulative_df, label = get_cumelative_count(dfs[3], 'Rank in ' + nw)
    set_up_subplot ((1,3), cumulative_df, label, pp, nw, st_xlim)
    
    # PhosphoPICK vs. GPS
    # Y kinases
    # rank of a pp kinase in gps predictions
    cumulative_df, label = get_cumelative_count(dfs[4], 'Rank in ' + gps)
    set_up_subplot ((2,0), cumulative_df, label, pp, gps, y_xlim)
    # rank of a gps kinase in pp predictions
    cumulative_df, label = get_cumelative_count(dfs[4], 'Rank in ' + pp)
    set_up_subplot ((2,1), cumulative_df, label, gps, pp, y_xlim)
    # ST kinases
    # rank of a pp kinase in gps predictions
    cumulative_df, label = get_cumelative_count(dfs[5], 'Rank in ' + gps)
    set_up_subplot ((2,2), cumulative_df, label, pp, gps, st_xlim)
    # rank of a gps kinase in pp predictions
    cumulative_df, label = get_cumelative_count(dfs[5], 'Rank in ' + pp)
    set_up_subplot ((2,3), cumulative_df, label, gps, pp, st_xlim)

#     plt.subplots_adjust(wspace=0.05, hspace=0.7)
    plt.subplots_adjust(wspace=0.05, hspace=0.3)
    plt.suptitle('    '+rstr.xeger('\|-{18}Tyrosine Kinases-{19}\|'+' '+'\|-{12}Serine/Threonine  Kinases-{13}\|'), 
                 fontsize=27)
    
    plt.savefig(output, format='svg')
    plt.show()


def get_cumelative_count(dfs, ref_p):
    """
    get the cumelative fraction kinases of the same-kinase rank
    
    Parameters 
    ---------
    dfs: list of dataframes of  rank summaries
        [df_low, df_med, df_high] 
        or 
        [df_low_low, df_low_med, df_low_high,
         df_med_low, df_med_med,df_med_high,
         df_high_low, df_high_med, df_high_high]
         
    ref_p: str
        reference predictor
    
    Returns
    -------
    cumulative_count_dfs: list of dataframes of the Rank and cumelative fraction kinases
    label: list
        label setting for plots
        
    """
    
    cumulative_count_dfs = []
    color = ['c','m','y','k']
    
    for df in dfs:
        index = 0
        for i in range(1, df.at[0,'number of overlapped kinase'] +1):
            if i == 1:
                cumulative_count = pd.DataFrame(columns=['Rank', 'cumulativeCounts'])
                cumulative_count.at[index,'Rank'] = i
                count = df[df[ref_p]==i][ref_p].count()
                cumulative_count.at[index,'cumulativeCounts'] = count/df.at[0,'number of overlapped kinase']
            else:
                cumulative_count.at[index,'Rank'] = i
                count = df[df[ref_p]==i][ref_p].count()
                cumulative_count.at[index,'cumulativeCounts'] = count/df.at[0,'number of overlapped kinase'] + cumulative_count.at[index-1,'cumulativeCounts']
         
            index +=1

        cumulative_count_dfs.append(cumulative_count)
    
    random_rank = pd.DataFrame(columns=['Rank', 'cumulativeCounts'])
    for j in range(0, df.at[0,'number of overlapped kinase']):
        random_rank.at[j,'Rank'] = j+1
        random_rank.at[j,'cumulativeCounts'] = (j+1)/df.at[0,'number of overlapped kinase']       
    cumulative_count_dfs.append(random_rank)        
       
    num_of_dfs = len(dfs)  
        
    if  num_of_dfs == 3:
        label = [['low_low', color[0], 1], 
                 ['med_med', color[1], 1],
                 ['high_high', color[2], 1],
                 ['random_chance', color[3], 1]]
        
    elif num_of_dfs == 9:
        label = [['low_low',color[0],0.2], ['low_med',color[0],0.6], ['low_high',color[0],1],
                 ['med_low',color[1],0.2], ['med_med',color[1],0.6], ['med_high',color[1],1],
                 ['high_low',color[2],0.2], ['high_med',color[2],0.6], ['high_high',color[2],1],
                 ['random_chance', color[3], 0.8]]
        
    return cumulative_count_dfs, label

def set_up_subplot (position, cumulative_df, label, p_kin, p_ref, xlim):
    """
    setup and plot subplot
    
    Parameters 
    ---------
    position: position of the subplot in the grid
    cumulative_df: list of dataframes of the Rank and cumelative fraction kinases
        [df_low, df_med, df_high] 
        or 
        [df_low_low, df_low_med, df_low_high,
         df_med_low, df_med_med,df_med_high,
         df_high_low, df_high_med, df_high_high]
    label: list
        label setting for plots [label, color, alpha]
         
    p-kin: str
        prodictor of the ranked kinase
    p_ref: str
        prodictor of the kinase ranked in 
    xlim: int
        
    """
    
    # setup the grid for subplots
    ax = plt.subplot2grid((3,4), position)
    ax.set_title(p_kin + ' : ' + p_ref, fontsize=22)
    
    i = 0
    # plot multiple cdfs in one subplot
    for j in cumulative_df:
        plt.plot(j['Rank'], j['cumulativeCounts'], marker='o', markersize=3, linewidth=1.5, label = label[i][0], color = label[i][1], alpha = label[i][2])
        i +=1
    # set yTicks
    n_yTick = 10
    ax.yaxis.set_major_locator(ticker.MaxNLocator(n_yTick))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    # only show x axis title for the bottom row subplots
    if position == (2,0) or position == (2,1) or position == (2,2) or position == (2,3):
        plt.xlabel('Rank')
    # only show y axis title for the left most subplots
    if position == (0,0) or position == (1,0) or position == (2,0):
        plt.ylabel('Fraction Kinases')
    else:
        ax.yaxis.set_major_formatter(plt.NullFormatter()) 
    # only show legends for the right most subplots
    if position == (0,3) or position == (1,3) or position == (2,3):
        # set costomized legends
        if len(cumulative_df) == 4:
            leg_line = [Line2D([0], [0], color='c', lw=1.5, marker='o', markersize=3),
                        Line2D([0], [0], color='m', lw=1.5, marker='o', markersize=3),
                        Line2D([0], [0], color='y', lw=1.5, marker='o', markersize=3),
                        Line2D([0], [0], color='k', lw=1.5, marker='o', markersize=3)]
            leg = ax.legend(leg_line, ['low', 'med' , 'high', 'random \n chance'], fontsize=18, loc='upper left', bbox_to_anchor=(1, 1.1), handletextpad=0.1, frameon=False)
        else:
            leg_line1 = [Line2D([0], [0], color='c', lw=1.5, marker='o', markersize=3),
                        Line2D([0], [0], color='m', lw=1.5, marker='o', markersize=3),
                        Line2D([0], [0], color='y', lw=1.5, marker='o', markersize=3)]
            leg1 = ax.legend(leg_line1, ['low', 'med' , 'high'], fontsize=18, loc='upper left', bbox_to_anchor=(1, 1.1), handletextpad=0.1, frameon=False)
            leg1.set_title(p_ref, prop = {'size':17})
            leg1._legend_box.align = 'left'
            leg_line2 = [Line2D([0], [0], color='gray', lw=1.5, marker='o', markersize=3, alpha = 0.2),
                        Line2D([0], [0], color='gray', lw=1.5, marker='o', markersize=3, alpha = 0.6),
                        Line2D([0], [0], color='gray', lw=1.5, marker='o', markersize=3, alpha = 1), 
                        Line2D([0], [0], color='k', lw=1.5, marker='o', markersize=3)]
            leg2 = ax.legend(leg_line2, ['low', 'med' , 'high', 'random \n chance'], fontsize=18, loc='upper left', bbox_to_anchor=(1, 0.65), handletextpad=0.1, frameon=False)
            leg2.set_title(p_kin, prop = {'size':17})
            leg2._legend_box.align = 'left'
            ax.add_artist(leg1)
    # x-axis label should be int only
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # remove the 1st x-axis tick lable for the subplots in the 2nd column from the left. (avoid tick label overlap)
    if position == (0,1) or position == (1,1) or position == (2,1):
        xticks = ax.xaxis.get_major_ticks()
        xticks[0].set_visible(False)
    plt.xlim(1,xlim)
    plt.ylim(0,0.9)
