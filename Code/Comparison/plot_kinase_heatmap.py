import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_kinase_heatmap(heatmap, filename, idx, col, kinase):
    """
    plot and save the heatmap of given Jaccard's Index matrix
    
    Parameters 
    ----------
    heatmap: dataframe of Jaccard's Index matrix
    filename: str
        output filename
    idx: str
        predictor name of the kinases represented by the indexes
    col: str
        predictor name of the kinases represented by the columns
    kinase: str
        kinase type
    
    """
    with sns.axes_style("white"):
        if kinase == 'y':
            sns.set(font_scale = 3)
            y_labels = True # show all kinase labels
        elif kinase == 'st':
            sns.set(font_scale = 2)
            y_labels, x_labels = reset_label (heatmap, 's') # only show kinases fit the criteria 

        if -1 in heatmap.values: #between predictors 
            # plot all kinases predicted by either of the predictor
            # add titles for the indexes and columns
            heatmap.index.name = idx
            heatmap.columns.name = col
            # plot heatmap with sorted Jaccard's Index matrix
            sort_heatmap (heatmap, kinase)
            plt.savefig(filename + '.eps', format='eps', bbox_inches="tight")
            if kinase == 'y':
                sort_heatmap (heatmap, 'st') # print another heatmap only show kinases that fit the criteria 
            
            # plot overlapped kinases predicted by both of the predictor
            masked_heatmap = heatmap.loc[~np.all(heatmap == -1, axis=1), ~np.all(heatmap == -1, axis=0)]
            masked_heatmap.index.name = idx
            masked_heatmap.columns.name = col
            sort_heatmap (masked_heatmap,  kinase)
            plt.savefig(filename + '_overlap.eps', format='eps', bbox_inches="tight")
            if kinase == 'y':
                sort_heatmap (masked_heatmap, 'st')
        else: # within predictor
            ax = sns.clustermap(heatmap, 
                            xticklabels = False, yticklabels = y_labels, 
                            cbar_kws={"orientation": "horizontal"}, cbar_pos=(0.1975, -0.05, 0.7175, 0.025), 
                            tree_kws=dict(linewidths=2), figsize = (50,50), vmin=0, vmax=1)
            ax.ax_col_dendrogram.set_visible(False)
            plt.savefig(filename + '.eps', format='eps', bbox_inches="tight")


def sort_heatmap (df, kinase):
    """
    sort the kinases in decreasing order of the same-kinase Jaccard's Index
    plot the heatmap in the sorted order
    
    Parameters 
    ----------
    df: dataframe
        Jaccard's Index matrix
    kinase: str
        kinase type
    
    Returns
    -------
    ax: ploted heatmap of given Jaccard's Index matrix
    """
    # sort the given Jaccard's Index matrix
    order = df.stack().sort_values(ascending=False).reset_index()
    order = order[order.iloc[:,0] == order.iloc[:,1]].reset_index(drop=True)
    order = order.iloc[:,0].tolist()
    df = df.reindex(order)
    df = df.reindex(columns = order)
    
    if kinase == 'y': 
        x_label = True
        y_label = True
    elif kinase == 'st':
        y_label, x_label = reset_label (df, 'p') # only show kinases fit the criteria 
    
    # plot the ordered Jaccard's Index matrix
    plt.figure(figsize = (30,30))
    ax = sns.heatmap(df, mask = df == -1, 
                     xticklabels = x_label, yticklabels = y_label, 
                     square = True, vmin = 0, vmax = 1)
    ax.tick_params(left=True, bottom=True)
    return ax

def reset_label (df, df_type):
    """
    for S/T kinases, only show label if the kinases:
        within predictor: have at lease one Jaccard's Index, other than the comparison with itself, higher than 0.49 
        between predictors: have at lease one Jaccard's Index is at the top 1%
    
    Parameters 
    ----------
    df: dataframe 
    df_type: str
        comparison type. 'p' for between predictors (paired); 's' for within predictor (single)
    
    Returns
    -------
    new_idx: list of labels that would show on the y-axis
    new_col: list of labels that would show on the x-axis
    """
    
    
    new_idx = []
    new_col = []
    if df_type == 'p': #between predictors (paired)
        for i in df.columns:
            if ((np.percentile(df[i], 99) <= df.at[i,i]) & (df[i] < 1)).any():
                new_col.append(i)
            else:
                new_col.append('')
        for j in df.index:
            if ((np.percentile(df.loc[j], 99) <= df.at[j,j]) & (df.loc[j] < 1)).any():
                new_idx.append(j)
            else:
                new_idx.append('')
    elif df_type == 's': # within predictor (single)
        thresh = 0.4900
        for i in df.columns:
            if ((df[i] >= thresh) & (df[i] < 1)).any():
                new_col.append(i)
            else:
                new_col.append('')
        for j in df.index:
            if ((df.loc[j] >= thresh) & (df.loc[j] < 1)).any():
                new_idx.append(j)
            else:
                new_idx.append('')

    return new_idx, new_col
        
