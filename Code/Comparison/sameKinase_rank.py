import pandas as pd

def kinase_rank(heatmap, df1, df2, predictor_idx, predictor_col):
    """
    Sort the Jaccard's Index for each kinases from the two compared predictors (index, columns) in given heatmap.  
    
    Parameters
    ----------
        heatmap: dataframe 
            Jaccard's Index matrix
        df1, df2: dataframe 
            prediction of the two compared predictors
        predictor_idx: str
            predictor represent by the idx of the heatmap
        predictor_col: str
            predictor represent by the cols of the heatmap
        
    Returns
    --------
        kin_rank: dataframe
            'Kinase':                                    Kinase Name
            'number of overlapped kinase'
            'Rank in ' + predictor_col:                  predictor_idx kinase rank in predictor_col
            'Rank in ' + predictor_idx:                  predictor_col kinase rank in predictor_idx
            'number of kinase in '+ predictor_col
            'number of kinase in '+ predictor_idx

    """
    # get list of unique kinases from the compaired predictors
    kinases1 = df1['Kinase Name'].drop_duplicates().tolist()
    kinases2  = df2['Kinase Name'].drop_duplicates().tolist()
    # get the number of kinases predicted by the two compaired predictors
    num_idx_kin = len(kinases1)
    num_col_kin = len(kinases2)
    # get the overlapped kinases
    overlap_kin = list(set(kinases1) & set(kinases2))
    # get the number of overlapped kinases
    num_overlap_kin =  len(overlap_kin)
    # create a dataframe for the ranking summary
    kin_rank = pd.DataFrame(columns=[
                                        'Kinase', # Kinase Name
                                        'number of overlapped kinase', # of the two predictors 
                                        'Rank in ' + predictor_col, # predictor_idx kinase rank in predictor_col
                                        'Rank in ' + predictor_idx, # predictor_col kinase rank in predictor_idx
                                        'number of kinase in '+ predictor_col, 
                                        'number of kinase in '+ predictor_idx, 
                                        ]
                                )
    # fill the dataframe columns
    kin_rank['Kinase'] = overlap_kin
    kin_rank['number of overlapped kinase'] = num_overlap_kin
    kin_rank['number of kinase in '+ predictor_col] = num_col_kin
    kin_rank['number of kinase in '+ predictor_idx] = num_idx_kin
    # iterate through the heatmap indexes 
    for i in heatmap.index:
        # sort values of each row
        li = heatmap.loc[i , : ].sort_values(ascending=False)
        li = li.reset_index()
        li = li.rename(columns={'index': 'kin'})
        # if the kinase that represented by the index is present in both compared predictors
        if li[i].nunique() != 1:
            # iterate through the row,
            for index, row in li.iterrows():
                # find the same kinase 
                if li.at[index,'kin'] == i:
                    # get the position in the sorted list
                    rank = index+1
                    # fill the rank in the dataframe
                    kin_rank.loc[kin_rank.Kinase == i, ['Rank in ' + predictor_col]] = rank
    # iterate through the heatmap cols
    for j in heatmap.columns:
        # sort values of each col
        li = heatmap[j].sort_values(ascending=False)
        li = li.reset_index()
        li = li.rename(columns={'index': 'kin'})
        # if the kinase that represented by the col is present in both compared predictors
        if li[j].nunique() != 1:
            # iterate through the col
            for index, row in li.iterrows():
                # find the same kinase 
                if li.at[index,'kin'] == j:
                    # get the position in the sorted list
                    rank = index+1
                    # fill the rank in the dataframe
                    kin_rank.loc[kin_rank.Kinase == j, ['Rank in ' + predictor_idx]] = rank
    
    print (len(kin_rank.index), 'overlapped kinases \n', 
             predictor_idx, ': ', num_idx_kin, 'kinases','\n', 
             predictor_col, ': ',num_col_kin, 'kinases', '\n')
    
    return kin_rank
