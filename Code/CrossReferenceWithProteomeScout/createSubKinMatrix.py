import pandas as pd
import os
import time

def createMatrix(fileName, output):
    """
    create and save the matrix from given dataframe
    
    Parameters
    ----------
    fileName : str
        input file name
    output : str
        output file name
        
    Returms
    --------
    df_matrix: dataframe
    
    """

    start = time.time()
    print ('reading input file....')
    df_input = pd.read_csv(fileName)
    end = time.time()
    print ('time: ', end - start)
    start = time.time()
    print ('creating matrix....')

    df_input['substrate_name'] = df_input['substrate_name'].fillna('-')
    # for prediction data that have scores for the predictions
    if 'score' in df_input.columns:
        df_matrix = df_input.pivot_table(index = ['substrate_id', 'substrate_name','substrate_acc','site','pep' ], columns = 'Kinase Name', values = 'score', fill_value = '-')
    # for annotation data that do not have scores (binary)
    else:
        df_input['found'] = 1 # create a new column, this would be the value (binary) of the matrix
        df_matrix = df_input.pivot_table(index = ['substrate_id', 'substrate_name','substrate_acc','site','pep' ], columns = 'Kinase Name', values = 'found', fill_value = 0)
             
    df_matrix.reset_index(inplace=True)
    end = time.time()
    print ('time: ', end - start)
    start = time.time()
    print ('saving file....')
    df_matrix.to_csv(output, chunksize = 1000000,index=False)
    end = time.time()
    print ('time: ', end - start)
    
    return df_matrix