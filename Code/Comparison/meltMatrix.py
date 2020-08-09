import pandas as pd

def meltMatrix (file):
    """
    create a list-like df given matrix dataframe
    
    Parameters
    ----------
    file : str
        input file name
    Returns
    ----------
    df: a list-like dataframe
    """
    df = pd.read_csv(file)
    df = pd.melt(df, 
                id_vars=['substrate_id','substrate_name','substrate_acc','site','pep'], 
                value_vars=list(df.columns[5:]),
                var_name='Kinase Name', 
                value_name='score')
    # remove kinase-substrate pairs with no score (not predicted)
    df = df[df['score']!='-']
    return df
