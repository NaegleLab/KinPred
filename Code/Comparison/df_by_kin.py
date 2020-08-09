import pandas as pd

def df_by_kin (df):
    """
    divided the given df by kinase type
    return two dfs only contain either Y kinases or S/T kinases

    Parameters
    ----------
        df: dataframe contains all kinases
    Returns:
        df_y:  Dataframe contains y kinases
        df_st: Dataframe contains y kinases
    """
    # file location (../../ for local, /Volumes/naegle_lab/Kinase Predictions/ for sammas)
    base = '/Volumes/naegle_lab/Kinase Predictions/'
    # globalKinaseMap location
    KinaseMap = base + 'Data/Map/globalKinaseMap.csv'
    
    # Kinase name type: Kinase Name or Preferred Name
    kin_name = 'Preferred Name' # use Preferred Name in all analysis
    kin = pd.read_csv(KinaseMap, usecols = [kin_name, 'Type'])

    # group kinases by kinase type in globalKinaseMap (annotation based on UniProt)
    y_kin = kin.loc[kin['Type'] == 'Pkinase_tyr']              # Y kinases
    st_kin = kin.loc[kin['Type'] == 'Pkinase']                 # S/T kinase
    dual_kin = kin.loc[kin['Type'] == 'Pkinase/Pkinase_tyr']   # dual specificity kinases

    df_y = df[(df['Kinase Name'].isin(y_kin[kin_name].append(dual_kin[kin_name])))] 
    # only keep dual specificity kinases at pY site
    df_y = df_y[df_y['site'].str.contains('Y')]
    
    df_st = df[df['Kinase Name'].isin(st_kin[kin_name].append(dual_kin[kin_name]))]
    # only keep dual specificity kinases at pS/T site
    df_st = df_st[df_st['site'].str.contains('S|T')]
    
    return df_y, df_st