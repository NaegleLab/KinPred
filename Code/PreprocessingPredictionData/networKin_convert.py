
import pandas as pd
import re
import os
import glob
#only need when testing the code
import time
import checkSite, getUniprotID

def convert_acc (filename):
    """
    Given a filename of raw NetworKIN results the file is converted to a dataframe and returned 
    
    Parameters
    ----------
    filename : str
        NetworKIN raw file
        
    Returns
    ----------
    NetworKIN: pandas df 
                '#Name'                                                     input fasta file deq identifier
                'Tree'                                                      type of prediction (we only need Kinase prediction 'KIN')
                'Target description'                                        substrate gene name
                'Position'                                                  aa + position in protein sequence
                'substrate_acc'                                             substrate uniprotID
                'Kinase/Phosphatase/Phospho-binding domain'                 kinase name used in NetworKIN
                'Kinase/Phosphatase/Phospho-binding domain description'     kinase name used in NetworKIN        
                'kinase_acc'                                                kinase uniprotID
                'NetworKIN score'                                           NetworKIN score 
                'Peptide sequence window'                                   +/- 5 AA
        
    """

    # define networKIN result df as df_raw
    df_raw = pd.read_csv(filename, usecols = ['#Name', 'Tree', 'Target description', 'Position', 'NetworKIN score', 'Kinase/Phosphatase/Phospho-binding domain', 'Kinase/Phosphatase/Phospho-binding domain description', 'Peptide sequence window'], sep='\t')
    # Filtering data: only keep the results of kinase prediction 
    df_raw = df_raw[df_raw.Tree == 'KIN']
    
    print ('getting unique sub')
    df_unique_sub = df_raw[['#Name']]
    df_unique_sub = df_unique_sub.drop_duplicates()
    print ('getting sub_acc')
    # get the substrate_acc from the input sequence identifier 
    df_unique_sub['substrate_acc'] = df_unique_sub.apply(lambda row :  row['#Name'].split("|")[1], axis = 1) 
    print ('merge')
    # marge the 'substrate_acc' column to df_raw
    df_raw = df_raw.merge(df_unique_sub, left_on=['#Name'], right_on=['#Name'], how = 'left')

    # correct the name in 'Kinase/Phosphatase/Phospho-binding domain description'
    correct_kinase = {'MST4' : 'STK26'}
    for key in correct_kinase:
        df_raw.loc[df_raw['Kinase/Phosphatase/Phospho-binding domain'] == key, 'Kinase/Phosphatase/Phospho-binding domain description']= correct_kinase[key]
    
    print ('getting unique kin')    
    df_unique_kin = df_raw[['Kinase/Phosphatase/Phospho-binding domain', 'Kinase/Phosphatase/Phospho-binding domain description']]
    df_unique_kin = df_unique_kin.drop_duplicates()
    print ('getting kin_acc')
    # get uniprotID for the kinases using 'Kinase/Phosphatase/Phospho-binding domain description' column
    df_unique_kin['kinase_acc'] = df_unique_kin.apply(lambda row :  getUniprotID.getUniprotID(row['Kinase/Phosphatase/Phospho-binding domain description'], 'gene name'), axis = 1) 
    print ('merge')
    # marge the 'kinase_acc' column to df_raw
    df_raw = df_raw.merge(df_unique_kin, left_on=['Kinase/Phosphatase/Phospho-binding domain', 'Kinase/Phosphatase/Phospho-binding domain description'], right_on=['Kinase/Phosphatase/Phospho-binding domain', 'Kinase/Phosphatase/Phospho-binding domain description'], how = 'left')

    return df_raw

def map_site (filename,  human_proteome):
    """
    Given NetworKIN dataframe file, positionally map the site to the given human proteome sequence dataframe
    
    Parameters
    ----------
    filename : str
        NetworKIN file
    human_proteome : str
        given human proteome sequence dataframe

        
    Returns
    ----------
    NetworKIN: pandas df
                'substrate_id'      unique IDs for the substrate phosphorylation site (substrate_acc + position)
                'substrate_acc'     substrate uniprotID
                'substrate_name'    substrate gene name
                'Position'          aa + position in protein sequence (orginal)
                'site'              aa + position in protein sequence (mapped)
                'pep'               +/- 5 AA
                'kinase_name'       kinase name used in NetworKIN
                'kinase_acc'        kinase uniprotID
                'score'             NetworKIN score 
        
    """
    df_temp = pd.read_csv(filename, usecols = ['Target description','substrate_acc', 'Kinase/Phosphatase/Phospho-binding domain description', 'kinase_acc', 'Position','Peptide sequence window', 'NetworKIN score'])
    # rename columns
    df_temp = df_temp.rename(columns={'Target description' : 'substrate_name',
                                        'Kinase/Phosphatase/Phospho-binding domain description' : 'kinase_name',
                                        'Peptide sequence window' : 'pep', 
                                        'NetworKIN score' : 'score'})
    # formatting pep seq (e.g. ----MsGSKSV --> ____MSGSKSV)
    df_temp['pep'] = df_temp.apply(lambda row :  re.sub(r"\-", "_", row['pep']).upper(), axis = 1) 

    df_unique_site = df_temp[['substrate_acc', 'Position', 'pep']]
    df_unique_site = df_unique_site.drop_duplicates()

    for index, row in df_unique_site.iterrows():
        id = df_unique_site.at[index,'substrate_acc']
        site = df_unique_site.at[index,'Position']
        pep = df_unique_site.at[index,'pep'] 
        # NetworKIN provide pep seq of +/- 5 AA at the site, so the position in peptide should be 6
        pos_in_pep = 6 
        
        new_site, site_confirm = checkSite.checkSite(id, site, pep, pos_in_pep, human_proteome)
        # if the site mapped to the referece human proteome seq at the exact position, the new_site = site, site_confirm = True
        # if the peptide mapped to the referece human proteome seq with a shift in position, new_site = site + shift, site_confirm = True
        if site_confirm == True:
            df_unique_site.at[index,'site'] = new_site
            # substrate_id = substrate_acc + position in protein seq
            df_unique_site.at[index,'substrate_id'] = id + '_' + str(int(re.search(r'(\d+)', site).group(1))) 
        # if the peptide can not mapped to the referece human proteome seq, site_confirm = False
        else:
            df_unique_site.at[index,'site'] = new_site
            df_unique_site.at[index,'substrate_id'] = 'outdated'
    # add 'site' and 'substrate_id' to the df
    df_temp = df_temp.merge(df_unique_site, left_on=['substrate_acc', 'Position', 'pep'], right_on=['substrate_acc', 'Position', 'pep'], how = 'left')

    return df_temp

def kin_convert_directory(load_directory, reference_filename, save_directory, convert_type):
    """
    Pulls all files from loading directory and save the output dataframe to the given saving directory
    
    Parameters
    ----------
    load_directory : str
        location to pull files from
    reference_filename : str
        filename of the referece human proteome sequence dataframe
        OR, 'na'
    save_directory : str
        location to save the NetworKIN result dataframe
    convert_type : str
        acc (convert substrate/kinase accession)
        OR, site (mapping site to reference sequence)
    """
    if not os.path.exists(save_directory):
                os.mkdir(save_directory) 
    # for converting substrate and kinase accession
    if convert_type == 'acc':
        all_files = glob.glob(load_directory + '*.txt')
        for filename in all_files:
            start = time.time()
            name = re.search(r'.+\/(.+).txt$', filename).group(1)
            print ("Formatting ", name, "...")
            df = convert_acc(filename)
   
            df.to_csv(save_directory + name + '.csv', index=False)
            
            end = time.time()
            print (f"Done. Time\t{(end-start):.3f}")
    # for mapping the site to the reference human proteome seq
    elif convert_type == 'site':
        HP_df = pd.read_csv(reference_filename, sep = '\t')

        all_files = glob.glob(load_directory + '*.csv')
        for filename in all_files:
            start = time.time()
            name = re.search(r'.+\/(.+)\.csv$', filename).group(1)
            print ("Formatting ", name, "...")
            df = map_site(filename, HP_df)

            df.to_csv(save_directory + name + '_mappedSite.csv', index=False)
            
            end = time.time()
            print (f"Done. Time\t{(end-start):.3f}")
