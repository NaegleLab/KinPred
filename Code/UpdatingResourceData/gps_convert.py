
import pandas as pd
import re
import os
import glob
#only need when testing the code
import time
import checkSite

def convert_acc(filename, valid_kinases):
    """
    Given a filename of raw GPS results the file is converted to a dataframe and returned 
    
    Parameters
    ----------
    filename : str
        gps raw file
    valid_kinases : pandas df
        gps predictor name
        
    Returns
    ----------
    gps: pandas df
        'substrate_acc'     substrate Uniprot id
        'residue'           S,T,Y
        'position'          position in sequence of phosphosite
        'site'               residue + position
        'Predictor'         full GPS predictor, e.g. TK/Src/SrcA/SRC
        'kinase'            name of the kinase in gps
        'kinase_acc'        uniprot id
        'pep'               +/- 7 peptide sequence surrounding phosphosite
        'score'             GPS score of phosphosite
        'threshold'         GPS cutoff threshold    
        
    """
    substrate_acc = ''
    # set up a dict for the df columns
    df_dict =  {
                'substrate_acc' : [],
                'residue'       : [],
                'position'      : [],
                'predictor'     : [],
                'pep'           : [],
                'score'         : [],
                'threshold'     : [],
                }
    with open(filename, 'r') as file:
        line = file.readline()
        while line:
            line = file.readline()
            line = line.rstrip('\n')
            # extract the 'substrate_acc' (substrate uniprotID) from the identifier line
            if line.startswith('>'):
                substrate_acc = line.split('|')[1]
            # extract other columns info from non-iditifier lines
            else:
                split = line.split('\t' or '\n')
                # add data to the dict
                if len(split) == 6:
                    df_dict['substrate_acc'].append(substrate_acc)
                    df_dict['position'].append(split[0])
                    df_dict['residue'].append(split[1])
                    df_dict['predictor'].append(split[2])
                    pep = split[3]
                    # conver the peptide format (PPGKHLVTEV***** --> PPGKHLVTEV_____)
                    pep = re.sub(r"\*", "_",pep)
                    df_dict['pep'].append(pep)
                    df_dict['score'].append(split[4])
                    df_dict['threshold'].append(split[5])
    # build the df using df_dict        
    df = pd.DataFrame(df_dict)
    # create the 'site' column (aa + position in protein seq)
    df['site'] = df['residue'] + df['position'].map(str)
    # add 'kinase_acc' (kinase uniprotID) by merge the valid_kinases file to the df
    df = pd.merge(df, valid_kinases, how = 'inner', on = 'predictor')
    df = df[['substrate_acc', 'site', 'residue', 'position', 'predictor','kinase', 'kinase_acc', 'pep', 'score', 'threshold']]
    
    return df

def map_site (filename,  human_proteome):
    """
    Given gps dataframe file, positionally map the site to the given human proteome sequence dataframe
    
    Parameters
    ----------
    filename : str
        gps raw file
    human_proteome : str
        given human proteome sequence dataframe

        
    Returns
    ----------
    gps: pandas df
        'substrate_id'      substrate_acc + position 
        'substrate_acc'     substrate Uniprot id
        'residue'           S,T,Y
        'position'          position in sequence of phosphosite
        'site'              residue + position 
        'mapped site'       residue + position in provided human proteome
        'predictor'         full GPS predictor, e.g. TK/Src/SrcA/SRC
        'kinase'            name of the kinase in gps
        'kinase_acc'        uniprot id
        'pep'               +/- 7 peptide sequence surrounding phosphosite
        'score'             GPS score of phosphosite
        'threshold'         GPS cutoff threshold 
        
    """
    start = time.time()
    print ('Reading input file...')
    df = pd.read_csv(filename)
    
    print ('Get unique substrate sites...')
    # get a list of unique substrate/site in GPS
    df_unique_sub = df[['substrate_acc', 'position', 'pep']]
    df_unique_sub = df_unique_sub.drop_duplicates()
    print ('Map unique substrate sites...')
    for index, row in df_unique_sub.iterrows():
        id = df_unique_sub.at[index, 'substrate_acc']
        pos = str(df_unique_sub.at[index, 'position'])
        pep = df_unique_sub.at[index, 'pep']
        pos_in_pep = 8
        # map the site by pep search
        new_site, site_confirm = checkSite.checkSite(id, pos, pep, pos_in_pep, human_proteome)
        
        if site_confirm == True:
            df_unique_sub.at[index,'mapped site'] = new_site
            # assign 'substrate_id' (substrate_acc + pos in proteion seq)
            df_unique_sub.at[index,'substrate_id'] = id + '_' + str(re.search(r'(\d+)', new_site).group(1))
        else:
            df_unique_sub.at[index,'mapped site'] = new_site
            df_unique_sub.at[index,'substrate_id'] = 'outdated'

    # add the 'substrate_id' to GPS df (df_subMap)
    df = df.merge(df_unique_sub, left_on=['substrate_acc', 'position', 'pep'], right_on=['substrate_acc', 'position', 'pep'], how = 'left')
            
    return df

def gps_convert_directory(load_directory, reference_filename, save_directory, convert_type):
    """
    Pulls all files from loading directory and save the output dataframe to the given saving directory
    
    Parameters
    ----------
    load_directory : str
        location to pull files from
    reference_filename : str
        filename of valid_kinase csv file
        OR, filename of the referece human proteome sequence dataframe
    save_directory : str
        location to save the GPS result dataframe
    convert_type : str
        acc (convert substrate/kinase accession)
        OR, site (mapping site to reference sequence)
    """
    # for converting substrate and kinase accession
    if convert_type == 'acc':
        valid_kinases = pd.read_csv(reference_filename)

        all_files = glob.glob(load_directory + "*.gps.fasta")
        for filename in all_files:
            start = time.time()
            name = re.search(r'.+\/(.+)\.gps\.fasta$', filename).group(1)
            print ("Formatting ", name, "...")
            df = convert_acc(filename, valid_kinases)
   
            df.to_csv(save_directory + name + '.csv', index=False)
            
            end = time.time()
            print (f"Done. Time\t{(end-start):.3f}")
    # for mapping the site to the reference human proteome seq
    elif convert_type == 'site':
        HP_df = pd.read_csv(reference_filename, sep = '\t')

        all_files = glob.glob(load_directory + "*.csv")
        for filename in all_files:
            start = time.time()
            name = re.search(r'.+\/(.+)$', filename).group(1)
            print ("Formatting ", name, "...")
            df = map_site(filename, HP_df)
   
            df.to_csv(save_directory + name, index=False)
            
            end = time.time()
            print (f"Done. Time\t{(end-start):.3f}")
