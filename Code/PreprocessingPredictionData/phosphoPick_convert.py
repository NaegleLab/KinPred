
import pandas as pd
import re
import os
import glob
#only need when testing the code
import time
import checkSite, getUniprotID

def convert_acc (filename):
    """
    Given a filename of raw PhosphoPICK results the file is converted to a dataframe and returned 
    
    Parameters
    ----------
    filename : str
        PhosphoPICK raw file
        
    Returns
    ----------
    PhosphoPICK: pandas df 
                'Uniprot-Acc'       the uniprotID of the predicted substrate by PhosphoPICK
                'blastp-identity'   blastp-identity of the predicted substrate to the submitted seq
                'substrate_acc'     substrate uniprotID
                'position'          position in protein sequence
                'kinase'            kinase name used in PhosphoPICK
                'kinase_acc'        kinase uniprotID
                'combined-p-value'  PhosphoPICK score 
        
    """

    df_raw = pd.read_csv(filename, usecols = ['identifier','Uniprot-Acc','blastp-identity', 'kinase', 'site', 'combined-p-value'], sep='\t')
    df_raw = df_raw.rename(columns={'site': 'position'})
    # only keep results that have 100% seq match
    df_raw = df_raw[df_raw['blastp-identity'] == 100]
    # remove isoforms
    df_raw = df_raw[~df_raw['Uniprot-Acc'].str.contains('-', na=False)]
    
    print ('getting unique sub')
    df_unique_sub = df_raw[['identifier','Uniprot-Acc']]
    df_unique_sub = df_unique_sub.drop_duplicates()
    print ('getting sub_acc')
    start = time.time()
    # get the uniprotID of the protein submitted for PhosphoPICK prediction
    df_unique_sub['acc_from_identifier'] = df_unique_sub.apply(lambda row :  row['identifier'].split("|")[1], axis = 1) 
    # retrieve the current uniprotID for the predicted substrate
    df_unique_sub['substrate_acc'] = df_unique_sub.apply(lambda row :  getUniprotID.getUniprotID(str(row['Uniprot-Acc']),'other'), axis = 1) 
    # only keep perdictions that the perdicted substrate is the input protein (remove any data with outdated uniprotID, and perdicted portein has seq with 100% blastp-identity but is not the input portein)
    df_unique_sub = df_unique_sub[df_unique_sub['acc_from_identifier'] == df_unique_sub['substrate_acc']]
    print ('merge')
    df_raw = df_raw.merge(df_unique_sub, left_on=['identifier','Uniprot-Acc'], right_on=['identifier', 'Uniprot-Acc'], how = 'right')
    df_raw = df_raw.drop(columns = ['identifier', 'acc_from_identifier'])
    end = time.time()
    print ('done', end-start)

    print ('getting unique kin')    
    df_unique_kin = df_raw[['kinase']]
    df_unique_kin = df_unique_kin.drop_duplicates()
    print ('getting kin_acc')
    start = time.time()
    # retrieve the current uniprotID for the predicted kinases
    df_unique_kin['kinase_acc'] = df_unique_kin.apply(lambda row :  getUniprotID.getUniprotID(str(row['kinase']), 'gene name'), axis = 1) 
    print ('merge')
    df_raw = df_raw.merge(df_unique_kin, left_on=['kinase'], right_on=['kinase'], how = 'left')
    end = time.time()
    print ('done', end-start) 

    return df_raw

def map_site (filename,  human_proteome):
    """
    Given PhosphoPICK dataframe file, positionally map the site to the given human proteome sequence dataframe
    
    Parameters
    ----------
    filename : str
        phosphoPICK file
    human_proteome : str
        given human proteome sequence dataframe

        
    Returns
    ----------
    phosphoPICK: pandas df
                'Uniprot-Acc'       the uniprotID of the predicted substrate by PhosphoPICK
                'substrate_id'      unique IDs for the substrate phosphorylation site (substrate_acc + position)
                'substrate_acc'     substrate uniprotID
                'site'              aa + position in protein sequence
                'position'          position in protein sequence
                'pep'               +/- 7 AA
                'kinase'            kinase name used in PhosphoPICK
                'kinase_acc'        kinase uniprotID
                'combined-p-value'  PhosphoPICK score 
        
    """
    start = time.time()
    print ('Reading input file...')
    df_temp = pd.read_csv(filename)
    print ('Get unique substrate sites...')

    # get a list of unique substrate/site in PhosphoPICK
    df_unique_sub = df_temp[['substrate_acc', 'position']]
    df_unique_sub = df_unique_sub.drop_duplicates()
    
    for index, row in df_unique_sub.iterrows():
        id = df_unique_sub.at[index, 'substrate_acc']
        try: 
            # get the sequence of the given uniprotID
            seq = checkSite.getSequence(id, human_proteome)
            pos = int(df_unique_sub.at[index, 'position'])
            # get the substrate site (aa + pos)
            df_unique_sub.at[index, 'site'] = seq[pos - 1] + str(pos)
            # get the +/- 7 AA around the site
            pep = checkSite.getPep(pos, 8, seq)
            df_unique_sub.at[index, 'pep'] = pep
            df_unique_sub.at[index, 'substrate_id'] = id + '_' + str(pos)
        # except any outdated records such as: 
        # . deleted uniprotID
        # . updated UniprotID
        # . sequence of the given uniprotID changed causing out of range error
        except:
            df_unique_sub.at[index, 'substrate_id'] = 'outdated'

    # add the 'substrate_id' to phosphoPICK df 
    df_temp = df_temp.merge(df_unique_sub, left_on=['substrate_acc', 'position'], right_on=['substrate_acc', 'position'], how = 'left')
          
    return df_temp

def pick_convert_directory(load_directory, reference_filename, save_directory, convert_type):
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
        location to save the phosophoPICK result dataframe
    convert_type : str
        acc (convert substrate/kinase accession)
        OR, site (mapping site to reference sequence)
    """
    if not os.path.exists(save_directory):
                os.mkdir(save_directory) 
    # for converting substrate and kinase accession
    if convert_type == 'acc':
        all_files = glob.glob(load_directory + "*.txt")
        for filename in all_files:
            start = time.time()
            name = re.search(r'.+\/(.+)$', filename).group(1)
            print ("Formatting ", name, "...")
            df = convert_acc(filename)
   
            df.to_csv(save_directory + name, index=False)
            
            end = time.time()
            print (f"Done. Time\t{(end-start):.3f}")
    # for mapping the site to the reference human proteome seq
    elif convert_type == 'site':
        HP_df = pd.read_csv(reference_filename, sep = '\t')

        all_files = glob.glob(load_directory + "*.txt")
        for filename in all_files:
            start = time.time()
            name = re.search(r'.+\/(.+)\.txt$', filename).group(1)
            print ("Formatting ", name, "...")
            df = map_site(filename, HP_df)
   
            df.to_csv(save_directory + name + '_mappedSite.csv', index=False)
            
            end = time.time()
            print (f"Done. Time\t{(end-start):.3f}")
