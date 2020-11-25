import sys
sys.path.append('../../../ProteomeScoutAPI/')
from proteomeScoutAPI import ProteomeScoutAPI
sys.path.append('../../../KinaseActivity2019/')
from kinase_activity.src import experiment
sys.path.append('../PreprocessingPredictionData/')
import checkSite, humanProteomesReference
import createSubKinMatrix
import pandas as pd
import re
import os
from datetime import date
from Bio.SeqIO.FastaIO import SimpleFastaParser
from urllib.request import urlopen
from zipfile import ZipFile
from io import BytesIO

def getPScoutData():
    """
    download and unzip current ProteomeScout data

    Output:
    ------------------------------------------
    download and unzipped ProteomeScout data are saved in 
    `/Data/ProteomeScout_Update` directory:
        data.tsv
        citations.tsv
    
    """
    # ProteomeScout phosphorylation data zip file downloading address
    data_add = 'https://proteomescout.wustl.edu/compendia/proteomescout_phosphorylation.zip'
    # phosphorylation data file name in the zip file
    data_file = 'data'
    citation_file = 'citations'
    # output dir
    dir = '../../Data/Raw/ProteomeScout_'+date.today().strftime('%Y-%m-%d')

    resp = urlopen(data_add)
    zipfile = ZipFile(BytesIO(resp.read()))

    if not os.path.exists(dir):
        os.mkdir(dir) 

    file = open(os.path.join(dir, data_file+".tsv"), 'w')
    for line in zipfile.open(data_file+".tsv").readlines():
        file.write(line.decode('utf-8'))
    file.close()

    file = open(os.path.join(dir, citation_file+".tsv"), 'w')
    for line in zipfile.open(citation_file+".tsv").readlines():
        file.write(line.decode('utf-8'))
    file.close()

def getHumanPTMs(pscout_data, ref_proteome):
    """
    get the current human phosphoproteome data
    map the site to the reference human proteome

    Parameters
    ----------
    pscout_data : str
        path to the ProteomeScout data
    ref_proteome : str
        path to the reference human proteome

    Return 
    -----------
    df: dataframe of the ProteomeScout data that 
        mapped to the reference human proteome
    
    """
    dir = '../../Data/Map/HumanProteome/'
    HP_df = '../../Data/Map/HumanProteome/' + os.path.splitext(os.path.basename(ref_proteome))[0]+'.csv'
    PS_update = os.path.dirname(pscout_data)+'Human_PhosphOme_all.csv'

    
    # get the list of substrate protein uniprotID from the reference human proteome
    if not os.path.exists(dir):
        os.mkdir(dir) 
    human_proteome_df = humanProteomesReference.fastaToCSV(ref_proteome, HP_df )

    accessions = human_proteome_df['UniprotID']

    # for every human protein, get the phosphosites and append to a dataframe, can require they have at least some number
    # of pieces of evidence
    PTM_API = ProteomeScoutAPI(pscout_data)
    numEvidences = 1 
    array = []
    for acc in accessions: 
        PTMs = PTM_API.get_PTMs_withEvidenceThreshold(acc, numEvidences)
        #mods comes back as list of lists, each with information about a site, e.g. ('37', 'Y', 'Phosphotyrosine')
        seq = PTM_API.get_sequence(acc)
        if isinstance(PTMs, list):
            for PTM in PTMs:
                if PTM[2] in ['Phosphoserine', 'Phosphothreonine', 'Phosphotyrosine']:
                    site = int(PTM[0])
                    aaSite = "%s%d"%(PTM[1], site)
                    pep = experiment.get_aligned_peptide(aaSite, seq, 7) #e.g. set to 7 for a 15-mer
                    idx = [i for i, a in enumerate(pep) if a.islower()][0]
                    if idx < 7:
                        replace = 7-idx # number of '_' append to the beginning of the pep seq
                        pep = "_" * replace + pep[0:len(pep)]
                    elif idx == 7 and len(pep) < 15 :
                        replace = 15 - len(pep) # number of '_' append to the end of the pep seq
                        pep = pep[0:len(pep)]+"_"*replace
                    pep = pep.upper()
                    pos_in_pep = 8 #e.g. set to 8 for a 15-mer
                    new_site, site_confirm = checkSite.checkSite(acc, aaSite, pep, pos_in_pep, human_proteome_df) #e.g. set to 7 for a 15-mer
                    # if the site mapped to the referece human proteome seq at the exact position, the new_site = site, site_confirm = True
                    # if the peptide mapped to the referece human proteome seq with a shift in position, new_site = site + shift, site_confirm = True
                    if site_confirm == True:
                        aaSite = new_site
                        # substrate_id = substrate_acc + position in protein seq
                        substrate_id = acc + '_' + str(int(re.search(r'(\d+)', aaSite).group(1))) 
                    # if the peptide can not mapped to the referece human proteome seq, site_confirm = False
                    else:
                        aaSite = new_site
                        substrate_id = 'outdated'
                    
                    array.append([substrate_id, acc, aaSite, pep])
                
    df = pd.DataFrame(array, columns=['substrate_id', 'acc', 'site', 'pep'])

    df.to_csv(PS_update)
    return df

def XRefProteomeScout(pscout_data, ref_proteome, old_version):
    """
    Cross reference with ProteomeScout
    the final formatted file only contains data with confirmed phosphorylation site
    Create final substrate/kinase matrices for PhosphoPICK, GPS, and NetworKIN 

    Parameters
    ----------
    pscout_data : str
        path to the ProteomeScout data
    ref_proteome : str
        path to the reference human proteome
    old_version: str (YYYY-MM-DD)
        input (unfiltered) file version

    Output
    -----------
    final standard formated files saved in 
    `../../Data/Final/` directory
    
    """
    
    # preprocessed prediction data 
    in_dir = '../../Data/Formatted/'
    PP_file = in_dir + 'PhosphoPICK/PhosphoPICK_formatted_' + old_version + '.csv'  # PhosphoPICK files      
    GPS_file = in_dir + 'GPS/GPS_formatted_' + old_version + '.csv'                 # GPS files
    NW_file = in_dir + 'NetworKIN/NetworKIN_formatted_' + old_version + '.csv'      # NetworKIN files

    # output dir
    out_dir = '../../Data/Final/'
    pp = out_dir+'PhosphoPICK/PhosphoPICK_' + date.today().strftime('%Y-%m-%d')
    gps = out_dir+'GPS/GPS_' + date.today().strftime('%Y-%m-%d')
    nw = out_dir+'NetworKIN/NetworKIN_' + date.today().strftime('%Y-%m-%d')
    

    pscout_df = getHumanPTMs(pscout_data, ref_proteome)
    file_list = [PP_file, GPS_file, NW_file]
    out_list = [pp, gps, nw]

    i = 0
    for file in file_list:
        df = pd.read_csv(file)
        df = df[df['substrate_id'].isin(pscout_df['substrate_id'])]
        if not os.path.exists(os.path.dirname(out_list[i])):
            os.mkdir(os.path.dirname(out_list[i])) 
        df.to_csv(out_list[i] + '.csv')
        df_matrix = createSubKinMatrix.createMatrix(out_list[i] + '.csv', out_list[i] + '_matrix.csv')
        i += 1 