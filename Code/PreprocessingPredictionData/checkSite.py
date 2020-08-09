import pandas as pd
from Bio import SeqIO
import re
from re import finditer
from urllib.request import urlopen
import urllib.error
from io import StringIO
import requests
import time

def checkSite(id, site, pep, pos_in_pep, HP_df):
    """
    positionally map the given substrate site to the reference humnan proteome seq
    
    Parameters
    ----------
    id : str
        substrate_acc (uniprotID)
    site : str or int
        str : aa + postion in protein seq
        int: postion in protein seq
    pep : str
        given pep seq around the site
    pos_in_pep : int
        position in the peptide
    HP_df : df
        the df of the reference human proteome seq
        
    Returns
    -------
    new_site: str
    site_confirm: bool
    """

    site_confirm = False
    new_site = ""

    if (site != 'nan'):  # skip NaN cells
        site_prep = re.search(r'(\d+)', site)
        # site contains aa and/or pos
        if site_prep:
            aa = pep[pos_in_pep-1]
            pos = int(site_prep.group(1))
        # site contains only aa
        else:
            site_prep = re.search(r'(\w)', site)
            aa = site_prep.group(1)
            new_site = "error" + aa  # no position
            site_confirm = False
            return new_site, site_confirm

        # remove "_" at the start and end of pep
        pep = re.sub(r"_+$", "",pep)
        if re.match(r'^_+',pep):
            x = re.match(r'^_+',pep)
            count= x.end()
            pep = re.sub(r"^_+", "",pep)
            pos_in_pep = pos_in_pep - count
        # get the seq of the given substrate_acc
        ref_seq = getSequence(id, HP_df)
        
        # if the substrate_acc exist
        if ref_seq != "":
            try:
                # try to get pep around the site in the reference seq
                ref_pep = getPep(pos, pos_in_pep, ref_seq)
            except:
                new_site = "(posOutOfRange)" + site
                site_confirm = False
                print (id, ": ", " corrected: " , new_site, " site: ", site, " seq len: ", len(ref_seq), 'pep: ', pep)
                return new_site, site_confirm

            if pep == ref_pep:
                site_confirm = True
                site = aa + str(pos)
                return site, site_confirm
            # if the pep around the site in the reference seq does not match the given pep 
            else:
                # check if the seq shifts
                site_shift = seqShift (pep,ref_seq,pos)
                # the given pep is not found
                if site_shift == -1:
                    new_site = "(seqNotFound)" + site
                    site_confirm = False
                    print (id, ": ", " corrected: " , new_site, " site: ", site, " seq len: ", len(ref_seq), 'pep: ', pep)
                else:
                    new_pos = pos_in_pep + site_shift
                    new_site = aa + str(new_pos)
                    site_confirm = True
                    #print (id, ": ", " site: ", site, " corrected: " , new_site, " shifted: ", new_pos - pos, " seq len: ", len(ref_seq), 'pep: ', pep)
        else:
            new_site = "(idNotFound)" + site
            site_confirm = False
            print (id, ": ", " corrected: " , new_site, " site: ", site, 'pep: ', pep)
    else:
        new_site = ""  # empty cells
        site_confirm = True
    
    return new_site, site_confirm

def seqShift (pattern, seq, pos):
    """
    check if the position of the given seq shifted in the reference human proteome seq
    returns the number of shifted postion that has the smallest differece from the given postion  
    
    Parameters
    ----------
    pattern : str
        given pep seq
    seq : str
        position in the peptide
    pos : int
        given positon in the sequence
        
    Returns
    ---------
    site_shift: int
        number of postion shifted
        
    """
    site_shift = -1
    old_site_shift = 0
    new_site_shift = 0
    pos_compare = 99999999999
    pattern = pattern.upper()
    for match in finditer(pattern, seq):
        new_site_shift = match.start()
        if abs(new_site_shift - pos) < pos_compare:
            pos_compare = abs(new_site_shift - pos)
            old_site_shift = new_site_shift
        site_shift = old_site_shift

    return site_shift

def getSequence(id, HP_df):
    """
    get the protein seq from the reference human prteome by given uniprot accession  
    
    Parameters
    ----------
    id : str
        protein uniprot accession
    HP_df : dataframe
        reference human prteome
        
    Returns
    -------
    seq: str
        protein sequence of given uniprot accession
    """
    try:
        seq = HP_df.loc[HP_df['UniprotID'] == id, 'sequence'].iloc[0]
    except:
        seq = ''  
        return seq  

    return seq
    
def getPep (site, pos_in_pep, ref_seq):
    """
    get the pep seq arount the given site from the reference sequence  
    
    Parameters
    ----------
    site : int
        pos in protein seq
    pos_in_pep : int
        position in the peptide (defines the number of aa will get around the site)
    ref_seq : str
        reference sequence
        
    Returns
    -------
    pep_seq: str
        pep seq arount the given site
    """
    # get the start of the pep seq in the provided seq (ref_seq)
    start = site - pos_in_pep 
    # get the end of the pep seq in the provided seq (ref_seq)
    end = site + (pos_in_pep - 1)

    # if the seq is shorter than 15 aa
    if start < 0 and end > len(ref_seq):
        replace = abs(start) # number of '_' append to the beginning of the pep seq
        pep_seq = "_" * replace + ref_seq[0:end]
        replace = end - len(ref_seq) # number of '_' append to the end of the pep seq
        pep_seq = pep_seq + "_" * replace
    # if the site is at the very beginning of the protein seq (within the first 7 aa)
    elif start < 0:
        replace = abs(start) # number of '_' append to the beginning of the pep seq
        pep_seq = "_" * replace + ref_seq[0:end]
    # if the site is at the very end of the protein seq (within the last 7 aa)
    elif end > len(ref_seq):
        replace = end - len(ref_seq) # number of '_' append to the end of the pep seq
        pep_seq = ref_seq[start:len(ref_seq)]+"_"*replace
    # if the site is in the middle of the protein seq
    else:
        pep_seq = ref_seq[start:end]

    return pep_seq