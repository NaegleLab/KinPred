import sys 
sys.path.append('../../../ProteomeScoutAPI/')
import proteomeScoutAPI
import re
import pandas as pd

def pull_proteome_scout_phosphorylation(filename, phospho_sites = ['Y', 'S','T']):
    """
    Finds all phosphorylation sites in the given Proteome Scout flat file
    Parameters
    -----------
    filename : str
        ProteomeScout flatfile
    proteins : list
        proteins search 
    """
    ps = proteomeScoutAPI(filename)

    df_full = pd.DataFrame(columns = ['protein_id', 'amino_acid_id'])

    for ID in PTM_API.uniqueKeys():
        for site in phospho_sites:
            df_full = df_full.append( phosphorylation_data(ps, ID, site), ignore_index=True)
    
    return df_full


def phosphorylation_data(ps_api, protein_id, filter, site_amino = 'Y'):
        """
        Builds dataframe that finds all amino acid and shows whether 
        site has evidence that it has undergone Phosphorylation

        Parameters
        ----------
        protein_id : string
           Protein id 
        site_amino : char
            eg. Tyrosine : Y
        num_sequence : int
            number to look forwards and backwards to generate sequence

        Returns
        -------
        pandas dataframe
            Protein 
                provided in method
            Amino_ID
                Amino acid + location in sequence (Y23)
            Sequence
                Sequence of Amino Acids surrounding core Amino Acid
        """
        num_sequence = filter[2]
        df_seq = build_amino_sequence(ps_api, protein_id, site_amino, num_sequence)

        if df_seq.empty:
            return df_seq

        #get phosphosites by evidence threshold
        if filter[0] == 'number':
            phos = ps_api.get_PTMs_withEvidenceThreshold(protein_id, filter[1])
        #get phosphosites by evidence source
        elif filter[0] == 'type':
            modList = ps_api.get_PTMs_withEvidence(protein_id)
            experimentList = filter[1]
            phos = []
            if isinstance(modList, list):
                for mod in modList:
                    PTM = mod['mod']
                    evidence = mod['evidence']
                    for e in evidence:
                        if e in experimentList:
                            phos.append(PTM)

        df_phos = pd.DataFrame(phos, columns= ['Location', 'AminoAcid', 'Modification'])
        df_phos['site'] = df_phos['AminoAcid'] + df_phos['Location']
        df_phos['phosphorylation'] = 1
        # # print df

        df_comb = pd.merge(df_seq, df_phos, how = 'left', on= 'site')

        df_comb.fillna(0, inplace = True)
        df_comb['phosphorylation'] = df_comb['phosphorylation'].astype(int)

        df_comb = df_comb[['substrate_id', 'site', 'peptide', 'phosphorylation']]
        
        return df_comb

def build_amino_sequence(ps_api, protein_id, amino_acid, num_sequence):
        """
        Generates dataframe of amino sequences surrounding desired amino acid. 
        Data is pulled from ProteomeScout flatfile using ProteomeScoutAPI
        
        Parameters
        ----------
        ps_api
            proteomeScoutAPI object 
        protein_id : string
           Protein id 
        amino_acid_id : char
            eg. Tyrosine : Y
        num_sequence : int
            number to look forwards and backwards to generate sequence
        
        Returns
        -------
        pandas dataframe 
            substrate_id 
                provided in method
            site
                Amino acid + location in sequence (Y23)
            sequence
                Sequence of Amino Acids surrounding core Amino Acid
        """
        sequence = ps_api.get_sequence(protein_id)
        columns = ['substrate_id', 'site', 'peptide']
        df = pd.DataFrame(columns=columns)
        seq_len = len(sequence)
        for i in range(len(sequence)):
            if sequence[i] == amino_acid:
                
                min = i - num_sequence
                max = i + num_sequence + 1

                min = 0 if min < 0 else min
                max = seq_len if  max > seq_len else max
                df = df.append({columns[0]: protein_id, 
                                columns[1]: amino_acid + str(i + 1), 
                                columns[2]: sequence[min:max]}, ignore_index=True)
        return df