import pandas as pd
import requests
import re
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import date


def downloadHumanProteomes (saveFileAs):
    """
    Downloads the Human Proteomes (canonical) from Unipro.org, and saves as fasta formate at the given dir/name.
    
    Parameters
    ----------
    saveFileAs : str
        dir/name of the downloaded fasta file

    """
    # set the url for the query
    url = 'http://www.uniprot.org/uniprot/'
    # set the parameters for the query
    payload = {
        "query": 'reviewed:yes AND organism:"Homo sapiens (Human) [9606]"',
        "format": "fasta"
    }
    # get the result of the query 
    result = requests.get(url, params=payload)
    # write the result in a fasta file
    if result.ok:
        file = open(saveFileAs, 'w')
        file.write(result.text)
        file.close()
    else:
        print('Something went wrong: ', result.status_code)

def downloadFasta (id_df, fileName):
    """
    Downloads the protein seq by uniprotID from Unipro.org, and saves as fasta formate at the given dir/name.
    
    Parameters
    ----------
    id_df: str
        uniprotIDs
     out_dir : str
        saving dir of the downloaded fasta file

    """
    file = open(fileName, 'w')
    # set the url for the query
    url = 'http://www.uniprot.org/uniprot/'
    
    for index, row in id_df.iterrows():
        id = id_df.at[index, 'UniprotID']
        # set the parameters for the query
        payload = {
            "query": id,
            "format": "fasta"
        }
        # get the result of the query 
        result = requests.get(url, params=payload)
        # write the result in a fasta file
        if result.ok:
            file.write(result.text)
        else:
            print('Something went wrong: ', result.status_code)
    file.close()


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    Parameters
    ----------
    iterator : sequences in the file that needed to split
    batch_size : int
        
    Returns
    ----------
    lists of batch_size entries

    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def splitFasta(fileName, outputDir, size):
    """Splits and saves the input file as multiple files with given size(number of sequences)

    Parameters
    ----------
    fileName : str
        the file that needs split
    outputDir : str
        the dir where the splited files saved at
    size : int
        number of senquences in each file (the last file may be shorter)

    """

    record_iter = SeqIO.parse(open(fileName, 'r'),"fasta")

    for i, batch in enumerate(batch_iterator(record_iter, size)):
        file = outputDir + 'group_%i.fasta' % (i + 1)
        with open(file, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, file))

def fastaToCSV (input, output):
    """Convert the input fasta file into a dataframe and returned

    Parameters
    ----------
    input : str
        the sequence file in fasta formate
    output : str
        the sequence file in csv formate
    
    Returns
    ----------
    human proteomes pandas df
        'UniprotID'
        'Gene Name' 
        'Entry Name'
        'sequence'

    """
    # convert the HumanProteome fasta file into a tab delimited csv
    with open (output,"w") as seq_df:
        # create header/column names
        seq_df.write('UniprotID' + '\t' + 'Gene Name' + '\t' + 'Entry Name' + '\t' + 'sequence' + '\n')
        with open(input) as fasta_file:  
            for title, sequence in SimpleFastaParser(fasta_file):
                identifier = re.match(r'sp\|(\S+)\|(\S+_HUMAN).+GN\=(\S+).+', title)
                # if there is a "Gene name"
                if identifier:  
                    uniprotID = identifier.group(1)
                    geneName = identifier.group(3)
                    entryName= identifier.group(2)
                # if there is no "Gene Name"
                else:
                    identifier = re.match(r'sp\|(\S+)\|(\S+_HUMAN).+', title)
                    uniprotID = identifier.group(1)
                    geneName = "N/A"
                    entryName= identifier.group(2)

                seq_df.write(uniprotID + '\t' + geneName + '\t' + entryName + '\t' + sequence + "\n")
    
    seq_df = pd.read_csv(output, sep = "\t")
    return seq_df