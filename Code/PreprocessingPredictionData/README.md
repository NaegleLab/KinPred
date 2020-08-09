- **CreateHumanProteomeDF.ipynb**:
  1. Download the current human proteomes as fasta file from UniProt.org. Convert it to a tab delimited csv file. The final output file will be used as the reference for substrate and kinase ID mapping. The dataframe contains the following columns: UniprotID, Gene Name, Entry Name, sequence. 
  2. Split the downloaded fasta file into smaller files with 1000 sequences in each (the last one has 367 sequences). The split fasta files are submitted to predictiors to retrieve raw predictions.
- **humanProteomesReference.py**
  - *downloadHumanProteomes*: Downloads the Human Proteomes (canonical) from Unipro.org, and saves as fasta formate at the given dir/name.
  - *downloadFasta*: Downloads the protein seq by uniprotID from Unipro.org, and saves as fasta formate at the given dir/name.
  - *batch_iterator, splitFasta*: Splits and saves the input file as multiple files with given size(number of sequences)
  - *fastaToCSV*: Convert the input fasta file into a dataframe and returned
- **FormattingGPS.ipynb**: Pre-process the GPS raw prediction into a standard format across predictors.  The standard formatted file includes all predictions on human kinase-substrate interaction in which the predicted sites are mapped to the fixed version of human proteome.  The standard formatted file contains information of unique IDs for the predicted phosphorylation site (substrate protein accession + position in protein seq), gene name for the substrates, Uniprot accessions for the substrates, site (aa + position in protein seq), peptide sequences around the sites, scores, common kinase names use across all predictors
- **gps_convert.py**: pulls GPS raw prediction files (in fasta like format), convert them into list-like dataframe.  The module also includes functions to map the substrate/kinase uniprot accession and site position to the fixed version of human proteome.
- **FormattingPhosphoPICK.ipynb**: Pre-process the PhosphoPICK raw prediction into a standard format across predictors.
- **phosphoPick_convert.py**: pulls PhosphoPICK raw prediction files, map the substrate/kinase uniprot accession and site position to the fixed version of human proteome.
- **FormattingNetworKIN.ipynb**: Pre-process the NetworKIN raw prediction into a standard format across predictors.
- **networKin_convert.py**: pulls NetworKIN raw prediction files, map the substrate/kinase uniprot accession and site position to the fixed version of human proteome.
- **getUniprotID.py**: get the Uniprot accession by query the given 'id' in Uniprot.org. 'id' can be other types of accession, or protein name
- **checkSite.py**: positionally map the given substrate site to the reference humnan proteome seq
  - *seqShift*: check if the position of the given seq shifted in the reference human proteome seq, returns the number of shifted postion that has the smallest differece from the given postion 
  - *getSequence*: get protein seq by from ref human proteome by uniprot accession
  - *getPep*: get the pep seq arount the given site from the reference sequence

