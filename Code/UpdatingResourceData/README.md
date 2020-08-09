- **UpdateHumanProteomeVersion.ipynb**:
  1. Update the fixed version Human Proteome to the current UniProt version
  2. Update KinPred for subatrates with sequence changes in the updated Human Proteome
- **UpdataNewPhosphosites.ipynb**: update KinPred with newly discovered phosphosites
- **humanProteomesReference.py**
  - *downloadHumanProteomes*: Downloads the Human Proteomes (canonical) from Unipro.org, and saves as fasta formate at the given dir/name.
  - *downloadFasta*: Downloads the protein seq by uniprotID from Unipro.org, and saves as fasta formate at the given dir/name.
  - *batch_iterator, splitFasta*: Splits and saves the input file as multiple files with given size(number of sequences)
  - *fastaToCSV*: Convert the input fasta file into a dataframe and returned
- **gps_convert.py**: pulls GPS raw prediction files (in fasta like format), convert them into list-like dataframe.  The module also includes functions to map the substrate/kinase uniprot accession and site position to the fixed version of human proteome.
- **phosphoPick_convert.py**: pulls PhosphoPICK raw prediction files, map the substrate/kinase uniprot accession and site position to the fixed version of human proteome.
- **networKin_convert.py**: pulls NetworKIN raw prediction files, map the substrate/kinase uniprot accession and site position to the fixed version of human proteome.
- **getUniprotID.py**: get the Uniprot accession by query the given 'id' in Uniprot.org. 'id' can be other types of accession, or protein name
- **checkSite.py**: positionally map the given substrate site to the reference humnan proteome seq
  - *seqShift*: check if the position of the given seq shifted in the reference human proteome seq, returns the number of shifted postion that has the smallest differece from the given postion 
  - *getSequence*: get protein seq by from ref human proteome by uniprot accession
  - *getPep*: get the pep seq arount the given site from the reference sequence

