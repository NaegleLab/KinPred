All the raw data and output files of the study should store under "Data". The "Data" derictory should have the following structure:

- **Data**
  - **Raw**:  
    - **GPS**: 
      - 21 GPS5.0 raw prediction files (please see `Get GPS Result Files` section under `Instruction` in [FormattingGPS.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingGPS.ipynb) ), 
      - GPS kinase predictor cutoff score files (available at: https://doi.org/10.6084/m9.figshare.12749342.v2, Raw(unfiltered) Data Files for KinPred)
    - **PhosphoPICK**: 
      - 21 PhosphoPICK raw prediction files (please see `Get PhosphPICK Result Files` section under `Instruction` in [FormattingPhosphoPICK.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingPhosphoPICK.ipynb) )
    - **NetworKIN**: 
      - 21 NetworKIN raw predictions files (please see `Get NetworKIN Result Files` section under `Instruction` in [FormattingNetworKIN.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingNetworKIN.ipynb) )
    - **HumanKinase**: 
      - globalKinaseMap.txt: manually created list of human kinases, this is the base of the final globalKinaseMap.csv
    - **HumanProteome**: 
    downloaded fasta files (available at: https://doi.org/10.6084/m9.figshare.12749342.v2, Raw(unfiltered) Data Files for KinPred)
    - **ProteomeScout2020-02-07**: 
    dowloaded ProteomScout reference files (available at: https://doi.org/10.6084/m9.figshare.12749342.v2, Raw(unfiltered) Data Files for KinPred)
  - **Temp**: all temp files generated programmatically. files for each predictor saved under separete dir 
    - **GPS**: 
      - **mappedAcc**: after map the substrate and kinase protein accessions 
      - **mappedSite**: after map the prediction site to the fix version of human proteome 
    - **PhosphoPICK**
      - **mappedAcc**
      - **mappedSite** 
    - **NetworKIN**: 21 NetworKIN raw predictions files
      - **mappedAcc**
      - **mappedSite** 
   - **Map**:
     - globalKinaseMap.csv: global kinase ontology that map between predictor-specific kinase names and a common kinase name. It contains kinase uniprot accession, common kinase name, preferred kinase names, description of kinase, kinase type,  reference kinase names in each predictor. (available at: https://doi.org/10.6084/m9.figshare.12749333.v1, Kinase Ontology). 
     - humanProteome dataframe converted from fasta file. It contains the uniprot accession, gene name, entry name, and sequence
   - **Formatted**: all human substrate-kinase predictions for the three predictors in a standered formate which contains a unique substrate id, substrate protein uniprot accession, substrate gene name, predicted site, peptide around the site, predicted kinase, and score. (to generate these files, please follow instructions in [FormattingGPS.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingGPS.ipynb), [FormattingPhosphoPICK.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingPhosphoPICK.ipynb), [FormattingNetworKIN.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingNetworKIN.ipynb))
      - **GPS**
      - **PhosphoPICK**
      - **NetworKIN**
  - **Final**: Final prediction matrix with prediction of know phosphosites (available at: https://doi.org/10.6084/m9.figshare.12749324.v1, Final Data - Matrix of kinase substrate edge weights by predictor).
    - **GPS**
    - **PhosphoPICK**
    - **NetworKIN**
  - **comparison**: output files for all the anaysis 
    - **Thresh**: filtered prediction by different thresholds, csv files
      - **GPS**
        - **all**: with all know phsophosites
        - **2exp**: with phsophosites that have at least 2 experimentail evidence
        - **3exp**: with phsophosites that have at least 3 experimentail evidence
      - **PhosphoPICK**
        - **all**
        - **2exp**
        - **3exp**
      - **NetworKIN**
        - **all**
        - **2exp**
        - **3exp**
    - **kinaseNetworkSize**: figures, svg files
    - **substrateNetworkSize**: figures, svg files
      - **box**
      - **contour**
      - **violin**
    - **kinaseNetworkSimilarity_withinPredictor**
      - **dataframe**: Jaccard's Index matries files, csv files
      - **heatmaps**: figure, eps files
      - **networks**: figure, eps files
    - **kinaseNetworkSimilarity_acrossPredictor**
      - **dataframe**: Jaccard's Index matries files, csv files
      - **heatmaps**: figure, eps files
      - **rank**: figure, svg files
        - **dataframe**: rank summaries, csv files
    - **randomizeNetwork**
      - **dataframe**: dataframe of the randomized networks, csv files
      - **kinaseNetworkSimilarity_acrossPredictor**: similarity comparsions between randomized networks, figure, svg files
        - **rank**: rank summaries, csv files; figure, svg files 
    - **studyBias**: figure, svg files
