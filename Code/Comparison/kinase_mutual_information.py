import pandas as pd
import numpy as np

def between_predictors(network1, network2, kinase_column, substrate_column, output):
    """
    Finds mutual information shared between kinases based on the substrate phosphorylated
    Mutual Information is defined as the intersection substrates between two kinases
    A substrate is defined as the substrate accession and site, i.e. P54760_596.
    Normalization is performed by comparing intersection of kinases vs union of the two kinases
    This the the Jaccard Index. Jaccard Distance can be calcualted by taking 1 - JI
    Parameters
    ----------
    network1 : pandas dataframe
        The network to analyze for mutual kinase information
    network2 : pandas dataframe
        The network to analyze for mutual kinase information
    kinase_column : str
        Column in network that contiains kinase information
    substrate_column : str
        Column in network that contains substrate information
    
    Returns
    --------
    heatmap : pandas dataframe
        Normalized mutual information into Jaccard Index. 
        size of intersection of two kinase networks / size of union of two kinase networks.
    """
    kinases1 = list(network1[kinase_column].unique())
    
    kinases2 = list(network2[kinase_column].unique())
    
    all_kin = list(set(kinases1) | set(kinases2)) 
    num_all_kin = len(all_kin)
    
    overlapped_kin = list(set(kinases1) & set(kinases2))
    
    normalized = np.zeros((num_all_kin, num_all_kin))
    
    for i in range(num_all_kin):
        for j in range(num_all_kin): 
            if all_kin[i] in overlapped_kin and all_kin[j] in overlapped_kin:
                substrates1 = set(network1[network1[kinase_column] == all_kin[i]][substrate_column])
                substrates2 = set(network2[network2[kinase_column] == all_kin[j]][substrate_column])
                normalized[i][j] = len(substrates1.intersection(substrates2)) / len(substrates1.union(substrates2))
            else:
                normalized[i][j] = -1
    
    heatmap = pd.DataFrame(normalized, index = all_kin, columns = all_kin)
    heatmap.to_csv(output + '.csv')
                                                    
    return heatmap

def within__predictors(network, kinase_column, substrate_column, output):
    """
    Finds mutual information shared between kinases based on the substrate phosphorylated
    Mutual Information is defined as the intersection substrates between two kinases
    A substrate is defined as the substrate accession and site, i.e. P54760_596.
    Normalization is performed by comparing intersection of kinases vs union of the two kinases
    This the the Jaccard Index. Jaccard Distance can be calcualted by taking 1 - JI
    Parameters
    ----------
    network : pandas dataframe
        The network to analyze for mutual kinase information
    kinase_column : str
        Column in network that contiains kinase information
    substrate_column : str
        Column in network that contains substrate information
    
    Returns
    --------
    heatmap : pandas dataframe
        Number of substrates that overlap between kinases
    normalized : pandas dataframe
        Normalized mutual information into Jaccard Index. 
        size of intersection of two kinase networks / size of union of two kinase networks.
    heatlist : pandas dataframe
        intersction of kinase networks
    """
    kinases = list(network[kinase_column].unique())
    num_kinases = len(kinases)
    
    heatlist = [[set() for i in range(num_kinases)] for j in range(num_kinases)]                                                           
    heatmap = np.zeros((num_kinases, num_kinases))
    normalized = np.zeros((num_kinases, num_kinases))
    
    
    for i in range(num_kinases):                                                                                    # Get kinase-substrate network of each kinase
        substrates = set(network[network[kinase_column] == kinases[i]][substrate_column])
        heatlist[i][i] = substrates
    
    for i in range(num_kinases):                                                                                    # Iterate through each row in heatmaps
        for j in range(num_kinases):                                                                                # Iterate through each column in heatmaps
            heatlist[i][j] = heatlist[i][i].intersection(heatlist[j][j])                                            # Find intersection of kinase networks
            heatmap[i][j] = len(heatlist[i][j])                                                                     # Get size of intersection
            normalized[i][j] = len(heatlist[i][j]) / len(heatlist[i][i].union(heatlist[j][j]))                      # Normalize via size of the union of two kinase networks
    
    heatmap = pd.DataFrame(heatmap, index = kinases, columns = kinases)
    normalized = pd.DataFrame(normalized, index = kinases, columns = kinases)
    heatlist = pd.DataFrame(heatlist, index = kinases, columns = kinases)

    normalized.to_csv(output + '.csv')

    return heatmap, normalized, heatlist