import requests
import re
from requests.exceptions import ConnectionError
import time
from io import StringIO 
import pandas as pd

def getUniprotID(id, id_type):
    """
    get the UniprotID by query the given 'id' in Uniprot.org
    
    Parameters
    ----------
    id : str
        the search term 
    id_type : str
        the search type:
            gene name (gene/protein name)
            other (other type of accession)
    """
    reviewed = '+reviewed:yes'
    # set the search terms
    
    if (id != 'nan'):  # skip NaN cells
        if id_type == "gene name":
            search_id = 'gene_exact:'+id
            # search in the reviewed entries
            search = search_id+reviewed
            id_matched = queryGeneName(id,search)
            if id_matched :
                new_id = id_matched
            else:
                # search in the unreviewed entries
                id_matched = queryGeneName(id,search_id)
                if id_matched :
                    new_id = id_matched
                    print (id + ' unreviewed ID ' + new_id)
                else:
                    new_id = "(no hit in human)"
                    print (id + ' no hit in human')

        else:
            # search in the reviewed entries
            id_matched = queryID(id+reviewed)
            if id_matched :
                    new_id = id_matched.group(1)
            else:
                # search in the unreviewed entries
                id_matched = queryID(id)
                if id_matched :
                    new_id = id_matched.group(1)
                    print (id + ' unreviewed ID ' + new_id)
                else:
                    new_id = "(no hit in human)"
                    print (id + ' no hit in human')
            
    else:
        new_id = "" #empty cells
       
    return new_id

def queryGeneName(id,search_id):
    try:
        query= search_id + '+organism:9606'   # search for reviewed in human
        params={
            'query': query,
            'format': 'tab',
            'columns':'id,genes(PREFERRED)',
            'sort': 'score'
        }
        response = requests.get('https://www.uniprot.org/uniprot', params = params, timeout=120)
        data = response.text
        if data:
            StringData = StringIO(data)
            df = pd.read_csv(StringData, sep ="\t")
            if df.loc[df['Gene names  (primary )'] == id].empty:
                id_matched = df.at[0,'Entry'] 
            else:
                id_matched = df.loc[df['Gene names  (primary )'] == id, 'Entry'].values[0]
        else:
            id_matched = ""

    except ConnectionError:    
        time.sleep(1)
        return queryGeneName(id,search_id)

    return id_matched

def queryID(search_id):
    try:
        query = search_id + '+organism:9606' #search for highest scored unreviewed in human
        params = {
            'query': query,
            'format' : 'list',
            'sort': 'score' 
        }
        response = requests.get('https://www.uniprot.org/uniprot', params = params, timeout=120)
        ids = str(response.content)
        id_list = ids.split('\\n')                        
        id1 = id_list[0]
        id_matched = re.search(r'b\'(\w\S{5})',id1)
    except ConnectionError:    
            time.sleep(1)
            return queryID(search_id)

    return id_matched
