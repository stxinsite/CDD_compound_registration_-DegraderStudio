import pandas as pd
import requests
import json
import os


def com_reg(file_path, mol_df, mytoken, projects, hypothesis, creators):
    """
    register compounds into CDD, usage only for Psivant
    Args:
        file_path(str): path of the sdf file
        mol_df(dataframe): load sdf with smiles
        mytoken: open the CDD token file and close after using
        projects: 'PLT_BRD4', 'PLT_KRAS_Degrader', 'PLT_SMARCA2', 'TEST'
        hypothesis: Design cycle number
        creators: Name of the designer

    Returns:
        mol_df(dataframe): smiles and SiTX#
        registration_dict(dictionary): dictionary of index and smiles
        
    Raises:
        ValueError: Error message raising from requests
    """
    registration_dict = {}
    base_url = "https://app.collaborativedrug.com/api/v1/vaults/4686/"
    url = base_url + "batches"
    token_dic = {}
    token_dic['X-CDD-token'] = mytoken
    headers = token_dic
    file_name = os.path.split(file_path)[-1]
    synonym_root = file_name[:-4].strip().replace(' ','_')
    for idx in mol_df.index.tolist():
        synonym = '_'.join([synonym_root, str(mol_df.iloc[idx]['ID'])]) 
        data_package = {"molecule": {"synonyms":[synonym],
                                     "smiles":mol_df['Smiles'][idx]},
                        "projects":[projects],
                        "batch_fields": {"Hypothesis": hypothesis,
                                         "Ideated By": creators}
                        }
        response = requests.post(url, headers=headers, data=json.dumps(data_package))
        if response.status_code==200:
            print(f"{response.json()['molecule_batch_identifier']} registered")
            registration_dict[mol_df.iloc[idx]['Smiles']] = [response.json()['molecule']['name'], response.json()['name']]
        else:
            print(response.json()['error'])
    return mol_df, registration_dict
    

def save_to_csv(mol_df, registration_dict, file_path):
    """
    save the registered #SiTX to csv
    Args:
        mol_df(dataframe): pandas dataframe after reading a 
        registration_list(list): record the time point when the session starts
        file_path(str): open the data file and close after using
    Returns:
        A csv in the same path with sdf
    Raises:
        ValueError: Some compounds in the mol_df list are not successfully registered
    """
    if len(mol_df)!=len(registration_dict):
        print('Some compounds in the list are not successfully registered.')   
    mol_df['molecule_batch_identifier'] = mol_df['Smiles'].map(registration_dict)
    split_df = pd.DataFrame(mol_df['molecule_batch_identifier'].tolist(), columns=['SiTX', 'Batch'])
    df = pd.concat([mol_df, split_df], axis=1)
    df = df.drop('molecule_batch_identifier', axis=1)
    df.to_csv(file_path[:-4]+'.csv', index=False)