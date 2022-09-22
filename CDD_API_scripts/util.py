import pandas as pd
import requests
import json
import os
from rdkit import Chem
from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D
from rdkit.Chem.rdmolfiles import SmilesWriteParams, CXSmilesFields, MolToMolBlock
from rdkit.Chem.MolStandardize import rdMolStandardize

def uncharge(mol):
    uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
    res = uncharger.uncharge(mol)
    res.UpdatePropertyCache(strict=False)

    return res


def addEnhancedStereoAnnotations(m):
    gpLookup = {Chem.StereoGroupType.STEREO_OR:"or",
                Chem.StereoGroupType.STEREO_AND:"&",
                Chem.StereoGroupType.STEREO_ABSOLUTE:"abs",
               }
    sgs = m.GetStereoGroups()
    for i,sg in enumerate(sgs):
        typ = gpLookup[sg.GetGroupType()]
        for at in sg.GetAtoms():
            nt = ""
            if at.HasProp("atomNote"):
                nt += at.GetProp("atomNote")+","
            nt += f"{typ}{i+1}"
            at.SetProp("atomNote",nt)


def standardizeMolecule(mol):
    params = SmilesWriteParams()
    params.isomericSmiles = False
    params.kekuleSmiles = False
    params.rootedAtAtom = -1
    params.canonical = True
    params.allBondsExplicit = False
    params.allHsExplicit = False

    unchargedMol = uncharge(mol)
    addEnhancedStereoAnnotations(unchargedMol)

    return Chem.MolToCXSmiles(unchargedMol, params, CXSmilesFields.CX_ENHANCEDSTEREO)


def com_reg(file_path, protac, mytoken, projects, hypothesis, creators):
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
    for idx in range(len(protac)):
        synonym = '_'.join([synonym_root, str(idx+1)]) 
        smile = standardizeMolecule(protac[idx])
        data_package = {"molecule": {"synonyms":[synonym],
                                    "smiles":smile},
                        "projects":[projects],
                        "batch_fields": {"Hypothesis": hypothesis,
                                        "Ideated By": creators}
                        }
        response = requests.post(url, headers=headers, data=json.dumps(data_package))
        if response.status_code==200:
            print(f"{response.json()['molecule_batch_identifier']} registered")
            registration_dict[smile] = [response.json()['molecule']['name'], response.json()['name']]
        else:
            print(response.json()['error'])
    return protac, registration_dict
    

def save_to_csv(protac, registration_dict, file_path):
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
    mol_df = pd.DataFrame(registration_dict.items(), columns=['Smiles', 'molecule_batch_identifier'])
    split_df = pd.DataFrame(mol_df['molecule_batch_identifier'].tolist(), columns=['SiTX', 'Batch'])
    df = pd.concat([mol_df, split_df], axis=1)
    df = df.drop('molecule_batch_identifier', axis=1)
    df.to_csv(file_path[:-4]+'.csv', index=False)