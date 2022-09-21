# %%
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools
import yaml
import util
import argparse

# %%
def main():
    """
    Procedures:
        1 load Songlu's CDD token
        2 load the config for registration
        3 identify data file path(sdf)
        4 call the registration function 
        5 save the #SiTX to a csv file
    """
    parser = argparse.ArgumentParser(description='compounds registration into CDD, only for DegraderStudio projects')
    parser.add_argument('file_path', type=str)
    parser.add_argument('projects', type=str, help="Should be one of 'PLT_BRD4/PLT_KRAS_Degrader/PLT_SMARCA2/TEST")
    parser.add_argument('hypothesis', type=str, help='Design cycle #')
    parser.add_argument('creators', type=str, help='Name of the designer')
    args = parser.parse_args()
    file_path = args.file_path
    projects = args.projects
    hypothesis = args.hypothesis
    creators = args.creators

    # %%
    project_list = ['PLT_BRD4', 'PLT_KRAS_Degrader', 'PLT_SMARCA2', 'TEST']
    assert projects in project_list, 'projects name should match it in CDD and only DegraderStudio projects are allowed'

    # %%
    # instert CDD token and load it here
    with open('/CDD_API_scripts/key.yaml') as f:
        mytoken = yaml.safe_load(f)

    # %%
    # load the sdf file with smiles
    mol_df = rdkit.Chem.PandasTools.LoadSDF(file_path, idName='ID', molColName=None, includeFingerprints=False, 
                                            isomericSmiles=True, smilesName='Smiles', embedProps=False, removeHs=True, strictParsing=True)

    # %%
    mol_df, registration_dict = util.com_reg(file_path, mol_df, mytoken, projects, hypothesis, creators)
    util.save_to_csv(mol_df, registration_dict, file_path)
    print(f'HTTP 200, {len(mol_df)} molecules has been registered')

    # %%


if __name__ == "__main__":
    main()



# %%


