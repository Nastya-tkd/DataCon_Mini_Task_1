import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import time
from tqdm import tqdm


def get_pubchem_data(compound_name):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

   
    cid_url = f"{base_url}/compound/name/{compound_name}/cids/JSON"
    try:
        response = requests.get(cid_url, timeout=30)
        response.raise_for_status()
        cid = response.json()['IdentifierList']['CID'][0]
    except Exception as e:
        print(f"Error getting CID for {compound_name}: {str(e)}")
        return None

   
    pug_view_url = f"{base_url}/compound/cid/{cid}/JSON"
    try:
        response = requests.get(pug_view_url, timeout=30)
        response.raise_for_status()
        data = response.json()

        
        result = {
            'compound_name': compound_name,
            'pubchem_cid': cid,
            'canonical_smiles': None,
            'isomeric_smiles': None,
            'inchi_key': None,
            'molecular_weight': None,
            'iupac_name': None,
            'source': 'PubChem'
        }

        
        if 'PC_Compounds' in data and len(data['PC_Compounds']) > 0:
            compound = data['PC_Compounds'][0]

            
            if 'props' in compound:
                for prop in compound['props']:
                    urn = prop.get('urn', {})
                    label = urn.get('label', '')
                    name = urn.get('name', '')

                    if label == 'SMILES':
                        if name == 'Canonical':
                            result['canonical_smiles'] = prop.get('value', {}).get('sval')
                        elif name == 'Isomeric':
                            result['isomeric_smiles'] = prop.get('value', {}).get('sval')
                    elif label == 'InChIKey':
                        result['inchi_key'] = prop.get('value', {}).get('sval')
                    elif label == 'Molecular Weight':
                        result['molecular_weight'] = prop.get('value', {}).get('fval')
                    elif label == 'IUPAC Name':
                        result['iupac_name'] = prop.get('value', {}).get('sval')

        return result

    except Exception as e:
        print(f"Error getting properties for {compound_name}: {str(e)}")
        return None



drug_names = [
    "Aspirin", "Ibuprofen", "Paracetamol", "Metformin", "Atorvastatin",
    "Omeprazole", "Sertraline", "Warfarin", "Diazepam", "Lisinopril",
    "Simvastatin", "Amoxicillin", "Ciprofloxacin", "Prednisone", "Hydrochlorothiazide"
]


results = []
for name in tqdm(drug_names, desc="Fetching data"):
    data = get_pubchem_data(name)
    if data:
        results.append(data)
    time.sleep(1)  # Ограничение запросов


df = pd.DataFrame(results)
print(f"\nSuccessfully retrieved data for {len(df)} compounds")



def canonicalize_smiles(smiles):
    if pd.isna(smiles) or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None


if not df.empty:
    df['canonical_smiles'] = df['canonical_smiles'].apply(canonicalize_smiles)
    df['isomeric_smiles'] = df['isomeric_smiles'].apply(canonicalize_smiles)


    
    def calculate_mw(smiles):
        if pd.isna(smiles):
            return None
        try:
            mol = Chem.MolFromSmiles(smiles)
            return Descriptors.MolWt(mol) if mol else None
        except:
            return None


    if 'molecular_weight' not in df.columns or df['molecular_weight'].isna().any():
        df['molecular_weight'] = df['canonical_smiles'].apply(calculate_mw)


output_file = "pubchem_compounds.csv"
df.to_csv(output_file, index=False)
print(f"Data saved to {output_file}")

if not df.empty:
    print("\nFirst few records:")
    print(df.head())
else:
    print("\nNo data was retrieved. Please check your internet connection and try again.")
