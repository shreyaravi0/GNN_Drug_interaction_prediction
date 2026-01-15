# predictor.py
import torch, base64, io, requests
from rdkit import Chem
from rdkit.Chem import Draw
from torch_geometric.data import Batch
from model_def import InteractionPredictor
from utils import mol_to_graph_data_obj

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = InteractionPredictor(hidden_channels=128)
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
model_path = os.path.join(BASE_DIR, "gnn_model.pth")
checkpoint = torch.load(model_path, map_location=device)
model.load_state_dict(checkpoint["model_state_dict"])
model.to(device)
model.eval()

def mol_image_base64(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode()

def get_iupac(smiles):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/JSON"
    try:
        r = requests.get(url, timeout=5)
        return r.json()["PropertyTable"]["Properties"][0]["IUPACName"]
    except:
        return "Unknown molecule"

def get_cid(smiles):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON"
    try:
        r = requests.post(url, data={"smiles": smiles}, timeout=5)
        return r.json()["IdentifierList"]["CID"][0]
    except:
        return None

def get_description(cid):
    if not cid: return "No description available."
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/description/JSON"
    try:
        r = requests.get(url, timeout=5)
        data = r.json()
        if "InformationList" in data and "Information" in data["InformationList"]:
            infos = data["InformationList"]["Information"]
            for info in infos:
                if "Description" in info:
                    return info["Description"]
        return "No description available."
    except:
        return "No description available."

def get_drug_info(smiles):
    cid = get_cid(smiles)
    return {
        "smiles": smiles,
        "iupac": get_iupac(smiles),
        "description": get_description(cid),
        "image": mol_image_base64(smiles)
    }

def predict(smiles1, smiles2):
    d1 = mol_to_graph_data_obj(smiles1)
    d2 = mol_to_graph_data_obj(smiles2)
    if d1 is None or d2 is None:
        return None

    b1 = Batch.from_data_list([d1]).to(device)
    b2 = Batch.from_data_list([d2]).to(device)

    with torch.no_grad():
        prob = torch.sigmoid(model(b1, b2)).item()

    return {
        "drug1": get_drug_info(smiles1),
        "drug2": get_drug_info(smiles2),
        "prediction": {
            "label": "Interacting" if prob > 0.5 else "Non-interacting",
            "probability": prob
        }
    }
