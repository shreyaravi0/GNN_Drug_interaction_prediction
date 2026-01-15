from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from predictor import predict

app = FastAPI(title="Drug Interaction API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class DrugPair(BaseModel):
    smiles1: str
    smiles2: str

@app.post("/predict")
def predict_interaction(data: DrugPair):
    result = predict(data.smiles1, data.smiles2)
    if result is None:
        return {"error": "Invalid SMILES"}
    return result
