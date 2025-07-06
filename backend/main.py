# backend/main.py
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from green_core import GreenAssess

app = FastAPI(title="Enviro-Mol API")

class MolIn(BaseModel):
    smiles: str

@app.post("/api/assess")
def assess(data: MolIn):
    try:
        ga = GreenAssess(data.smiles)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"eri": ga.eri, "scores": ga.scores, "desc": ga.descriptors}
