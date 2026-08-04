"""
StrainQC-nf launcher — thin FastAPI service that re-exposes the StrainQC API contract,
backed by `nextflow run` instead of Celery/Slurm-over-SSH. Runs on the HPC (or any host
that can launch the pipeline). The React frontend talks to this unchanged.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from app.config import settings
from app.endpoints import analysis, evidence, upload

settings.ensure_dirs()

app = FastAPI(title=settings.PROJECT_NAME, openapi_url=f"{settings.API_V1_STR}/openapi.json")

# TODO(auth): allow_origins=["*"] is fine for a single-lab instance, NOT a shared one.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_credentials=True,
    allow_methods=["*"], allow_headers=["*"],
)

# Serve pipeline outputs (evidence HTML, coverage, unmapped, hc-SNP CSV) like the old app.
app.mount("/results_files", StaticFiles(directory=str(settings.RESULTS_DIR)), name="results_files")

app.include_router(upload.router, prefix=f"{settings.API_V1_STR}/upload", tags=["upload"])
app.include_router(analysis.router, prefix=f"{settings.API_V1_STR}/analysis", tags=["analysis"])
app.include_router(evidence.router, prefix=f"{settings.API_V1_STR}/analysis", tags=["analysis"])


@app.get("/")
def root():
    return {"message": "StrainQC-nf launcher", "pipeline_dir": str(settings.PIPELINE_DIR)}
