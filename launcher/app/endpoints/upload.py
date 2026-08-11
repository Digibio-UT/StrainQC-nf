"""Upload endpoints — same contract as the old app (backend-agnostic)."""
import os
import shutil
import uuid

from fastapi import APIRouter, File, HTTPException, UploadFile

from app.config import settings

router = APIRouter()


@router.post("/")
async def upload_file(file: UploadFile = File(...)):
    if not file.filename.endswith(settings.ALLOWED_UPLOAD_EXT):
        raise HTTPException(status_code=400, detail="Invalid file format.")

    settings.ensure_dirs()
    file_id = str(uuid.uuid4())
    ext = file.filename.split(".")[-1]
    if file.filename.endswith(".fastq.gz"):
        ext = "fastq.gz"
    elif file.filename.endswith(".fq.gz"):
        ext = "fq.gz"

    saved_as = f"{file_id}.{ext}"
    dest = settings.UPLOAD_DIR / saved_as
    try:
        with open(dest, "wb") as buf:
            shutil.copyfileobj(file.file, buf)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Could not save file: {e}")

    return {"filename": file.filename, "file_id": file_id,
            "saved_as": saved_as, "message": "File uploaded successfully"}


@router.get("/files")
async def list_files(type: str = None):
    files = []
    # Curated reference library first (friendly names), then uploads. Skip the library for
    # non-reference listings.
    dirs = []
    if type == "reference":
        dirs.append(settings.REFERENCE_DIR)
    dirs.append(settings.UPLOAD_DIR)

    seen = set()
    for d in dirs:
        if not d.exists():
            continue
        for p in sorted(d.iterdir()):
            if p.name.startswith(".") or not p.is_file() or p.name in seen:
                continue
            if type == "reference" and not p.name.endswith((".gbk", ".gb", ".fasta", ".fna", ".fa")):
                continue
            seen.add(p.name)
            files.append({"filename": p.name, "size": os.path.getsize(p)})
    return files
