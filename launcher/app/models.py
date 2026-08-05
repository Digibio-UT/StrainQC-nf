"""Pydantic models shared by the launcher endpoints (ported from the old app's
models/annotations.py + models/consensus.py)."""
from typing import List, Optional

from pydantic import BaseModel


class FeatureEdit(BaseModel):
    original_id: Optional[str] = None       # existing feature id; None for a new add
    new_name: Optional[str] = None
    new_product: Optional[str] = None
    new_locus_tag: Optional[str] = None
    action: str = "edit"                     # "edit" | "delete" | "add"
    # For action == "add":
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[int] = 1
    seq_id: Optional[str] = None


class MutationProposal(BaseModel):
    mutation_id: str
    seq_id: str
    position: int
    type: str                                # 'SNP' | 'DEL' | 'INS' | ...
    frequency: float
    gene_name: Optional[str] = None
    locus_tag: Optional[str] = None
    original_product: Optional[str] = None
    mutation_detail: str
    proposed_name: str
    is_coding: bool
    selected: bool = False
    context_left: Optional[str] = None
    context_right: Optional[str] = None
    coverage_context: Optional[List[int]] = None
    nearby_features: Optional[List[dict]] = None
