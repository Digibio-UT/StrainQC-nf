# Corrected API Contract — StrainQC (LOCAL source of truth)

> Produced 2026-06-25 from the LOCAL Mac copy at `/Users/marko/Documents/DTU_project/SlurmApp`.
> Supersedes the server-snapshot contract in `HANDOFF_frontend_source_of_truth.md`.
> Build the Nextflow/Seqera re-platform against THIS spec.

## Source-of-truth verdict

| Layer | Authoritative copy | Evidence |
|---|---|---|
| Frontend | **LOCAL** | FileUpload 276→**495** lines; ResultsView 1251→**2091** lines; **GenomeBrowser.tsx (680 lines) is local-only** — not in server snapshot. Local strictly larger ⇒ ahead, not behind. |
| Backend API | **LOCAL** | Local git history contains the entire Manta arc (`23708a1`…`fd068f8`) **plus the final revert `30ce748`**. Server snapshot shows mid-Manta state (2758 lines, `manta_vcf` in results); local is post-revert (2621 lines, **zero** Manta refs in `backend/app`). |
| Pipeline / executor | **LOCAL** | Same revert; `grep -ri manta backend/app` = empty. No remote configured (`git remote -v` empty) → copies sync manually; local holds the newest commits. |

No `git` remote exists, so the two checkouts can only diverge by manual copy. Local is ahead on all three.

## Drift vs. the server snapshot (handoff §"API contract as observed on the SERVER copy")

- **`manta_vcf` is GONE.** Reverted in `30ce748`. Results.json no longer emits it; the Nextflow contract must NOT include it (or mark removed). This is the single biggest backend contract change.
- **4 routes the server marked "INCOMPLETE" are now fully enumerated** (annotations GET/POST, consensus_gbk, consensus/preview, consensus/download). See below.
- **`POST /annotations/{job_id}` takes a single `FeatureEdit`**, not `AnnotationEdits{edits:[]}`. (The `AnnotationEdits` model exists but is unused by the route.)
- **`consensus/preview` default `freq_threshold=0.0`** (server guessed 0.9). `consensus/` and `consensus_gbk/` default 0.9.
- **Frontend sends `resolution=2000`** to `/coverage` (route default is 1000).
- **`/evidence/{job_id}`** is `get_gatk_evidence` — reads JSON under `gatk/` and renders HTML. breseq evidence is NOT served here; the frontend links breseq evidence as a **static file** under `/results_files/...` instead.
- **`consensus_gbk/{job_id}` (GET)** appears to be legacy/unused by the frontend — the UI uses `POST consensus/download/{job_id}` instead. Candidate to drop in the rewrite.

---

## Full backend route surface (16 routes)

Base prefix `/api/v1`. CORS `*`, credentials on. Static mount `/results_files` → `RESULTS_DIR`.
Routers: `upload` → `/api/v1/upload`, `analysis` → `/api/v1/analysis`.

### Upload
| Method | Path | Request | Response |
|---|---|---|---|
| POST | `/upload/` | multipart `file`; accepts `.fastq/.fq(.gz)`, `.fasta/.fna/.fa`, `.gbk/.gb` | `{filename, file_id, saved_as, message}`; saved as `{uuid}.{ext}` |
| GET | `/upload/files?type=reference` | `type` optional (`reference` ⇒ only gbk/gb/fasta/fna/fa) | `[{filename, size}]` |

### Analysis
| Method | Path | Request | Response |
|---|---|---|---|
| POST | `/analysis/start` | `{file_ids:[str], reference_genome_id?, description?, options?:dict}` | `{job_id, task_id, status:"queued", message}` |
| GET | `/analysis/status/{job_id}` | — | `{job_id, status, error?, steps?}`; status ∈ completed/failed/processing/unknown |
| GET | `/analysis/results/{job_id}` | — | full results.json (see schema below); seq_id version suffixes stripped |
| GET | `/analysis/evidence/{job_id}?path=` | `path` = rel path under `gatk/` (e.g. `evidence/x.json`) | `text/html` (GATK read-evidence viewer) |
| GET | `/analysis/coverage/{job_id}?resolution=1000&chrom=` | `resolution` (FE sends 2000), `chrom` optional | `{job_id, total_bases, data:[{chrom,pos,depth}]}`; cached per-chrom |
| GET | `/analysis/active_jobs` | — | `[{job_id, task_id, worker, state, started_at}]` |
| POST | `/analysis/stop/{job_id}` | — | `{message}` (revokes Celery task) |
| GET | `/analysis/genome_features/{job_id}` | — | `{genome_id, length, features:[GenomeFeature]}` (BioPython GBK parse + annotation edits) |
| GET | `/analysis/consensus/{job_id}?freq_threshold=0.9&min_qual=0` | thresholds | FASTA (`text/plain`) |
| POST | `/analysis/annotations/{job_id}` | **single `FeatureEdit`** (see model) | persisted to `annotations.json` |
| GET | `/analysis/annotations/{job_id}` | — | saved edits (`{edits:[...]}` or `[]`) |
| GET | `/analysis/consensus_gbk/{job_id}?freq_threshold=0.9&min_qual=0` | thresholds | GBK file *(legacy; FE doesn't call — drop candidate)* |
| GET | `/analysis/consensus/preview/{job_id}?freq_threshold=0.0&min_qual=0` | thresholds | `List[MutationProposal]` |
| POST | `/analysis/consensus/download/{job_id}?freq_threshold=0.9&min_qual=0` | body `List[MutationProposal]` (approved) | GBK file (blob) |

---

## Pydantic models

**`AnalysisRequest`**: `file_ids:List[str]`, `reference_genome_id?:str`, `description?:str`, `options?:dict`.
`options` is free-form dict carrying: `variant_caller("breseq"|"gatk")`, `polymorphism_mode:bool`, `ploidy?`, plus
breseq knobs (`breseq_limit_fold_coverage`, `breseq_base_quality_cutoff`, `breseq_min_mapping_quality`,
`breseq_polymorphism_frequency_cutoff`, `breseq_no_junction_prediction`, `breseq_consensus_min_cov_each_strand`,
`breseq_targeted_sequencing`) and GATK knobs (`gatk_min_base_quality`, `gatk_min_mapping_quality`,
`gatk_max_reads_per_align_start`, `gatk_dont_use_soft_clipped`, `gatk_call_conf`, `gatk_callable_depth`).

**`FeatureEdit`**: `original_id?`, `new_name?`, `new_product?`, `new_locus_tag?`, `action="edit"|"delete"|"add"`,
and for add: `start?`, `end?`, `strand?=1`, `seq_id?`.

**`MutationProposal`**: `mutation_id`, `seq_id`, `position:int`, `type`, `frequency:float`, `gene_name?`,
`locus_tag?`, `original_product?`, `mutation_detail`, `proposed_name`, `is_coding:bool`, `selected=False`,
`context_left?`, `context_right?`, `coverage_context?:List[int]`, `nearby_features?:List[dict]`.

---

## results.json schema (the GET /results payload the frontend consumes)

```
{
  job_id, status,
  variant_caller?: "breseq"|"gatk",
  polymorphism_mode?: bool,
  reference_file?,
  steps?: { <step>: <state-string> },
  mutations?: [ { type, seq_id, position, mutation, freq, annotation, gene,
                  description, details:[str], evidence_links?:[{text,url}] } ],
  read_stats?: { total_reads, mapped_percentage, unmapped_percentage, unmapped_reads },
  coverage_file?, hc_snps_file?, unmapped_files?:[str],
  mash_results?:        [ {identity, shared_hashes, median_multiplicity, p_value, query_id, comment} ],
  mash_unmapped_results?: [ <same as mash_results> ],
  contaminants?:        [ {name, rank, percentage, reads} ],
  blast_results?:       [ {title, count, identity, evalue, queries} ],
  snp_blast_results?:   [ <BlastHit> ],   snp_query_file?
  // manta_vcf  ← REMOVED in 30ce748; do NOT emit
}
```

---

## Frontend → API call map (the real contract to preserve)

axios over Vite proxy → `http://127.0.0.1:8000`; proxies `/api` and `/results_files`. No axios baseURL; paths hardcoded with `/api/v1`.

- **FileUpload.tsx** (495): `GET /upload/files?type=reference` (reads `filename`); `POST /upload/` ×N (reads `file_id`); `POST /analysis/start` (reads `job_id`; builds the `options` object above).
- **ResultsView.tsx** (2091): polls `GET /status/{id}` @2s (`status`,`steps`,`error`); `GET /results/{id}` (AnalysisResult); `GET /coverage/{id}?resolution=2000&chrom=` (`data[]`); `GET /genome_features/{id}` (`features`); `POST /annotations/{id}` (edit/delete/add single FeatureEdit) then refetch; `GET /consensus/{id}?freq_threshold&min_qual` (blob FASTA); `GET /consensus/preview/{id}?...` (proposals); `POST /consensus/download/{id}?...` (blob GBK); `GET /evidence/{id}?path=` (HTML in iframe, GATK).
  Static links: `/results_files/{id}/{gatk|breseq/data}/{file}` (unmapped), `/results_files/{id}/{hc_snps_file}` (CSV), `/results_files/{id}/breseq/output/{href}` (breseq evidence).
- **GenomeBrowser.tsx** (680, local-only): pure SVG presentation, **no API calls**; receives `features`/`variants` from ResultsView, emits `onEvidenceOpen(url)`.
- **App.tsx** (17): thin shell.

Frontend TS interfaces: `Mutation, Contaminant, BlastHit, MashHit, AnalysisResult, CoverageData, GenomeFeature` (ResultsView) and `GenomeFeature, VariantMarker, GenomeBrowserProps` (GenomeBrowser).

---

## Implications for the Nextflow re-platform

1. Contract is **frozen against local**, minus `manta_vcf`. The pipeline's final "report" process must emit `results.json` exactly as above on any backend.
2. **Interactive/on-demand routes re-read pipeline artifacts** (BAM/VCF/coverage/GBK): `coverage`, `evidence`, `genome_features`, `consensus*`. These are the egress-cost concern (open strategy Q2) — pre-compute coverage-downsample + evidence JSON + features in-pipeline so the API only serves files.
3. `consensus_gbk` (GET) is dead — drop it; keep `consensus/preview` + `consensus/download`.
4. Executor layer (`worker.py`, `slurm_runner.py`) is the part Nextflow + profiles replaces; `wrappers.py`→processes, `parsers.py`→report process. No Manta to port.
5. Pre-rewrite hygiene (carry forward, don't fix in working app): `allow_origins=["*"]`, no auth, file lookup by `startswith(file_id)`, placeholder `reference.fasta` fallback in `/start`.
