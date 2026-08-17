# StrainQC-nf frontend

React + Vite web UI for StrainQC-nf. It talks to the **launcher** (`../launcher`) over HTTP —
upload reads + reference, start a run, then view mutations, coverage, the genome browser,
contamination and evidence. This is the same UI as the original app, reused unchanged apart
from small fixes.

## Run it (development)

```bash
cd frontend
npm install
npm run dev          # Vite dev server, prints a http://localhost:5173 URL
```

The dev server proxies `/api` and `/results_files` to `http://127.0.0.1:8000` (see
`vite.config.ts`), i.e. **the launcher must be reachable at `localhost:8000`**:

- Launcher running locally → nothing else to do.
- Launcher on an HPC/remote host → forward its port first, e.g.
  `ssh -N -L 8000:localhost:8000 <user>@<host>` (VS Code Remote-SSH auto-forwards it too).

Then open the printed URL. Quick check the launcher is reachable: `curl -s localhost:8000/`
should return the launcher's JSON.

## Production build

Dev mode is convenient but slower for the genome browser. For a snappy build:

```bash
npm run build        # outputs dist/
npm run preview      # serves dist/ locally
```

`Dockerfile` + `nginx.conf` build a static nginx image serving `dist/` — useful for a hosted
deployment. (A future option is to have the launcher serve `dist/` directly, so there's no
separate frontend server; not wired yet.)

## Note for remote setups

A few result links (unmapped reads, hc-SNP CSV) are currently hardcoded to
`http://localhost:8000/...`, so they rely on the launcher being at `localhost:8000` (the
tunnel above satisfies this). Everything else uses the relative `/api` + `/results_files`
paths through the Vite proxy.
