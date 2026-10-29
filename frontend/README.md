# StrainQC-nf frontend

React + Vite web UI for StrainQC-nf. It talks to the **launcher** (`../launcher`) over HTTP —
upload reads + reference, start a run, then view mutations, coverage, copy number, the genome
browser, contamination and evidence. This is the same UI as the original app, reused unchanged
apart from small fixes.

## Requirements

- **Node.js 20.19+ or 22.12+** (Vite 7 does not run on older versions; check with
  `node --version`). npm comes with Node.
- A running launcher (the backend) reachable at `localhost:8000`, see below.

There is no `requirements.txt` here: that is the launcher's (Python) file. For the frontend,
`package.json` lists the dependencies and `package-lock.json` pins their exact versions;
`npm ci` installs exactly those.

## Run it (development)

```bash
cd frontend
npm ci               # install the pinned dependencies (first time, or after a pull)
npm run dev          # Vite dev server, prints a http://localhost:5173 URL
```

The dev server proxies `/api` and `/results_files` to `http://127.0.0.1:8000` (see
`vite.config.ts`), i.e. **the launcher must be reachable at `localhost:8000`**. Quick check:
`curl -s localhost:8000/` should return the launcher's JSON.

## Connecting to a launcher

| Where the launcher runs | What to do |
|---|---|
| On this machine | nothing |
| On an HPC login node (e.g. UT HPC) | forward its port: `ssh -N -L 8000:localhost:8000 <user>@<host>` (VS Code Remote-SSH forwards it too) |
| On a cloud VM in Seqera mode (e.g. DTU BRIGHT) | the same tunnel to the VM: `ssh -N -i <key> -L 8000:127.0.0.1:8000 <user>@<vm>`; the VM's firewall must allow SSH from your address |

The UI is the same in every case. In Seqera mode (`STRAINQC_LAUNCH_MODE=seqera`, see
`../launcher/app/config.py`) a submitted job becomes a run in the site's Seqera workspace, and
its results appear in the UI when the run finishes.

## Production build

Dev mode is convenient but slower for the genome browser. For a snappy build:

```bash
npm run build        # outputs dist/
npm run preview      # serves dist/ locally
```

`Dockerfile` + `nginx.conf` build a static nginx image serving `dist/` and proxying `/api` and
`/results_files` to the launcher.

## Hosting the built UI next to the launcher (users need no Node.js)

On a VM that already runs the launcher on `127.0.0.1:8000`, nginx can serve `dist/` and proxy
to the launcher, so users need only a browser and one SSH tunnel:

1. Build on any machine with Node (`npm run build`) and copy `dist/` to the VM.
2. Install nginx and add a site based on `nginx.conf`, with `listen 127.0.0.1:8080;`, `root`
   pointing at the copied `dist/`, and both `proxy_pass` targets set to `http://127.0.0.1:8000`.
   Keep `client_max_body_size` (FASTQ uploads), and disable nginx's default site so nothing
   listens publicly.
3. Users: `ssh -N -L 8080:127.0.0.1:8080 <user>@<vm>`, then open http://localhost:8080.

To update the UI, rebuild and copy `dist/` again; the launcher does not need a restart.
