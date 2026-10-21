# RNA Frustratometer

A research project on local energetic frustration in RNA structures. Work in progress.

## Setup

Requires Python 3.12, uv, Node.js 24 and Bun 1.4.3.

```sh
uv sync --locked
bun install --cwd frontend --frozen-lockfile
```

## Frontend

From the repository root:

```sh
bun run --cwd frontend dev
```

Open http://localhost:5173 (or the URL printed by Vite).
Tailwind CSS is compiled automatically by Vite; no separate watcher is needed.

## Backend

Start FastAPI in a separate terminal from the repository root:

```sh
uv run --locked fastapi dev backend/app.py
```

Health check: http://127.0.0.1:8000/api/v1/health.
Interactive API docs: http://127.0.0.1:8000/docs.

## Checks

```sh
uv run pre-commit install --hook-type pre-commit --hook-type pre-push
uv run pre-commit run --all-files
uv run pre-commit run --all-files --hook-stage pre-push
```
