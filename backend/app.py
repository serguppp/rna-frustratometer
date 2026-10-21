from importlib.metadata import version

from fastapi import FastAPI

app = FastAPI(
    title="RNA Frustratometer API",
    version=version("rna-frustratometer"),
)


@app.get("/api/v1/health", tags=["health"])
def health() -> dict[str, str]:
    return {"status": "ok"}
