# RNA Frustratometer

Python package for studying local energetic frustration in RNA structures.
Work in progress. The web application lives in
[rna-frustratometer-web](https://github.com/serguppp/rna-frustratometer-web).

## Development

Requires Python 3.12 and uv.

```sh
uv sync --locked
uv run pre-commit install --hook-type pre-commit --hook-type pre-push
uv run pre-commit run --all-files
uv run pre-commit run --all-files --hook-stage pre-push
uv build
```

## Branches

- `main`: shared package setup and structure-loading utilities.
- `add-frustration-briq`: the current BRiQ-based algorithm and visualizations.
- `add-frustration-cgrnasp`: development branch for fixed-geometry cgRNASP scoring.

The cgRNASP branch starts from the BRiQ implementation; cgRNASP is not implemented
yet. Both branches retain the shared development history.

BRiQ analysis requires the custom `briq:memory` Docker image in addition to the
Python dependencies. It is not included in this repository.
