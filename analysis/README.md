# Annotation pipeline

This directory implements the HG002 Chromosome 9 Variant Annotation. Start with the [project README](../README.md) or the [complete report](../reports/experiment-report.html).

## Scope and inputs

The completed experiment is HG002 GRCh38 chr9:128,000,001–129,000,000, pinned to GIAB NIST v4.2.1 and its matching high-confidence BED. The interval was fixed before examining annotation consequences. Its 1,307 included alleles yielded 38,589 transcript annotations. The ClinVar overlay is pinned to 2026-11-08. The whole chromosome was not annotated.

`scope.json` records the coordinates, eligibility rules and bounded REST budget. Source paths, hashes and acquisition timestamps are retained in the locally generated manifests. The original acquired manifests retain the paths and settings actually used at acquisition; moving the project does not rewrite historical evidence.

## Environment

- Python 3.10+; the observed run used Python 3.13.3 on Apple Silicon with 8 logical CPUs and 8 GiB RAM.
- `requirements.txt`: pysam 0.23.3, providing embedded bcftools/HTSlib 1.21 functionality.
- Report generation is optional and uses `requirements-report.txt`; PDF generation needs the matching Chromium browser.
- Browser checks use Playwright; Node is only needed for optional JavaScript syntax checks.

Use `python -m pip` through the active virtual environment. Local environments and large acquired inputs are excluded from Git.

## Commands

Run from the repository root:

```sh
python3 -m venv analysis/.venv
analysis/.venv/bin/python -m pip install -r analysis/requirements.txt

analysis/.venv/bin/python analysis/pipeline.py acquire
analysis/.venv/bin/python analysis/pipeline.py prepare
analysis/.venv/bin/python analysis/pipeline.py annotate
analysis/.venv/bin/python analysis/clinvar.py acquire
analysis/.venv/bin/python analysis/clinvar.py overlay
analysis/.venv/bin/python analysis/pipeline.py analyze
analysis/.venv/bin/python analysis/pipeline.py render
analysis/.venv/bin/python analysis/pipeline.py verify
```

Acquisition uses indexed VCF extraction; full remote VCF downloads and full-file hashes are not claimed. Reference FASTA and the matching BED are cached locally. Annotation uses no more than 100 alleles per request and a declared 2,000-allele pilot limit. Do not expand this REST workflow to a whole chromosome.

The `analyze`, `render` and `verify` commands make no network requests. They require an authentic successful local snapshot. A fresh clone includes source and finished shareable reports, while scientific input caches and per-allele explorer data must be acquired again.

## Scientific rules

- Split only ALT alleles present in the sample genotype; retain original records and ALT indexes.
- Require supported ACGT sequence with REF/ALT lengths no greater than 50 bases.
- Require the complete anchored REF span inside both the analyzed interval and a high-confidence region before and after normalization.
- Verify reference bases against GRCh38; normalize with bcftools and deduplicate on chromosome/position/REF/ALT.
- Keep every reported transcript annotation; regulatory and other feature annotations remain separate.
- Select a representative per allele–gene pair: MANE Select, canonical, then accession order, with ties and fallback counts visible.
- Compare term sets and the pinned Ensembl molecular-consequence ranking; do not interpret molecular severity as pathogenicity.
- Match ClinVar only after consistent normalization and exact allele identity. Retain review status, conditions, conflicts and unmatched counts.

See the [data dictionary](DATA_DICTIONARY.md) for each export, coordinate convention and denominator.

## Safe snapshots and outputs

A successful annotation snapshot records raw chunk hashes and input/settings/reference fingerprints. Chunks resume only when these match. Complete snapshot promotion is atomic; immutable prepared input copies remain with the snapshot. Analytical exports are staged before promotion, with the prior outputs retained locally for recovery. Avoid concurrent pipeline invocations.

Generated local evidence includes allele summaries, full transcript rows, representative selections, gene summaries, consequence spectra, 25 kb bins, source-allele mappings, exclusions, cases and clinical matches. The browser loads a compact overview and fetches transcript evidence locally only after allele selection.

The public report files are under `reports/`. The generated explorer is under `analysis/report/` after reproduction. Its default 2D view explains filters and uses a fixed count scale; the optional 3D component bundles Three.js 0.180.0 with its licence.

## Rebuild the shareable reports

First acquire/rebuild both the main analysis and the small teaching validation:

```sh
python3 validation/run.py fetch
python3 validation/run.py report
python3 validation/run.py verify
```

Install the optional report dependencies in a suitable environment, then:

```sh
python -m pip install -r analysis/requirements-report.txt
python -m playwright install chromium
MPLCONFIGDIR=/tmp/vep-report-matplotlib python analysis/build_experiment_report.py
python analysis/render_pdf.py
```

The HTML embeds its CSS and five vector scientific figures, including two 3D views. The builder also exports the figure SVGs for the project README. The PDF renderer works offline. Source assertions and run scope are explicit throughout. External links are optional citations, not rendering dependencies.

## Verification

```sh
analysis/.venv/bin/python -m unittest discover -s analysis/tests -v
python3 -m unittest discover -s validation/tests -v
analysis/.venv/bin/python analysis/pipeline.py verify
```

To exercise the explorer, start the local server from the root and run the browser check in an environment with Playwright/Chromium:

```sh
python3 -m http.server 8000 --bind 127.0.0.1
# In another terminal:
python analysis/tests/browser_check.py
```

The server URL is http://127.0.0.1:8000/analysis/report/index.html. Browser checks verify automatic filters, combined/reset/empty states, region and allele selection, valid zero-transcript results, keyboard focus, mobile layout and WebGL failure. Chromium was checked; the matching Safari/WebKit executable was unavailable. No formal accessibility certification is claimed.

## Resource limits

Preparation took approximately 0.21 seconds and 97 MB peak RSS excluding downloads; 14 sequential REST batches took approximately 440 seconds. These observations do not predict offline chromosome throughput.

The standard Ensembl 116 GRCh38 cache archive is 27,644,657,162 bytes compressed, larger than initially available storage before extraction. Docker's daemon was unavailable and no local VEP/cache was found. Published images support linux/arm64, so architecture is not the blocker. A provisioned offline backend, streaming adapter and chromosome-scale verification remain future work.

Raw upstream snapshots remain local pending redistribution review. Scientific releases, source references and dependency versions are retained for reproducibility. Nothing in this workflow creates a remote repository, pushes or publishes the project.
