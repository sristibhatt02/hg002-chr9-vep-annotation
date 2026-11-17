# Teaching-example validation

Five public Ensembl chromosome 9 examples provide a small regression dataset for reference checks, anchored deletion reconciliation and transcript annotation. They are separate from the HG002 study totals; this directory is a validation component of the same project.

```sh
python3 validation/run.py fetch
python3 validation/run.py report
python3 validation/run.py verify
python3 -m unittest discover -s validation/tests -v
```

The fixed inputs are in `data/training-variants.grch38.vcf`. All five REF intervals must match GRCh38. Authentic Ensembl REST output, hashes and reference evidence are acquired locally under `results/`; the validation HTML is generated under `report/`. Those generated files are excluded from Git. No synthetic fixtures are used as scientific results.

The observed validation snapshot contains five accounted-for inputs and 216 transcript rows. COQ4 examples on ENST00000300452.8 contrast missense and intronic consequences. The finished project report is [here](../reports/experiment-report.html).
