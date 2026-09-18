# Continuous reconstruction-impact notebook

`01_reconstruction_impact.ipynb` is the small scientific workbench for a
single Raw/SVC sample. It keeps the story visible from input declaration to
integrated observations:

1. configuration and input contract;
2. independent Raw/SVC overview;
3. independent partitions and visible intermediate labels;
4. optional Raw-defined membership comparison;
5. full Raw `Level1` anatomy context;
6. common physical windows with independent side units and draws;
7. SVC State and positive Gain regions;
8. native-graph Moran tables;
9. side-specific pathway activity;
10. an integrated reading with unavailable reasons preserved.

The source notebook is intentionally unexecuted. From this directory, the
default parameter cell looks for `../data/example/sample.yaml` and writes
exploratory tables and figures to `../output/notebook/example`. Override either
location without editing the notebook:

```bash
SAMPLE_YAML=/path/to/sample.yaml \
OUTPUT_DIR=/path/to/notebook-output \
jupyter nbconvert --to notebook --execute notebooks/01_reconstruction_impact.ipynb \
  --output 01_reconstruction_impact.executed.ipynb
```

The repository's small synthetic fixture is created by the root example
script and is suitable for checking the call sequence. A real sample must
declare the expression scale from upstream evidence; the notebook does not
infer normalization or reverse a transformation. Missing Level2 only affects
stages that request it. Missing the optional AUCell provider leaves the other
stages readable and is recorded in the integrated observations.

This notebook calls the same public methods and stage functions used by the
analysis package, but its output is exploratory evidence. The batch runner is
the owner of formal saved results and the static report; report rendering
reads those saved results and never executes notebook calculations.
