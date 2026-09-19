# Pathway activity

Pathway activity is a side-specific gene-set score. Raw and SVC are scored on
their own native units; there is no observation pairing and no donor or
projection object in this analysis.

## Calculation

The declared GMT resource is read as a gene-set mapping. For each selected
gene set, the workflow records requested genes, genes present in the object's
gene axis, and the coverage fraction before asking the optional OmicVerse
`single.geneset_aucell` provider to calculate AUCell scores. The provider,
provider version, scorer, score key, AUC threshold, and seed are saved with the
metadata. The default workflow AUC threshold is 0.05 and the default seed is
42; the notebook keeps both visible beside the call.

The provider interprets `AUC_threshold` as a detection quantile used to derive
the rank cutoff, rather than as a final score cutoff. The accepted values are
0.01, 0.05, 0.1, 0.5, and 1.0. Metadata records the effective rank fraction
and rank cutoff, so a dense or uniformly detected input can be reported as
unavailable when no valid cutoff lies between two and `n_genes - 1`. The
workflow does not replace that failure with a different threshold or scoring
method.

Coverage is a property of the object and may differ between Raw and SVC. A
gene set with no overlap, an unavailable provider, is marked unavailable with its reason. An unexpected provider failure remains an execution error, including its traceback. The workflow does not substitute a rank
score, fill missing scores with zero, or turn coverage into a biological
validation claim. A score distribution and spatial field are descriptive
outputs that need the declared resource and coverage to be interpreted.

The notebook uses the formal sample/project gene-set parameters, with explicit overrides displayed. Synthetic DEMO_PROGRAM belongs only to the example configuration. Impact scores each side’s fixed scoring cohort once per program and derives all spatial summaries from saved unit scores; changing Region count never triggers another scoring call.

## Availability boundary

Pathway absence does not invalidate reconstruction impact or Moran results.
The batch result can be partial, and the notebook continues to the integrated
observations cell with a clear unavailable entry. The read-only report consumes
saved score tables and metadata only; it does not invoke OmicVerse or read an
H5AD during rendering.
