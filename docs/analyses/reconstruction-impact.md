# Reconstruction impact

This analysis asks how the SVC object changes representation and local spatial
diversity relative to Raw. Raw and SVC are loaded as independent objects. Their
unit IDs, order, gene sets, and native coordinate graphs may differ. The
analysis does not infer a correspondence from position or from row order.

## Scientific sequence

The notebook exposes the stages in this order:

1. read the sample declaration and inspect the two objects;
2. partition each side on its own expression matrix;
3. optionally compare memberships on a concrete Raw-defined cohort;
4. retain the full Raw `Level1` field as anatomy context;
5. place both sides in the same physical window grid while keeping their unit
   membership and rarefaction draws independent;
6. summarize local diversity and classify State and Gain separately.

The continuous notebook exposes these stages as separate cells. Its calls are
`load_sample`, `prepare_side`, `partition_scope`, optional
`compute_partitions` plus `compare_membership`, `raw_anatomy_context`,
`native_window_diversity`, `select_region_threshold` plus
`flag_region_windows`, the native Moran wrapper, and the AUCell wrapper. It
does not call `reconstruction_impact.run` from the notebook, so intermediate
objects, parameters, plots, and unavailable reasons remain visible beside the
calculation that produced them.

Partitioning is a representation diagnostic. The working copy may be library
size normalized and log transformed for PCA, neighbors, and Leiden; the input
H5AD is never changed. The resolution, feature count, seed, window size,
minimum support, and number of draws are explicit parameters. A partition
summary records the units and genes actually retained so a small fixture does
not look like a full analysis.

Membership change is the only paired component. It starts from the Raw parent
cohort, intersects the SVC IDs for that parent, recomputes the SVC partition on
those IDs, and applies Hungarian matching to the two label vectors. The
assignment table is therefore a statement about that cohort and comparison
edge. It must not be used to redefine the native SVC population or to make
Moran, pathway, and local-diversity analysis require paired IDs. An
empty intersection or unavailable membership partition is reported as
unavailable.

The notebook leaves this comparison disabled by default. Setting
`RUN_MEMBERSHIP=1` opts into the explicit Raw-defined edge; otherwise all
independent stages continue without a pairing request.

The saved same-resolution comparison uses the configured Leiden resolution on
both sides. A separate optional `matched_k_resolutions` sweep can seek an SVC
resolution whose cluster count exactly matches the Raw parent count. Its sweep,
selected controlled-K comparison, and any no-exact-K result are saved under
`membership_matched_k_*`; these artifacts do not redefine the same-resolution
comparison.

Anatomy is a full Raw `Level1` context table. It is not inferred from Level2,
from a partition, or from the selected local subset. A window may be labelled
Tumor, Normal, Interface, or Other from the configured broad labels; absent
labels leave the anatomy stage unavailable rather than silently changing the
meaning of the field.

## Local diversity and regions

The two sides use the full Raw coordinate origin and the configured
`spatial.microns_per_coordinate` scale with one physical window side length, so
a window ID refers to the same coordinate grid. Each side contributes its own
native units to that grid and draws `min_units` units independently. The
workflow API accepts the sampled side as `adata` and keeps `n_window_draws`
and the side-specific random seed explicit. For a draw with
class proportions (p_k), the reported quantities are

\[
H = -\sum_k p_k\log p_k,\qquad
N_{\mathrm{eff}} = e^H,\qquad
E = N_{\mathrm{eff}}/K_{\mathrm{obs}}.
\]

Windows below the support threshold retain their coordinates and unit count,
but their diversity values are unavailable. Gain is computed only where both
sides have valid metrics, as

\[
\Delta N_{\mathrm{eff}} = N_{\mathrm{eff,SVC}} - N_{\mathrm{eff,Raw}}.
\]

State and Gain answer different questions. State uses the SVC `N_eff` field
from each non-`All` SVC parent Leiden partition and an accepted stable threshold
to identify high-diversity SVC windows. Gain uses
the positive part of the cross-side delta after the common window grid has been
formed. A threshold that does not meet the support or bootstrap stability rule
is recorded as `no_stable_threshold`; it is not replaced by a percentile chosen
after looking at the map.

When an upstream Raw `Level2` column is present, the workflow also saves a Raw
Level2 local-diversity baseline for each Raw parent and a common-window table
beside the corresponding native SVC parent Leiden field. The saved
`descriptive_delta_neff` compares two different label systems, so it is a
reference for reading rather than a Gain calculation or a region-selection
input. Missing Raw Level2 makes only these baseline artifacts unavailable.

The notebook may display a continuous field and a support map even when no
stable State or Gain region can be declared. Region summaries must therefore
carry a status and a reason. They are descriptive evidence, not a claim that
SVC improves biology.

## Defaults shown in the notebook

The current workflow exposes these starting values beside the calls: random
seed 42 (with a distinct SVC seed for independent sampling), Leiden resolution
0.5, up to 2,000 partition features, six native spatial neighbors for Moran,
200 local diversity draws, four minimum units per window, and 500 threshold
bootstrap draws where the threshold stage is available. The notebook uses an
explicit 32-micron fixture window while the batch API default is 100 microns;
real analyses should set `window_side_microns` deliberately. It is converted
through `spatial.microns_per_coordinate` from the sample declaration.
These are parameters, not hidden constants, and changing them changes the
scientific question and the saved provenance.

## Availability boundary

Missing Level2 does not block the reconstruction-impact story because anatomy
uses full Raw `Level1`. Only a stage that explicitly requests Level2 is
unavailable. Missing optional pathway providers likewise leaves the impact
tables and plots usable. The notebook reports the provider or input reason and
does not write zero scores or placeholder regions.

The batch runner owns published tables, `result.json`, and the read-only HTML
report. The exploratory notebook writes only its requested output directory;
its intermediate tables and figures are evidence for reading and debugging,
not a second formal batch result.
