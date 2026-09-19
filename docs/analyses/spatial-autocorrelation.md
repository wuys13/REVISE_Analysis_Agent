# Spatial autocorrelation

Spatial autocorrelation is calculated per side on native units. Raw and SVC
each receive their own graph from their own coordinates; a shared observation
ID list or shared graph is not an implicit prerequisite. A comparison of two
Moran tables is a separate, explicit gene-level reading step.

## Calculation

The default graph is a deterministic symmetric binary six-nearest-neighbor
graph with self edges removed. Ties follow input order. For each gene, the
calculation works on a copy of the object, requires the fixed finite,
nonnegative, unlogged linear `.X` contract, applies library-size normalization
to a target sum of 10,000 and `log1p`, and then computes

\[
I = \frac{n}{W}\frac{\sum_{ij}w_{ij}(x_i-\bar{x})(x_j-\bar{x})}
                 {\sum_i(x_i-\bar{x})^2},
\]

where (W=\sum_{ij}w_{ij}). The input `.X` remains unchanged. Results retain
the gene ID, Moran value, unit count, edge count, status, and reason.

Genes with constant expression have no defined denominator and are reported as
`not_computable`. Too few units or a graph without edges is reported as
`insufficient_support`. Non-finite coordinates or expression values are input
errors. These states are meaningful output, not zeros.

The workflow default is six neighbors and seed 42 for deterministic sampling
when a sample-size limit is requested. The notebook displays the native Raw
and SVC tables and their distributions independently. It does not call a
paired Moran implementation or claim that a difference is a biological gain.
The notebook-facing call is `spatial_autocorrelation.compute_moran` once per
side; the full objects are used for this stage even when the partition view
uses an exploratory unit sample.

## Full and shared gene views

The default per-side tables retain each object's available genes. If a later
comparison is restricted to shared gene IDs, the intersection and its
denominator must be recorded in that comparison table. Shared genes do not
justify sharing a spatial graph: spatial support is still native to each
object.

## Availability boundary

Missing coordinates, duplicate IDs, insufficient support, constant genes, and
missing optional dependencies are reported with a status and reason. A saved
batch result may therefore be `partial` while the other side or another
analysis remains usable. The static report reads these saved tables; it never
opens H5AD files or recomputes Moran's I.

Impact runs this native calculation per configured scope, retains each side’s full gene table, and writes a separate shared-gene comparison with actual unit and gene denominators. Unknown expression identity leaves the component unavailable; legacy unknown scale remains unavailable, and old log-scale declarations are rejected.
