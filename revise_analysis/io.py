"""Sample configuration and independent Raw/SVC loading."""
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

from .analyses._shared import normalize_labels


LINEAR_EXPRESSION_SCALE = "untransformed_nonnegative"
_COMPATIBLE_LINEAR_SCALES = {"untransformed", LINEAR_EXPRESSION_SCALE}
_UNKNOWN_SCALE = "unknown"


def _expression_declaration(config: dict, side: str) -> dict:
    """Return a validated source declaration without inferring from values."""
    expression = config.get("expression", {})
    if not isinstance(expression, dict):
        raise ValueError("expression must be a mapping")
    per_side = expression.get(side, {})
    if not isinstance(per_side, dict):
        raise ValueError(f"expression.{side} must be a mapping")
    identity = per_side.get("identity", "unknown")
    if not isinstance(identity, str) or not identity.strip():
        raise ValueError(f"expression.{side}.identity must be a nonempty string")
    matrix = per_side.get("matrix", "X")
    if matrix != "X":
        raise ValueError("The single-object interface currently consumes matrix X only")
    legacy_scales = []
    for location, owner in (("expression.scale", expression),
                            (f"expression.{side}.scale", per_side)):
        if "scale" not in owner:
            continue
        scale = owner["scale"]
        if not isinstance(scale, str) or not scale.strip():
            raise ValueError(f"{location} must be a nonempty string when declared")
        if scale not in _COMPATIBLE_LINEAR_SCALES | {_UNKNOWN_SCALE}:
            raise ValueError(
                f"{location} {scale!r} is incompatible with the formal "
                "nonnegative linear X contract"
            )
        legacy_scales.append(scale)
    return {"identity": identity, "legacy_scale_unknown": _UNKNOWN_SCALE in legacy_scales,
            "matrix": matrix}


def _validate_expression_config(config: dict) -> None:
    if not isinstance(config, dict):
        raise ValueError("config must be a mapping")
    for side in ("raw", "svc"):
        _expression_declaration(config, side)


@dataclass
class Sample:
    sample_id: str
    raw: Any
    svc: Any
    config: dict
    source: Path

    def __post_init__(self):
        # Direct construction is a supported test/notebook path and must obey
        # the same formal input contract as load_sample().
        _validate_expression_config(self.config)

    @property
    def broad_key(self):
        return self.config.get("columns", {}).get("broad", "Level1")

    @property
    def subtype_key(self):
        return self.config.get("columns", {}).get("subtype", "Level2")

    @property
    def reconstruction_key(self):
        return self.config.get("columns", {}).get("reconstruction", "SVC_cluster")

    def expression(self, side: str) -> dict:
        """Return the fixed computation contract and source-declared identity."""
        if side not in {"raw", "svc"}:
            raise ValueError("side must be raw or svc")
        declaration = _expression_declaration(self.config, side)
        return {"identity": declaration["identity"],
                "scale": LINEAR_EXPRESSION_SCALE,
                "matrix": declaration["matrix"]}

    def expression_unavailable(self, side: str) -> str | None:
        if side not in {"raw", "svc"}:
            raise ValueError("side must be raw or svc")
        source = _expression_declaration(self.config, side)
        declaration = self.expression(side)
        if declaration["identity"] == "unknown":
            return f"{side}: expression matrix identity is unknown"
        if source["legacy_scale_unknown"]:
            return f"{side}: legacy expression scale is explicitly unknown"
        if getattr(self, side).X is None:
            return f"{side}: expression matrix .X is absent"
        return None

    @property
    def spatial_key(self):
        return self.config.get("spatial", {}).get("key", "spatial")

    def labels(self, side: str, key: str | None = None):
        adata = getattr(self, side)
        selected_key = key or self.broad_key
        values = adata.obs[selected_key]
        # Reconstruction labels are scientific input identities.  They are
        # copied verbatim; categorical cell-type labels use the one canonical
        # slash-to-underscore representation and retain pandas missing values.
        if selected_key == self.reconstruction_key:
            return values.copy()
        return normalize_labels(values)


def safe_segment(value: str) -> str:
    """A sample/analysis identity is one directory name, never a path."""
    if not isinstance(value, str) or not value.strip() or value in {".", ".."}:
        raise ValueError("Identity must be a nonempty directory name")
    if any(token in value for token in ("/", "\\", "\0")) or value.startswith("."):
        raise ValueError(f"Unsafe identity: {value!r}")
    return value


def read_yaml(path: str | Path) -> dict:
    path = Path(path).expanduser().resolve()
    with path.open(encoding="utf-8") as stream:
        value = yaml.safe_load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a YAML mapping: {path}")
    if value.get("schema_version") != 1:
        raise ValueError("schema_version must be 1")
    return value


def read_sample_config(sample_yaml: str | Path) -> tuple[Path, dict]:
    source = Path(sample_yaml).expanduser().resolve()
    config = read_yaml(source)
    safe_segment(config.get("sample_id"))
    files = config.get("files")
    if not isinstance(files, dict) or not all(isinstance(files.get(k), str) for k in ("raw", "svc")):
        raise ValueError("files.raw and files.svc must name H5AD files")
    _validate_expression_config(config)
    return source, config


def resolve_path(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def load_sample(sample_yaml: str | Path) -> Sample:
    """Read native objects without pairing, transformation or annotation.

    X has one fixed linear contract; identity is source-declared, never inferred.
    Scientific requirements such as coordinates are checked by their consumers.
    """
    import anndata as ad

    source, config = read_sample_config(sample_yaml)
    objects = {}
    for side in ("raw", "svc"):
        objects[side] = ad.read_h5ad(resolve_path(config["files"][side], source.parent))
        obj = objects[side]
        if not obj.obs_names.is_unique or not obj.var_names.is_unique:
            raise ValueError(f"{side}: observation and gene IDs must be unique")
    return Sample(config["sample_id"], objects["raw"], objects["svc"], config, source)


def read_gene_sets(path: str | Path, names: list[str] | None = None) -> dict[str, list[str]]:
    """Read an explicitly selected GMT; callers own missing-resource behavior."""
    sets = {}
    with Path(path).open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[0] or not any(parts[2:]):
                raise ValueError(f"Invalid GMT row {line_number} in {path}")
            if parts[0] in sets:
                raise ValueError(f"Duplicate GMT gene set {parts[0]!r} at row {line_number}")
            sets[parts[0]] = list(dict.fromkeys(g for g in parts[2:] if g))
    if names is not None:
        if not isinstance(names, list) or not names:
            raise ValueError("gene_set_names must be a nonempty list")
        missing = set(names) - sets.keys()
        if missing:
            raise ValueError(f"Gene sets absent from {Path(path).name}: {sorted(missing)}")
        sets = {name: sets[name] for name in names}
    if not sets:
        raise ValueError(f"No gene sets in {path}")
    return sets
