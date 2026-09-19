"""Sample configuration and independent Raw/SVC loading."""
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

@dataclass
class Sample:
    sample_id: str
    raw: Any
    svc: Any
    config: dict
    source: Path

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
        """Declared matrix identity and transformation, never inferred from .X."""
        if side not in {"raw", "svc"}:
            raise ValueError("side must be raw or svc")
        declaration = self.config.get("expression", {})
        per_side = declaration.get(side, {})
        # An older shared scale remains readable, but does not assert identity.
        return {"identity": per_side.get("identity", "unknown"),
                "scale": per_side.get("scale", declaration.get("scale", "unknown")),
                "matrix": per_side.get("matrix", "X")}

    def expression_unavailable(self, side: str) -> str | None:
        declaration = self.expression(side)
        if declaration["identity"] == "unknown":
            return f"{side}: expression matrix identity is unknown"
        if declaration["scale"] == "unknown":
            return f"{side}: expression transformation state is unknown"
        if declaration["scale"] not in {"untransformed_nonnegative", "log1p", "log1p_nonnegative"}:
            return f"{side}: unsupported expression scale {declaration['scale']!r}"
        if getattr(self, side).X is None:
            return f"{side}: expression matrix .X is absent"
        return None

    @property
    def spatial_key(self):
        return self.config.get("spatial", {}).get("key", "spatial")

    def labels(self, side: str, key: str | None = None):
        adata = getattr(self, side)
        values = adata.obs[key or self.broad_key].astype("string")
        return values.replace(self.config.get("label_aliases", {}))


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
    expression = config.get("expression", {})
    if not isinstance(expression, dict):
        raise ValueError("expression must be a mapping")
    for side in ("raw", "svc"):
        declaration = expression.get(side, {})
        if not isinstance(declaration, dict):
            raise ValueError(f"expression.{side} must be a mapping")
        for field in ("identity", "scale"):
            value = declaration.get(field, "unknown")
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"expression.{side}.{field} must be a nonempty string")
        if declaration.get("matrix", "X") != "X":
            raise ValueError("The single-object interface currently consumes matrix X only")
    return source, config


def resolve_path(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def load_sample(sample_yaml: str | Path) -> Sample:
    """Read native objects without pairing, transformation or annotation.

    Expression scale is a declared upstream contract, not inferred from values.
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
