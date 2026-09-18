"""Independent post-reconstruction analysis tools."""

__version__ = "0.1.0"

def load_sample(*args, **kwargs):
    from .io import load_sample as implementation
    return implementation(*args, **kwargs)

def run_analysis(*args, **kwargs):
    from .runner import run_analysis as implementation
    return implementation(*args, **kwargs)

def run_batch(*args, **kwargs):
    from .batch import run_batch as implementation
    return implementation(*args, **kwargs)
