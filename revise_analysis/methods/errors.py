"""Expected missing scientific prerequisites, distinct from execution failures."""


class PrerequisiteUnavailable(ValueError):
    """The method cannot run on this declared input or support."""
