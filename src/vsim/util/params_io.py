import json
from typing import Dict, Any
from dataclasses import replace

from vsim.model.params import ModelParams

def load_params_json(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("Parameter file must be a JSON object (dict).")
    return data

def override_params(base: ModelParams, updates: Dict[str, Any]) -> ModelParams:
    allowed = set(base.__dataclass_fields__.keys())
    unknown = set(updates.keys()) - allowed
    if unknown:
        raise KeyError(f"Unknown parameter keys in file: {sorted(unknown)}")
    return replace(base, **updates)
