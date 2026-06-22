import yaml

REQUIRED_KEYS = ["input_path", "file_type", "functional", "damping", "run_settings", "output_unit"]

def validate_required(cfg):
    missing = [key for key in REQUIRED_KEYS if key not in cfg]
    if missing:
        raise KeyError(f"Missing required keys in YAML: {missing}")

def apply_defaults(cfg):
    # Ensure run_settings exists
    if "run_settings" not in cfg or cfg["run_settings"] is None:
        cfg["run_settings"] = {}
    rs = cfg["run_settings"]
    rs.setdefault("pairwise", True)
    rs.setdefault("triples", True)
    rs.setdefault("c11", False)
    rs.setdefault("verbose_conv", True)
    rs.setdefault("extrapolate_triples", True)
    rs.setdefault("output_unit", "Ry")
    rs.setdefault("r2tol", 0.99)
    return cfg

def normalize(cfg):
    # Ensure input_path is always a list
    if isinstance(cfg["input_path"], str):
        cfg["input_path"] = [cfg["input_path"]]
    return cfg

def parse_config(yaml_path):
    with open(yaml_path) as f:
        cfg = yaml.safe_load(f)

    validate_required(cfg)
    cfg = apply_defaults(cfg)
    cfg = normalize(cfg)

    return cfg
