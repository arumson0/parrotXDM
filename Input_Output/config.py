from typing import List, Optional, Dict

class Config:
    def __init__(self, yaml_dict: Dict):
        # --- Required top-level fields ---
        self.input_path: List[str] = yaml_dict["input_path"]
        self.file_type: str = yaml_dict["file_type"]
        self.functional: str = yaml_dict["functional"]
        self.damping: str = yaml_dict["damping"]
        self.output_unit: str = yaml_dict["output_unit"]

        # --- Overrides (now also at the top level) ---
        self.override_a1: Optional[float] = yaml_dict.get("override_a1")
        self.override_a2: Optional[float] = yaml_dict.get("override_a2")
        self.override_zdamp: Optional[float] = yaml_dict.get("override_zdamp")

        # --- Run settings (nested block) ---
        run_settings = yaml_dict.get("run_settings", {})
        self.verbose_conv: bool = run_settings.get("verbose_conv", True)
        self.pairwise: bool = run_settings.get("pairwise", True)
        self.triples: bool = run_settings.get("triples", True)
        self.boolc11: bool = run_settings.get("c11", True)
        self.extrapolate_triples: bool = run_settings.get("extrapolate_triples", True)
        self.r2tol: float = run_settings.get("r2tol", 0.99)

    def __repr__(self):
        return (
            f"Config(input_path={self.input_path}, file_type={self.file_type}, "
            f"functional={self.functional}, damping={self.damping}, "
            f"pairwise={self.pairwise}, triples={self.triples}, "
            f"extrapolate_triples={self.extrapolate_triples}, "
            f"overrides={{a1={self.override_a1}, a2={self.override_a2}, zdamp={self.override_zdamp}}})"
        )
