from dataclasses import dataclass
@dataclass
class ModelParams:
    k_I: float = 0.05
    k_dXI: float = 1.0
    k_A: float = 0.008
    k_B: float = 0.0001
    k_XA: float = 0.0
    k_XB: float = 0.03
    k_da: float = 0.01
    k_E: float = 0.006
    k_Xmax: float = 0.4
    k_Xhalf: float = 1.0
    k_Xdeg: float = 0.001
    k_fX: float = 0.2
    k_d: float = 1.0

SCALED_KEYS = {
    "k_I", "k_A", "k_B", "k_XA", "k_XB", "k_da", "k_E", "k_Xmax", "k_Xhalf", "k_Xdeg", "k_fX"
}
