def fmt_float(x: float) -> str:
    """
    Convert float to filesystem-safe string.
    Example:
      0.0  -> "0"
      2.5  -> "2p5"
      1e-3 -> "0p001"
    """
    return f"{x:g}".replace(".", "p")

