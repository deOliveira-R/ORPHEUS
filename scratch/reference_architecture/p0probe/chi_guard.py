import numpy as np, dataclasses
from orpheus.derivations.common.xs_library import get_mixture, get_xs
B = get_mixture("B","2g")
print("B chi", B.chi, "SigP", B.SigP)
try:
    m = dataclasses.replace(B, chi=np.array([1.0,0.0]))
    print("replace accepted: chi", m.chi)
except Exception as e:
    print("refused:", type(e).__name__, str(e)[:200])
