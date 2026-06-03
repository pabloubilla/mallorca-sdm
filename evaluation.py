import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from typing import List, Dict

def per_species_auc(y_true: pd.DataFrame, y_score: np.ndarray, species: List[str]) -> Dict[str, float]:
    scores: Dict[str, float] = {}
    for i, sp in enumerate(species):
        try:
            auc = roc_auc_score(y_true[sp].values, y_score[:, i])
        except ValueError:
            auc = np.nan
        scores[sp] = auc
    return scores