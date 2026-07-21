"""
HocusPockets: Holistic Computacional Utility for Scoring Pockets

Authors: Machado, L.F.; Torres, P.H.M.; Spelta, G.I.; Oliveira Filho, J.A., Pascutti, P.G.

[Federal University of Rio de Janeiro]

Contact info:
Carlos Chagas Filho Institute of Biophysics
LMDM - Laboratory for Molecular Modeling and Dynamics
Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
Health Sciences Centre - block G1-19
E-mail address: luizafuruno@ufrj.br
"""

import logging
import os
import sys
from typing import Optional


def setup_logging(
    output_path: Optional[str] = None, verbose: bool = False
) -> logging.Logger:
    logger = logging.getLogger("hocuspockets")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)

    if verbose and output_path:
        fh = logging.FileHandler(os.path.join(output_path, "hocuspockets.log"))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(fh)

    return logger
