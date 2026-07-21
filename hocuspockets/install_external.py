# Description
###############################################################################
"""
This module provides functionality to install external software dependencies required by the HocusPockets pipeline.
"""

# Imports
###############################################################################
import logging
import subprocess
import urllib.request
from pathlib import Path
from hocuspockets.constants import MOPAC_BINARY

logger = logging.getLogger(__name__)

# License
###############################################################################
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


class InstallExternalSoftware:
    """
    Class to handle the installation of external software dependencies.
    """

    def __init__(self, path: str = Path.home()):
        self.path = path

    def install(self) -> None:
        """
        Install the required external software dependencies.
        """
        self.__install_fpocket()
        self.__install_mopac()
        self.__install_ConCavity()
        self.__install_primordia()

    def __install_fpocket(self) -> None:
        """
        Install fpocket software.
        """
        try:
            subprocess.run(
                [
                    "conda",
                    "run",
                    "-n",
                    "HocusPockets",
                    "conda",
                    "install",
                    "-y",
                    "fpocket",
                ],
                check=True,
            )
        except subprocess.CalledProcessError as e:
            logger.error("Failed to install fpocket: %s", e)

    def __install_mopac(self) -> None:
        """
        Install mopac software.
        """
        try:
            subprocess.run(
                [
                    "conda",
                    "run",
                    "-n",
                    "HocusPockets",
                    "conda",
                    "install",
                    "-y",
                    MOPAC_BINARY,
                ],
                check=True,
            )
        except subprocess.CalledProcessError as e:
            logger.error("Failed to install mopac: %s", e)

    def __install_ConCavity(self) -> None:
        """
        Install ConCavity software.
        """
        try:
            urllib.request.urlretrieve(
                "https://compbio.cs.princeton.edu/concavity/concavity_distr.tar.gz",
                f"{self.path}/concavity_distr.tar.gz",
            )
        except Exception as e:
            logger.error("Failed to download ConCavity: %s", e)
        try:
            subprocess.run(
                ["tar", "-xzf", f"{self.path}/concavity_distr.tar.gz", "-C", self.path],
                check=True,
            )
            subprocess.run(["make", "-C", f"{self.path}/concavity_distr"], check=True)
        except subprocess.CalledProcessError as e:
            logger.error("Failed to install ConCavity: %s", e)
        try:
            with open(f"{self.path}/.bashrc", "a") as f:
                f.write(
                    f"\nexport PATH={self.path}/concavity_distr/bin/x86_64/:$PATH\n"
                )
        except Exception as e:
            logger.error("Failed to update .bashrc for ConCavity: %s", e)

    def __install_primordia(self) -> None:
        """
        Install Primordia software.
        """
        try:
            urllib.request.urlretrieve(
                "https://github.com/igorChem/PRIMoRDiA1.0v/archive/master.zip",
                f"{self.path}/primordia.zip",
            )
        except Exception as e:
            logger.error("Failed to download PRIMoRDiA: %s", e)
        try:
            subprocess.run(
                ["unzip", f"{self.path}/primordia.zip", "-d", self.path], check=True
            )
        except subprocess.CalledProcessError as e:
            logger.error("Failed to install PRIMoRDiA: %s", e)

        try:
            with open(f"{self.path}/.bashrc", "a") as f:
                f.write(f"\nexport PATH={self.path}/PRIMoRDiA1.0v-master/:$PATH\n")
        except Exception as e:
            logger.error("Failed to update .bashrc for PRIMoRDiA: %s", e)


def install_external_software() -> None:
    """
    Install the required external software dependencies for HocusPockets.
    """
    try:
        installer = InstallExternalSoftware()
        installer.install()
    except Exception as e:
        logger.error("Installation of external software failed: %s", e)
