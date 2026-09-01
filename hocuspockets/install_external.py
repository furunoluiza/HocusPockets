# Description
###############################################################################
"""
This module provides functionality to install external software dependencies required by the HocusPockets pipeline.
"""

# Imports
###############################################################################
import logging
import os
import subprocess
import urllib.error
import urllib.request
from pathlib import Path

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
                    "-c",
                    "conda-forge",
                    "mopac",
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
        except (OSError, urllib.error.URLError) as e:
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
        except OSError as e:
            logger.error("Failed to update .bashrc for ConCavity: %s", e)

    def __install_primordia(self) -> None:
        """
        Build and install PRIMoRDiA 1.50 from its CMake project.

        The executable is installed below ``self.path`` so this installer does
        not require root privileges. The install directory is added to both
        the current process and the user's bash startup file.
        """
        source_dir = Path(self.path) / "PRIMoRDiA_1.50v"
        build_dir = source_dir / "build"
        install_dir = source_dir / "install"
        executable = install_dir / "bin" / "PRIMoRDiA_1.50v"

        try:
            source_dir.parent.mkdir(parents=True, exist_ok=True)
            if not source_dir.exists():
                subprocess.run(
                    [
                        "git",
                        "clone",
                        "--depth",
                        "1",
                        "https://github.com/bardenChem/PRIMoRDiA_1.50v.git",
                        str(source_dir),
                    ],
                    check=True,
                )

            subprocess.run(
                [
                    "conda",
                    "run",
                    "-n",
                    "HocusPockets",
                    "conda",
                    "install",
                    "-y",
                    "-c",
                    "conda-forge",
                    "cmake",
                    "eigen=3.4",
                    "boost-cpp",
                    "zlib",
                    "make",
                ],
                check=True,
            )
            conda_prefix = subprocess.run(
                ["conda", "run", "-n", "HocusPockets", "printenv", "CONDA_PREFIX"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            subprocess.run(
                [
                    "conda",
                    "run",
                    "-n",
                    "HocusPockets",
                    "cmake",
                    "-S",
                    str(source_dir),
                    "-B",
                    str(build_dir),
                    "-DCMAKE_BUILD_TYPE=Release",
                    f"-DCMAKE_INSTALL_PREFIX={install_dir}",
                    f"-DCMAKE_INSTALL_RPATH={conda_prefix}/lib",
                ],
                check=True,
            )
            subprocess.run(
                [
                    "conda",
                    "run",
                    "-n",
                    "HocusPockets",
                    "cmake",
                    "--build",
                    str(build_dir),
                    "--parallel",
                ],
                check=True,
            )
            subprocess.run(
                [
                    "conda",
                    "run",
                    "-n",
                    "HocusPockets",
                    "cmake",
                    "--install",
                    str(build_dir),
                ],
                check=True,
            )
            if not executable.is_file():
                raise FileNotFoundError(f"PRIMoRDiA executable not found: {executable}")

            executable.chmod(executable.stat().st_mode | 0o111)
            bin_dir = str(executable.parent)
            os.environ["PATH"] = f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}"
            bashrc = Path(self.path) / ".bashrc"
            path_line = f"export PATH={bin_dir}:$PATH"
            existing = bashrc.read_text() if bashrc.exists() else ""
            if path_line not in existing:
                with bashrc.open("a") as f:
                    f.write(f"\n{path_line}\n")
            logger.info("PRIMoRDiA 1.50 installed at %s", executable)
        except (OSError, subprocess.CalledProcessError) as e:
            logger.error("Failed to install PRIMoRDiA 1.50: %s", e)


def install_external_software() -> None:
    """
    Install the required external software dependencies for HocusPockets.
    """
    try:
        installer = InstallExternalSoftware()
        installer.install()
    except (OSError, subprocess.CalledProcessError) as e:
        logger.error("Installation of external software failed: %s", e)
