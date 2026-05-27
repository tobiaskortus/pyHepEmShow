import os
import sys
import subprocess
import glob
import shutil
import sysconfig
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _BuildExt

VARIANT = "msc_no_ad_fix_yes"
ROOT = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.join(ROOT, "hepemshow", f"build_{VARIANT}")
G4HEPEM_INSTALL = os.path.join(ROOT, "externals", "g4hepem", f"install_{VARIANT}")


def _so_files():
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
    return glob.glob(os.path.join(BUILD_DIR, f"pyhepemshow*{ext_suffix}"))


def _g4hepem_libs():
    """Shared libraries that must be bundled alongside the extension."""
    lib_dir = os.path.join(G4HEPEM_INSTALL, "lib")
    if not os.path.isdir(lib_dir):
        lib_dir = os.path.join(G4HEPEM_INSTALL, "lib64")
    return (
        glob.glob(os.path.join(lib_dir, "libg4HepEm*.dylib"))
        + glob.glob(os.path.join(lib_dir, "libg4HepEm*.so*"))
    )


def _init_submodules():
    git_dir = os.path.join(ROOT, ".git")
    if not os.path.exists(git_dir):
        return
    submodule_dirs = [
        os.path.join(ROOT, "externals", "g4hepem"),
        os.path.join(ROOT, "externals", "CoDiPack"),
    ]
    needs_init = any(
        not os.path.isdir(d) or not os.listdir(d) for d in submodule_dirs
    )
    if needs_init:
        subprocess.run(
            ["git", "submodule", "update", "--init", "--recursive"],
            check=True, cwd=ROOT,
        )


def _ensure_built():
    _init_submodules()
    # First time: configure (cmake) + build
    if not os.path.isdir(BUILD_DIR):
        script = os.path.join(ROOT, "externals", "configure.sh")
        subprocess.run(
            ["bash", script, VARIANT, "", sys.executable],
            check=True, cwd=ROOT,
        )
    # Subsequent times: just build (make)
    else:
        subprocess.run(["make", "-j", "pyhepemshow"], check=True, cwd=BUILD_DIR)
    return _so_files()


class CMakeBuildExt(_BuildExt):
    def run(self):
        so_files = _ensure_built()
        libs = _g4hepem_libs()

        for f in so_files + libs:
            shutil.copy2(f, ROOT)

        os.makedirs(self.build_lib, exist_ok=True)
        for f in so_files + libs:
            shutil.copy2(f, self.build_lib)


setup(
    cmdclass={"build_ext": CMakeBuildExt},
    ext_modules=[Extension("pyhepemshow", sources=[])],
    packages=["hepemshow"],
    package_data={"hepemshow": ["data/hepem_data.json"]},
    include_package_data=True,
)
