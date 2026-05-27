import os
import sys
import subprocess
import glob
import shutil
import sysconfig
import tomllib

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _BuildExt
from setuptools.command.build_py import build_py as _BuildPy

VARIANT = "msc_no_ad_fix_yes"
ROOT = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.path.join(ROOT, "hepemshow", f"build_{VARIANT}")
G4HEPEM_INSTALL = os.path.join(ROOT, "externals", "g4hepem", f"install_{VARIANT}")


with open(os.path.join(ROOT, "pyproject.toml"), "rb") as _f:
    BUILD_TYPE = tomllib.load(_f)["tool"]["pyhepemshow"]["build_type"]


def _so_files():
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
    return glob.glob(os.path.join(BUILD_DIR, f"hepemshow*{ext_suffix}"))


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
            ["bash", script, VARIANT, "", sys.executable, BUILD_TYPE],
            check=True, cwd=ROOT,
        )
    # Subsequent times: update build type then make
    else:
        subprocess.run(
            ["cmake", ".", f"-DCMAKE_BUILD_TYPE={BUILD_TYPE}"],
            check=True, cwd=BUILD_DIR,
        )
        subprocess.run(["make", "-j", "hepemshow"], check=True, cwd=BUILD_DIR)
    return _so_files()


class CMakeBuildPy(_BuildPy):
    """Stage the JSON data file into the source tree before build_py collects package_data."""
    def run(self):
        src_json = os.path.join(ROOT, "hepemshow", "data", "hepem_data.json")
        dst_dir = os.path.join(ROOT, "pyhepemshow", "data")
        os.makedirs(dst_dir, exist_ok=True)
        shutil.copy2(src_json, os.path.join(dst_dir, "hepem_data.json"))
        super().run()


class CMakeBuildExt(_BuildExt):
    def run(self):
        so_files = _ensure_built()
        libs = _g4hepem_libs()

        # Copy the extension and shared libs into the pyhepemshow package directory
        pkg_src_dir = os.path.join(ROOT, "pyhepemshow")
        os.makedirs(pkg_src_dir, exist_ok=True)
        for f in so_files + libs:
            shutil.copy2(f, pkg_src_dir)

        build_lib_pkg = os.path.join(self.build_lib, "pyhepemshow")
        os.makedirs(build_lib_pkg, exist_ok=True)
        for f in so_files + libs:
            shutil.copy2(f, build_lib_pkg)


setup(
    cmdclass={"build_py": CMakeBuildPy, "build_ext": CMakeBuildExt},
    ext_modules=[Extension("pyhepemshow.hepemshow", sources=[])],
    packages=["pyhepemshow"],
    package_data={"pyhepemshow": ["data/hepem_data.json", "*.dylib", "*.so", "*.so.*"]},
    include_package_data=True,
)
