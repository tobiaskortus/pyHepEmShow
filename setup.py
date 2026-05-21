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


def _so_files():
    ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
    return glob.glob(os.path.join(BUILD_DIR, f"pyhepemshow*{ext_suffix}"))


def _ensure_built():
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

        for so in so_files:
            shutil.copy2(so, ROOT)

        os.makedirs(self.build_lib, exist_ok=True)
        for so in so_files:
            shutil.copy2(so, self.build_lib)


setup(
    cmdclass={"build_ext": CMakeBuildExt},
    ext_modules=[Extension("pyhepemshow", sources=[])],
)
