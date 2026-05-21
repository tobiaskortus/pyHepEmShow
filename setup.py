import os
import subprocess
from setuptools import setup
from setuptools.command.build import build as _Build

VARIANT = "msc_no_ad_fix_yes"
ROOT = os.path.dirname(os.path.abspath(__file__))
BINARY_REL = os.path.join("hepemshow", f"build_{VARIANT}", "HepEmShow")
BINARY_ABS = os.path.join(ROOT, BINARY_REL)


class Build(_Build):
    def run(self):
        if not os.path.isfile(BINARY_ABS):
            script = os.path.join(ROOT, "externals", "configure.sh")
            subprocess.run(["bash", script, VARIANT, ""], check=True, cwd=ROOT)
        super().run()


setup(
    cmdclass={"build": Build},
    data_files=[("bin", [BINARY_REL])],
)
