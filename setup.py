"""
Custom build script for PASTAQ package.
Handles CMake-based C++ extension building.
"""
import os
import sys
import re
import platform
import subprocess
import multiprocessing
from pathlib import Path

from packaging.version import Version
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name: str, sourcedir: str = "."):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def initialize_options(self):
        super().initialize_options()
        self.debug = getattr(self, "debug", False)

    def finalize_options(self):
        super().finalize_options()
        env_debug = os.environ.get("PASTAQ_DEBUG")
        if env_debug is not None:
            self.debug = self.debug or (env_debug == "1")

    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"], text=True)
        except OSError as e:
            raise RuntimeError("CMake is required to build PASTAQ C++ extension") from e

        m = re.search(r"version\s*([\d.]+)", out)
        if not m or Version(m.group(1)) < Version("3.14.0"):
            raise RuntimeError("CMake >= 3.14 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext: CMakeExtension):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cfg = "Debug" if self.debug else "Release"

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DPASTAQ_ENABLE_TESTS=OFF",
            "-DEIGEN_BUILD_DOC=OFF",
            "-DEIGEN_BUILD_TESTING=OFF",
            "-DDOWNLOAD_CATCH=OFF",
        ]

        generator = os.environ.get("CMAKE_GENERATOR", "")
        # Handle multi-config generators
        if platform.system() != "Windows" and any(g in generator for g in ("Xcode", "Ninja Multi-Config")):
            cmake_args.append(f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}")

        extra = os.environ.get("PASTAQ_CMAKE_ARGS", "")
        if extra:
            cmake_args.extend(extra.split())

        build_args = ["--config", cfg]
        if platform.system() == "Windows":
            cmake_args.append(f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}")
            if sys.maxsize > 2**32:
                cmake_args.append("-A")
                cmake_args.append("x64")
            build_args.extend(["--", "/m"])
        else:
            jobs = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL") or str(max(1, multiprocessing.cpu_count() - 1))
            build_args.extend(["--", f"-j{jobs}"])

        # Portable VERSION_INFO define
        env = os.environ.copy()
        version = self.distribution.get_version()
        is_msvc = platform.system() == "Windows" and ("MSC" in sys.version or generator.startswith("Visual Studio"))
        define_flag = "/D" if is_msvc else "-D"
        for var in ("CXXFLAGS", "CFLAGS"):
            env[var] = (env.get(var, "") + f" {define_flag}VERSION_INFO=\\\"{version}\\\"").strip()

        build_temp = Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)

        print(f"[PASTAQ] Configuring CMake in {build_temp}")
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=str(build_temp), env=env)

        print(f"[PASTAQ] Building extension to {extdir}")
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=str(build_temp))


# Configure extension: dotted path ensures correct placement under the pastaq package
setup(
    ext_modules=[CMakeExtension("pastaq.pastaq_cpp", sourcedir=".")],
    cmdclass={"build_ext": CMakeBuild},
)
