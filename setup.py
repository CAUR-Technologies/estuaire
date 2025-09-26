#!/usr/bin/env python3
"""Setuptools build script for the estuaire fast-marching extensions.

Notes:
- Under PEP 517 build isolation, declare build-time deps in pyproject.toml:
  numpy, Cython, and Mako (since templates render at build time).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from distutils import sysconfig


HERE = Path(__file__).resolve().parent
BUILD_DIR = HERE / "build"
INCLUDE_DIR = BUILD_DIR / "include"
PYX_DIR = BUILD_DIR / "estuaire"
TEMPLATES_DIR = HERE / "templates"
PXD_SRC_DIR = HERE / "src" / "estuaire"


def render_template(src: Path, dest: Path) -> None:
    # Lazy import so build isolation (or no-isolation) can satisfy dependency
    from mako.template import Template

    dest.parent.mkdir(parents=True, exist_ok=True)
    template = Template(filename=str(src))
    dest.write_text(template.render())


def ensure_generated_sources() -> None:
    # Header templates
    for name in ("nvect", "narray", "solver", "raytrace"):
        render_template(
            TEMPLATES_DIR / f"{name}.hpp.mako",
            INCLUDE_DIR / f"{name}.hpp",
        )

    # Cython sources
    PYX_DIR.mkdir(parents=True, exist_ok=True)
    (PYX_DIR / "__init__.py").touch()
    for pxd in PXD_SRC_DIR.glob("*.pxd"):
        shutil.copy(pxd, PYX_DIR / pxd.name)

    render_template(
        TEMPLATES_DIR / "solver.pyx.mako",
        BUILD_DIR / "solver.pyx",
    )
    render_template(
        TEMPLATES_DIR / "raytrace.pyx.mako",
        BUILD_DIR / "raytrace.pyx",
    )


COMMON_ARGS = {
    "include_dirs": [str(INCLUDE_DIR)],
    "language": "c++",
    "extra_compile_args": ["-std=c++14", "-Wno-unused-function"],
    # "extra_link_args": ["-stdlib=libc++"],
    "extra_link_args": [],
}


extensions = [
    Extension(
        "estuaire.solver",
        sources=[str(BUILD_DIR / "solver.pyx")],
        **COMMON_ARGS,
    ),
    Extension(
        "estuaire.raytrace",
        sources=[str(BUILD_DIR / "raytrace.pyx")],
        **COMMON_ARGS,
    ),
]


class BuildExt(build_ext):
    """Custom build_ext that renders templates and patches macOS flags."""

    def run(self) -> None:
        # Resolve a valid SDK path (or none) and set/unset SDKROOT accordingly
        os.environ.pop("SDKROOT", None)
        self.sdk_path = ""
        if shutil.which("xcrun"):
            try:
                p = subprocess.check_output(
                    ["xcrun", "--show-sdk-path"], text=True
                ).strip()
                if p and Path(p).exists():
                    self.sdk_path = p
                    os.environ["SDKROOT"] = self.sdk_path
            except subprocess.CalledProcessError:
                self.sdk_path = ""

        # Scrub -isysroot in sysconfig and compiler executables, then render sources
        self._patch_compiler_config()
        ensure_generated_sources()

        # If Cython is available, cythonize the .pyx sources and set include_path
        try:
            from Cython.Build import cythonize  # type: ignore
        except Exception:
            cythonize = None
        if cythonize is not None:
            self.extensions = cythonize(
                self.extensions,
                language_level=3,
                include_path=[str(PYX_DIR)],
            )

        # Ensure setuptools-specific bits exist, and enforce final -isysroot
        for ext in self.extensions:
            if not hasattr(ext, "_needs_stub"):
                ext._needs_stub = False
            if self.sdk_path:
                ext.extra_compile_args = list(getattr(
                    ext, "extra_compile_args", []
                )) + ["-isysroot", self.sdk_path]
                ext.extra_link_args = list(getattr(
                    ext, "extra_link_args", []
                )) + ["-isysroot", self.sdk_path]

        super().run()

    def finalize_options(self) -> None:
        super().finalize_options()
        # Lazy import so build isolation can install numpy first
        import numpy as np  # noqa: WPS433
        np_inc = np.get_include()
        for ext in self.extensions:
            if np_inc not in ext.include_dirs:
                ext.include_dirs.append(np_inc)

    def build_extensions(self) -> None:
        # Scrub any lingering -isysroot in compiler/linker executables
        self._scrub_compiler_executables()
        super().build_extensions()

    def _scrub_compiler_executables(self) -> None:
        comp = self.compiler
        # Distutils compilers keep command lists in these attributes
        attrs = [
            "compiler",
            "compiler_so",
            "compiler_cxx",
            "linker_so",
            "linker_exe",
            "archiver",
        ]
        for attr in attrs:
            cmd = getattr(comp, attr, None)
            if not cmd:
                continue
            parts = cmd.split() if isinstance(cmd, str) else list(cmd)
            i = 0
            while i < len(parts):
                if parts[i] == "-isysroot":
                    # Remove -isysroot <path> unless we have a valid SDK
                    if not self.sdk_path:
                        del parts[i : i + 2]
                        continue
                    # Replace path if it differs
                    if i + 1 < len(parts):
                        parts[i + 1] = self.sdk_path
                    else:
                        parts.append(self.sdk_path)
                    i += 2
                    continue
                i += 1
            setattr(comp, attr, " ".join(parts) if isinstance(cmd, str) else parts)

    def _patch_compiler_config(self) -> None:
        vars = sysconfig.get_config_vars()
        keys = ["CC", "CXX", "LDSHARED", "CFLAGS", "LDFLAGS", "CXXFLAGS"]
        for key in keys:
            val = vars.get(key)
            if not isinstance(val, str):
                continue
            parts = val.split()
            i = 0
            while i < len(parts):
                if parts[i] == "-isysroot":
                    if not self.sdk_path:
                        del parts[i : i + 2]
                        continue
                    if i + 1 < len(parts):
                        parts[i + 1] = self.sdk_path
                    else:
                        parts.append(self.sdk_path)
                    i += 2
                    continue
                i += 1
            vars[key] = " ".join(parts)


        super().run()

    def finalize_options(self) -> None:
        super().finalize_options()
        # Lazy import so build isolation can install numpy first
        import numpy as np  # noqa: WPS433
        np_inc = np.get_include()
        for ext in self.extensions:
            if np_inc not in ext.include_dirs:
                ext.include_dirs.append(np_inc)

    def build_extensions(self) -> None:
        compiler = self.compiler
        for attr in ("compiler", "compiler_so", "linker_so"):
            cmd = getattr(compiler, attr, None)
            if cmd is None:
                continue
            parts = cmd.split() if isinstance(cmd, str) else list(cmd)

            i = 0
            while i < len(parts):
                if parts[i] == "-isysroot":
                    if self.sdk_path:
                        if i + 1 < len(parts):
                            parts[i + 1] = self.sdk_path
                        else:
                            parts.append(self.sdk_path)
                        i += 2
                        continue
                    else:
                        del parts[i : i + 2]
                        continue
                i += 1

            if isinstance(cmd, str):
                setattr(compiler, attr, " ".join(parts))
            else:
                setattr(compiler, attr, parts)

        super().build_extensions()

    def _patch_compiler_config(self) -> None:
        vars = sysconfig.get_config_vars()
        keys = ["CC", "CXX", "LDSHARED", "CFLAGS", "LDFLAGS", "CXXFLAGS"]
        for key in keys:
            val = vars.get(key)
            if not isinstance(val, str):
                continue
            parts = val.split()
            i = 0
            while i < len(parts):
                if parts[i] == "-isysroot":
                    if self.sdk_path:
                        if i + 1 < len(parts):
                            parts[i + 1] = self.sdk_path
                        else:
                            parts.append(self.sdk_path)
                        i += 2
                        continue
                    else:
                        del parts[i : i + 2]
                        continue
                i += 1
            vars[key] = " ".join(parts)


setup(
    name="estuaire",
    version="1.0.0",
    description="Fast marching estuaire solver with ray tracing",
    author="J.-P. Mercier",
    packages=["estuaire", "estuaire.core"],
    include_package_data=True,
    ext_modules=extensions,
    cmdclass={"build_ext": BuildExt},
)
