#!/usr/bin/env python3
"""Build a native solver and freeze the GUI with PyInstaller.

Run this script *on each target operating system*. PyInstaller is not a
cross-compiler, so Windows packages are built on Windows, macOS packages on
macOS, and Linux packages on Linux.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from build_solver import build, find_solver

ROOT = Path(__file__).resolve().parent


def run(cmd: list[str]) -> None:
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("onedir", "onefile"), default="onedir",
                    help="onedir is recommended for first builds and scientific Python apps")
    ap.add_argument("--boost-root", default=os.environ.get("BOOST_ROOT", ""))
    ap.add_argument("--clean", action="store_true", help="remove previous CMake/PyInstaller outputs")
    ap.add_argument(
        "--skip-solver-build",
        action="store_true",
        help="reuse an already-built solver from build/ instead of invoking CMake again",
    )
    ap.add_argument("--name", default="AnalyticMagneticField")
    args = ap.parse_args()

    pyinstaller = shutil.which("pyinstaller")
    if pyinstaller is None:
        raise SystemExit("PyInstaller is not installed. Run: pip install -r requirements-packaging.txt")

    build_dir = ROOT / "build"
    if args.clean:
        for d in (ROOT / "dist", ROOT / "pyinstaller-build"):
            if d.exists():
                shutil.rmtree(d)

    if args.skip_solver_build:
        solver = find_solver(build_dir)
        if solver is None:
            raise SystemExit(
                "--skip-solver-build was requested, but no field_solver executable "
                "was found under build/. Run build_solver.py first."
            )
        print(f"Using existing solver: {solver}")
    else:
        solver = build(build_dir, args.boost_root, clean=args.clean)

    cmd = [
        pyinstaller,
        "--noconfirm",
        "--clean",
        "--windowed",
        "--name", args.name,
        "--distpath", str(ROOT / "dist"),
        "--workpath", str(ROOT / "pyinstaller-build"),
        "--specpath", str(ROOT / "pyinstaller-build"),
        "--add-binary", f"{solver}:bin",
    ]
    if args.mode == "onefile":
        cmd.append("--onefile")
    else:
        cmd.append("--onedir")

    # Gives macOS bundles a stable bundle identifier. This does not sign or
    # notarize the application; those are distribution/publishing steps.
    if sys.platform == "darwin":
        cmd.extend(["--osx-bundle-identifier", "org.analyticmagneticfield.app"])

    cmd.append(str(ROOT / "app.py"))
    run(cmd)

    print("\nPackage created under:", ROOT / "dist")
    if sys.platform == "darwin" and args.mode == "onedir":
        print(f"macOS app bundle: {ROOT / 'dist' / (args.name + '.app')}")
    elif os.name == "nt":
        if args.mode == "onefile":
            print(f"Windows executable: {ROOT / 'dist' / (args.name + '.exe')}")
        else:
            print(f"Windows folder bundle: {ROOT / 'dist' / args.name}")
    else:
        print(f"Linux bundle: {ROOT / 'dist' / args.name}")


if __name__ == "__main__":
    main()
