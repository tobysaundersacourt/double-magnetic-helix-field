#!/usr/bin/env python3
"""Cross-platform CMake build helper for the native field solver."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def solver_name() -> str:
    return "field_solver.exe" if os.name == "nt" else "field_solver"


def find_solver(build_dir: Path) -> Path | None:
    name = solver_name()
    candidates = [
        build_dir / name,
        build_dir / "Release" / name,
        build_dir / "RelWithDebInfo" / name,
        build_dir / "Debug" / name,
    ]
    for p in candidates:
        if p.is_file():
            return p
    matches = list(build_dir.rglob(name)) if build_dir.exists() else []
    return matches[0] if matches else None


def run(cmd: list[str]) -> None:
    print("$", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def build(build_dir: Path, boost_root: str = "", clean: bool = False) -> Path:
    cmake = shutil.which("cmake")
    if cmake is None:
        raise SystemExit("CMake was not found on PATH.")

    if clean and build_dir.exists():
        shutil.rmtree(build_dir)

    configure = [
        cmake,
        "-S", str(ROOT),
        "-B", str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    if boost_root:
        configure.extend([f"-DBOOST_ROOT={boost_root}", f"-DBoost_ROOT={boost_root}"])

    run(configure)
    run([cmake, "--build", str(build_dir), "--config", "Release", "--parallel"])

    solver = find_solver(build_dir)
    if solver is None:
        raise SystemExit(f"Build succeeded but {solver_name()} was not found under {build_dir}")
    print(f"Built solver: {solver}")
    return solver


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default=str(ROOT / "build"))
    ap.add_argument("--boost-root", default=os.environ.get("BOOST_ROOT", ""))
    ap.add_argument("--clean", action="store_true")
    args = ap.parse_args()
    build(Path(args.build_dir).resolve(), args.boost_root, args.clean)


if __name__ == "__main__":
    main()
