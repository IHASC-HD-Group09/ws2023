from typing import List
from pathlib import Path
import shlex
import subprocess
import argparse
import shutil
from itertools import combinations
import json


CXX_O_FLAGS: List[str] = [
    "-O0",
    "-O1",
    "-O2",
    "-O3",
]


CXX_FLAGS: List[str] = [
    "-funroll-loops",
    "-ffast-math",
    "-march=native",
]


def run_with_flags(flags: str, build_cwd: str) -> str:
    subprocess.run(
        shlex.split(f"cmake -DCMAKE_CXX_FLAGS='{flags}' .."),
        cwd=build_cwd,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        shlex.split("make"),
        cwd=build_cwd,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    result: str = subprocess.run(
        shlex.split("./main"),
        cwd=build_cwd,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return result

def save_result(result: str, flags: str) -> None:
    results: dict = {}
    try:
        with open(f"./results/results.json", "r") as f:
            results = json.load(f)
    except FileNotFoundError:
        pass
    results[flags] = result.splitlines()
    with open(f"./results/results.json", "w") as f:
        json.dump(results, f, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Builds and runs the CMake project in this directory.\n"
            "Build files will be placed in `./build` and be removed after running.\n"
            "Results will be placed in `./results`.\n"
        )
    )
    parser.add_argument(
        "--flags",
        type=str,
        nargs="+",
        help=f"Compiler flags, if none provided will run with every combination of {CXX_FLAGS}",
    )
    args = parser.parse_args()

    flags: List[str] = args.flags if args.flags else CXX_FLAGS
    do_combinations: bool = args.flags is None

    build_dir: Path = Path("./build")
    shutil.rmtree(build_dir, ignore_errors=True)
    build_dir.mkdir(exist_ok=True)

    results_dir: Path = Path("./results")
    shutil.rmtree(results_dir, ignore_errors=True)
    results_dir.mkdir(exist_ok=True)


    for o_flag in CXX_O_FLAGS:
        result: str = run_with_flags(o_flag, str(build_dir))
        save_result(result, o_flag)
        # run with all combinations of flags
        for i in range(1, len(flags) + 1):
            for combo in combinations(flags, i):
                combo_flags: str = f"{o_flag} {' '.join(combo)}"
                result: str = run_with_flags(combo_flags, str(build_dir))
                save_result(result, combo_flags)