#!/usr/bin/env python3
"""
This wrapper is intended to be invoked by cargo as a runner.

It's invoked like ./codecov_wrapper.py <program> <arguments>

It will ensure that $LLVM_PROFILE_FILE is properly set, based on the hash of the executable.
"""
import os
import sys
from pathlib import Path

# Add the root of the swanky repo to the path.
sys.path.append(str(Path(__file__).parent.parent.parent.parent))

from etc.ci.xattr_hash import cached_blake2b

args = sys.argv[1:]
exe = Path(args[0]).resolve()
dst = Path(os.environ["SWANKY_CODECOV_DST"]) / cached_blake2b(exe).hex()
dst.mkdir(exist_ok=True, parents=True)
try:
    (dst / "exe").symlink_to(exe)
except FileExistsError:
    pass
os.environ["LLVM_PROFILE_FILE"] = str(dst / "coverage-%p-%m.profraw")
if "SWANKY_CODECOV_NEXT_RUNNER" in os.environ:
    args.insert(0, os.environ["SWANKY_CODECOV_NEXT_RUNNER"])
os.execvp(args[0], args)
