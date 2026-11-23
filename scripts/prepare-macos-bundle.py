"""Relocate the bundled R framework and seal the Apple Silicon application."""
import os
from pathlib import Path
import subprocess
import sys


def run(*args):
    subprocess.run(args, check=True)


def relocate_r(atop):
    home = (atop / "runtime/R.framework/Resources").resolve(strict=True)
    lib = home / "lib"
    for link in (home / "fontconfig/fonts/conf.d").iterdir():
        if link.is_symlink() and link.readlink().is_absolute():
            target = home / "fontconfig/fontconfig/conf.avail" / link.name
            target.resolve(strict=True)
            link.unlink()
            link.symlink_to(os.path.relpath(target, link.parent))

    # Rscript starts a second R process for SVG export and parallel workers.
    launcher = home / "bin/R"
    text = launcher.read_text()
    prefix = "/Library/Frameworks/R.framework/Resources"
    if f'R_HOME_DIR="{prefix}"' in text:
        text = text.replace(f'R_HOME_DIR="{prefix}"',
                            'R_HOME_DIR=$(CDPATH= cd -P "$(dirname "$0")/.." && pwd)')
        launcher.write_text(text.replace(prefix + "/", "${R_HOME}/"))
    rscript = home / "bin/Rscript"
    binary = home / "bin/exec/Rscript"
    if not binary.exists():
        rscript.rename(binary)
        rscript.write_text('#!/bin/sh\nset -eu\nRHOME=$(CDPATH= cd -P "$(dirname "$0")/.." && pwd)\n'
                           'export RHOME\nexec "$RHOME/bin/exec/Rscript" "$@"\n')
        rscript.chmod(0o755)

    candidates = [atop / "analysis/lib/atop_permutation_engine.so", home / "bin/exec/R"]
    candidates.extend(home.rglob("*.so"))
    candidates.extend(home.rglob("*.dylib"))
    prefix = "/Library/Frameworks/R.framework/Versions/4.6/Resources/lib/"
    for file in candidates:
        if file.is_symlink():
            continue
        output = subprocess.check_output(["otool", "-arch", "arm64", "-L", file], text=True)
        changes = []
        seen = set()
        for line in output.splitlines()[1:]:
            name = line.strip().split(" (compatibility version")[0]
            if not name.startswith(prefix) or name in seen:
                continue
            seen.add(name)
            target = lib / Path(name).name
            if target.resolve(strict=True) == file.resolve() and file.suffix == ".dylib":
                changes.extend(["-id", "@loader_path/" + file.name])
            else:
                changes.extend(["-change", name, "@loader_path/" + os.path.relpath(target, file.parent)])
        if changes:
            run("install_name_tool", *changes, file)
            run("codesign", "--force", "--sign", "-", file)


def seal(app):
    frameworks = app / "Contents/Frameworks"
    # Sign nested code before sealing its enclosing bundle.
    paths = sorted((p for p in frameworks.rglob("*") if not p.is_symlink()
                    and (p.suffix in {".app", ".framework", ".dylib", ".xpc"}
                         or p.name in {"chrome_crashpad_handler", "ShipIt"})),
                   key=lambda p: len(p.parts), reverse=True)
    for path in paths:
        check = subprocess.run(["codesign", "--verify", "--strict", path], capture_output=True)
        if check.returncode:
            run("codesign", "--force", "--sign", "-", path)
    run("codesign", "--force", "--sign", "-", app)
    run("codesign", "--verify", "--deep", "--strict", app)


if __name__ == "__main__":
    app = Path(sys.argv[1]).resolve(strict=True)
    relocate_r(app / "Contents/Resources/atop")
    seal(app)
