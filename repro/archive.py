#!/usr/bin/env python3
"""Archive-only publisher: durable verification -> no-replace rename -> read-only.

This module never reads a scientific object, transforms values, or runs a model.
The caller must close all writers before publication. Published paths are never
moved, renamed, overwritten, or made writable by this module.
"""

from __future__ import annotations

import ctypes
import errno
import hashlib
import os
import stat
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def regular_files(directory: Path) -> list[Path]:
    result = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise RuntimeError(f"ARCHIVE_SYMLINK_FORBIDDEN:{path}")
        if path.is_file():
            result.append(path)
        elif not path.is_dir():
            raise RuntimeError(f"ARCHIVE_SPECIAL_FILE_FORBIDDEN:{path}")
    return result


def read_hashes(directory: Path, name: str = "hashes.sha256") -> dict[str, str]:
    hashes = {}
    for line in (directory / name).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            raise RuntimeError("ARCHIVE_BLANK_HASH_LINE")
        digest, relative = line.split("  ", 1)
        path = Path(relative)
        if (len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest)
                or path.is_absolute() or ".." in path.parts or relative in hashes
                or relative == name):
            raise RuntimeError("ARCHIVE_INVALID_OR_DUPLICATE_HASH_ENTRY")
        hashes[relative] = digest
    if not hashes:
        raise RuntimeError("ARCHIVE_EMPTY_HASH_MANIFEST")
    return hashes


def verify(directory: Path, required: set[str], closed: bool = True) -> dict[str, str]:
    if directory.is_symlink() or not directory.is_dir():
        raise RuntimeError("ARCHIVE_INVALID_DIRECTORY")
    files = regular_files(directory)
    observed = {str(p.relative_to(directory)) for p in files}
    if not required <= observed:
        raise RuntimeError(f"ARCHIVE_REQUIRED_FILES_MISSING:{sorted(required - observed)}")
    hashes = read_hashes(directory)
    if closed and set(hashes) != observed - {"hashes.sha256"}:
        raise RuntimeError("ARCHIVE_HASH_INVENTORY_NOT_CLOSED")
    for relative, expected in hashes.items():
        path = directory / relative
        if not path.is_file() or sha256(path) != expected:
            raise RuntimeError(f"ARCHIVE_HASH_MISMATCH:{relative}")
    return hashes


def write_hashes(directory: Path) -> None:
    files = [p for p in regular_files(directory) if p.name != "hashes.sha256"]
    destination = directory / "hashes.sha256"
    # All writers are closed; this is the only writable, unpublished manifest.
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for path in files:
            handle.write(f"{sha256(path)}  {path.relative_to(directory).as_posix()}\n")
        handle.flush()
        os.fsync(handle.fileno())


def rename_noreplace(source: Path, destination: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    try:
        renameat2 = libc.renameat2
    except AttributeError as error:
        raise RuntimeError("ATOMIC_NOREPLACE_RENAME_UNAVAILABLE") from error
    renameat2.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
                          ctypes.c_char_p, ctypes.c_uint]
    renameat2.restype = ctypes.c_int
    result = renameat2(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
    if result:
        code = ctypes.get_errno()
        if code == errno.EEXIST:
            raise FileExistsError(code, os.strerror(code), str(destination))
        raise OSError(code, os.strerror(code), str(destination))


def publish(staging: Path, destination: Path, required: set[str]) -> dict[str, object]:
    staging = staging.absolute()
    destination = destination.absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"PUBLISHED_TASK_IMMUTABLE:{destination}")
    if not staging.is_dir() or not destination.parent.is_dir():
        raise RuntimeError("ARCHIVE_PATH_MISSING")
    if staging.stat().st_dev != destination.parent.stat().st_dev:
        raise RuntimeError("ARCHIVE_NOT_ON_SAME_FILESYSTEM")
    files = regular_files(staging)
    for path in files:
        with path.open("rb") as handle:
            os.fsync(handle.fileno())
    directories = [p for p in staging.rglob("*") if p.is_dir()]
    for directory in sorted(directories, key=lambda p: len(p.parts), reverse=True):
        fsync_directory(directory)
    fsync_directory(staging)
    expected = verify(staging, required)
    # This is the sole directory rename. No-replace is enforced in the kernel.
    rename_noreplace(staging, destination)
    fsync_directory(staging.parent)
    fsync_directory(destination.parent)
    # chmod is deliberately after publication. Do not move destination again.
    for path in regular_files(destination):
        os.chmod(path, 0o444)
    for directory in sorted([p for p in destination.rglob("*") if p.is_dir()],
                            key=lambda p: len(p.parts), reverse=True):
        os.chmod(directory, 0o555)
        fsync_directory(directory)
    os.chmod(destination, 0o555)
    fsync_directory(destination)
    fsync_directory(destination.parent)
    observed = verify(destination, required)
    if observed != expected or staging.exists():
        raise RuntimeError("ARCHIVE_POSTPUBLICATION_INTEGRITY_FAILURE")
    if any(stat.S_IMODE(p.stat().st_mode) & 0o222
           for p in [destination, *destination.rglob("*")]):
        raise RuntimeError("ARCHIVE_READONLY_FAILURE")
    return {"atomic_noreplace_rename": "PASS", "required_files_and_sha256": "PASS",
            "postrename_readonly": "PASS", "published_file_count": len(observed) + 1,
            "task_hash_manifest_sha256": sha256(destination / "hashes.sha256"),
            "final_directory": str(destination)}
