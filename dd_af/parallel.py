"""Generic embarrassingly-parallel map helper for independent per-task work
(sampling replicas, batch pocket detection, ...), copied from
dd_docking/parallel.py (itself copied from dd_overlay/parallel.py) to avoid
a cross-package runtime dependency; kept identical since the pattern is
domain-agnostic.
"""
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, Iterator, List, TypeVar

T = TypeVar("T")
R = TypeVar("R")


def parallel_map(
    func: Callable[[T], R], items: List[T], n_jobs: int = 1, announce: bool = True,
) -> Iterator[R]:
    """Yield func(item) for each item, in completion order.

    n_jobs=1 runs sequentially in-process (no subprocess overhead, default).
    n_jobs<=0 uses all available CPU cores (os.cpu_count()).
    n_jobs>1 uses that many worker processes.

    `func` and each item must be picklable (module-level function, plain
    dataclasses/numpy arrays qualify).
    """
    if not items:
        return
    if n_jobs == 1:
        for item in items:
            yield func(item)
        return

    workers = n_jobs if n_jobs and n_jobs > 0 else os.cpu_count()
    if announce:
        print(f"[parallel] using {workers} worker processes for {len(items)} tasks", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(func, item) for item in items]
        for future in as_completed(futures):
            yield future.result()
