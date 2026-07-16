"""Retry compatibility helpers."""

from __future__ import annotations

import functools
import time
from collections.abc import Callable
from typing import ParamSpec, TypeVar

P = ParamSpec("P")
R = TypeVar("R")

try:  # pragma: no cover - exercised when tenacity is installed
    from tenacity import retry, stop_after_attempt, wait_exponential
except ModuleNotFoundError:  # pragma: no cover - simple fallback for minimal environments

    def stop_after_attempt(attempts: int) -> int:
        return attempts

    def wait_exponential(multiplier: float = 1, min: float = 1, max: float = 8) -> tuple[float, float, float]:
        return (multiplier, min, max)

    def retry(stop: int = 3, wait: tuple[float, float, float] = (1, 1, 8)) -> Callable[[Callable[P, R]], Callable[P, R]]:
        def decorator(func: Callable[P, R]) -> Callable[P, R]:
            @functools.wraps(func)
            def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
                last_exc: Exception | None = None
                multiplier, min_delay, max_delay = wait
                for attempt in range(1, stop + 1):
                    try:
                        return func(*args, **kwargs)
                    except Exception as exc:  # noqa: BLE001
                        last_exc = exc
                        if attempt == stop:
                            break
                        delay = max(min_delay, min(max_delay, multiplier * (2 ** (attempt - 1))))
                        time.sleep(delay)
                assert last_exc is not None
                raise last_exc

            return wrapper

        return decorator
