"""
Console reporting.

All user facing output of the pipeline goes through this module so that the
driver process and the per-chromosome worker processes share one visual style.
Only ASCII characters are emitted: the log is routinely redirected to a file on
clusters where the locale is not UTF-8.
"""

import sys
import time

LINE_WIDTH = 80


def _emit(text):
    """Write a single line to stdout, flushing so piped logs stay ordered."""
    sys.stdout.write(text + "\n")
    sys.stdout.flush()


# --------------------------------------------------------------------------- #
# Structural elements
# --------------------------------------------------------------------------- #
def banner(title):
    """Print a full-width section banner, e.g. the header of a pipeline stage."""
    _emit("")
    _emit("=" * LINE_WIDTH)
    _emit(title.upper().center(LINE_WIDTH))
    _emit("=" * LINE_WIDTH)


def rule(label=None, char="-"):
    """Print a horizontal rule, optionally carrying a left aligned label."""
    if label is None:
        _emit(char * LINE_WIDTH)
        return
    prefix = char * 3 + " " + label + " "
    _emit(prefix + char * max(LINE_WIDTH - len(prefix), 3))


# --------------------------------------------------------------------------- #
# Message levels
# --------------------------------------------------------------------------- #
def info(message):
    _emit("* " + message)


def detail(message):
    _emit("    " + message)


def ok(message):
    _emit("[OK] " + message)


def warn(message):
    _emit("[WARNING] " + message)


def step(chrom, message):
    """Per-chromosome progress line emitted by the genotyping subprocess."""
    _emit("[{}] {}".format(chrom, message))


def substep(chrom, message):
    _emit("[{}]   {}".format(chrom, message))


# --------------------------------------------------------------------------- #
# Formatting helpers
# --------------------------------------------------------------------------- #
def count(value):
    """Group digits in thousands: 2745732 -> '2,745,732'."""
    return "{:,}".format(int(value))


def bases(length):
    """Render a genomic length with a sensible unit, e.g. '28.5 Mb'."""
    length = int(length)
    for scale, unit in ((1_000_000_000, "Gb"), (1_000_000, "Mb"), (1_000, "kb")):
        if length >= scale:
            return "{:.1f} {}".format(length / scale, unit)
    return "{} bp".format(length)


def duration(seconds):
    """Render an elapsed time as '1h 02m 13s', '8m 12s' or '42.1s'."""
    seconds = float(seconds)
    if seconds < 60:
        return "{:.1f}s".format(seconds)
    minutes, sec = divmod(int(round(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return "{}h {:02d}m {:02d}s".format(hours, minutes, sec)
    return "{}m {:02d}s".format(minutes, sec)


class Timer:
    """Context manager returning the wall-clock time spent inside the block."""

    def __enter__(self):
        self._start = time.perf_counter()
        return self

    def __exit__(self, *exc_info):
        self.elapsed = time.perf_counter() - self._start
        return False

    @property
    def seconds(self):
        return getattr(self, "elapsed", time.perf_counter() - self._start)
