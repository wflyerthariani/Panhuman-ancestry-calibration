"""Region reads from a remote bgzipped, tabix-indexed VCF using plain HTTP range requests.

htslib can read such files itself, but against the gnomAD bucket its HTTP/2 connection is
reset part-way through a region ("Error in the HTTP2 framing layer"), every time, for the
dense HGDP+1KG files. This reader does the same job with short, independent, retried range
requests: it looks up the region's file offsets in the local copy of the .tbi index, downloads
from the first block, and stops at the first record past the region.

Local files do not need this; the caller uses tabix on them directly.
"""

from __future__ import annotations

import gzip
import struct
import zlib
from pathlib import Path
from typing import Iterator

from .mavedb import _request

FIRST_PIECE = 512 << 10    # bytes in the first range request of a read; doubles up to MAX_PIECE
MAX_PIECE = 8 << 20
MAX_BLOCK = 1 << 16        # a BGZF block is at most 64 KiB


class Tbi:
    """The parts of a tabix index needed for region lookups."""

    def __init__(self, path: Path):
        data = gzip.decompress(path.read_bytes())
        if data[:4] != b"TBI\x01":
            raise ValueError(f"{path.name} is not a tabix index")
        n_ref, = struct.unpack_from("<i", data, 4)
        l_nm, = struct.unpack_from("<i", data, 32)
        names = data[36 : 36 + l_nm].split(b"\x00")[:n_ref]
        pos = 36 + l_nm
        self.bins: dict[str, dict[int, list[tuple[int, int]]]] = {}
        self.linear: dict[str, tuple[int, ...]] = {}
        for name in names:
            bins = {}
            n_bin, = struct.unpack_from("<i", data, pos)
            pos += 4
            for _ in range(n_bin):
                b, n_chunk = struct.unpack_from("<Ii", data, pos)
                pos += 8
                vals = struct.unpack_from(f"<{2 * n_chunk}Q", data, pos)
                pos += 16 * n_chunk
                bins[b] = list(zip(vals[::2], vals[1::2]))
            n_intv, = struct.unpack_from("<i", data, pos)
            pos += 4
            self.linear[name.decode()] = struct.unpack_from(f"<{n_intv}Q", data, pos)
            pos += 8 * n_intv
            self.bins[name.decode()] = bins

    def start_offset(self, contig: str, beg: int, end: int) -> int | None:
        """Virtual offset to start reading from for 0-based half-open [beg, end); None if nothing is there."""
        bins, linear = self.bins.get(contig), self.linear.get(contig, ())
        if not bins:
            return None
        min_off = linear[min(beg >> 14, len(linear) - 1)] if linear else 0
        wanted = [0]
        for shift, first in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
            wanted.extend(range(first + (beg >> shift), first + ((end - 1) >> shift) + 1))
        starts = [max(c_beg, min_off) for b in wanted for c_beg, c_end in bins.get(b, ()) if c_end > min_off]
        return min(starts) if starts else None


def _blocks(uri: str, coffset: int) -> Iterator[bytes]:
    """Decompressed BGZF blocks from compressed offset `coffset` onwards."""
    buf, at, piece = b"", coffset, FIRST_PIECE
    while True:
        r = _request("GET", uri, headers={"Range": f"bytes={at}-{at + piece - 1}"})
        if not r.content:
            return
        buf += r.content
        at += len(r.content)
        o = 0
        while o + 18 <= len(buf):
            if buf[o : o + 4] != b"\x1f\x8b\x08\x04":
                raise ValueError("not a BGZF block boundary")
            size = struct.unpack_from("<H", buf, o + 16)[0] + 1   # BSIZE; the BC subfield comes first in bgzip output
            if buf[o + 12 : o + 14] != b"BC":
                raise ValueError("unexpected BGZF extra field")
            if o + size > len(buf):
                break
            yield zlib.decompress(buf[o : o + size], 31)
            o += size
        buf = buf[o:]
        if len(r.content) < piece:
            return
        piece = min(2 * piece, MAX_PIECE)


def header(uri: str) -> str:
    """The VCF header, through the #CHROM line."""
    text = b""
    for block in _blocks(uri, 0):
        text += block
        i = text.find(b"\n#CHROM")
        j = text.find(b"\n", i + 1) if i >= 0 else -1
        if j >= 0:
            return text[: j + 1].decode()
    raise ValueError("no #CHROM line found")


def region_lines(uri: str, tbi: Tbi, contig: str, beg: int, end: int) -> Iterator[bytes]:
    """VCF records overlapping 0-based half-open [beg, end), as tabix would return them."""
    voff = tbi.start_offset(contig, beg, end)
    if voff is None:
        return
    name, tail, first = contig.encode(), b"", True
    for block in _blocks(uri, voff >> 16):
        if first:
            block, first = block[voff & 0xFFFF :], False
        *lines, tail = (tail + block).split(b"\n")
        for line in lines:
            chrom, pos, _, ref, _ = line.split(b"\t", 4)
            start = int(pos) - 1
            if chrom != name or start >= end:
                return
            if start + len(ref) > beg:
                yield line + b"\n"


def write_regions(uri: str, tbi_path: Path, contig: str, intervals: list[tuple[int, int]], dest: Path, head: str) -> int:
    """Write the header and every record overlapping `intervals` to `dest` (gzip). Returns the record count."""
    tbi = Tbi(tbi_path)
    n = 0
    with gzip.open(dest, "wb", compresslevel=1) as out:
        out.write(head.encode())
        for beg, end in intervals:
            for line in region_lines(uri, tbi, contig, beg, end):
                out.write(line)
                n += 1
    return n
