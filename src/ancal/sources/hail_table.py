"""Minimal reader for small Hail native Tables keyed by a single string field.

Koenig et al. (2024) publish the HGDP+1KG related-sample list only as a Hail Table
(pca_preprocessing/related_sample_ids.ht). Installing Hail plus Java to read about
700 IDs would add roughly 700 MB, so this module decodes exactly one layout, and
refuses anything else:

  row type  Struct{<field>:String}, encoded +EBaseStruct{<field>:EBinary}
  buffers   LEB128 <- Blocking(65536) <- ZstdBlock(65536) <- StreamBlock

Encoding (Hail 0.2.x):
  block    := int32 LE len, int32 LE decompressed_len, zstd frame   (len = 4 + frame size)
  stream   := (0x01 row)* 0x00
  row      := missing-bits byte (1 optional field) + EBinary
  EBinary  := LEB128 length + bytes
"""

from __future__ import annotations

import gzip
import hashlib
import json
import struct

import requests
import zstandard

EXPECTED_BUFFER_SPEC = {
    "name": "LEB128BufferSpec",
    "child": {
        "name": "BlockingBufferSpec",
        "blockSize": 65536,
        "child": {
            "name": "ZstdBlockBufferSpec",
            "blockSize": 65536,
            "child": {"name": "StreamBlockBufferSpec"},
        },
    },
}


class UnsupportedHailTable(ValueError):
    pass


def _get(url: str) -> bytes:
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.content


def _decompress_blocks(data: bytes) -> bytes:
    out, pos, dctx = bytearray(), 0, zstandard.ZstdDecompressor()
    while pos < len(data):
        if pos + 8 > len(data):
            raise UnsupportedHailTable("truncated block header")
        # clen counts the 4-byte decompressed-length field plus the zstd frame
        clen, dlen = struct.unpack_from("<ii", data, pos)
        frame = data[pos + 8 : pos + 4 + clen]
        pos += 4 + clen
        block = dctx.decompress(frame, max_output_size=dlen)
        if len(block) != dlen:
            raise UnsupportedHailTable(f"block decompressed to {len(block)} bytes, expected {dlen}")
        out += block
    return bytes(out)


def _leb128(buf: bytes, pos: int) -> tuple[int, int]:
    result = shift = 0
    while True:
        b = buf[pos]
        pos += 1
        result |= (b & 0x7F) << shift
        if not b & 0x80:
            return result, pos
        shift += 7


def decode_string_rows(raw: bytes) -> list[str | None]:
    """Decode one partition file of a single-optional-string-field table."""
    buf, pos, rows = _decompress_blocks(raw), 0, []
    while True:
        if pos >= len(buf):
            raise UnsupportedHailTable("stream ended without end-of-rows marker")
        marker = buf[pos]
        pos += 1
        if marker == 0:
            break
        if marker != 1:
            raise UnsupportedHailTable(f"unexpected row marker {marker}")
        missing = buf[pos]
        pos += 1
        if missing & 1:
            rows.append(None)
            continue
        n, pos = _leb128(buf, pos)
        rows.append(buf[pos : pos + n].decode())
        pos += n
    if pos != len(buf):
        raise UnsupportedHailTable(f"{len(buf) - pos} trailing bytes after end-of-rows marker")
    return rows


def read_string_key_table(base_url: str) -> tuple[list[str], dict]:
    """Read every row of a remote string-keyed Hail Table. Returns (values, provenance info)."""
    base_url = base_url.rstrip("/")
    top = json.loads(gzip.decompress(_get(f"{base_url}/metadata.json.gz")))
    rows_meta = json.loads(gzip.decompress(_get(f"{base_url}/rows/metadata.json.gz")))

    codec = rows_meta["_codecSpec"]
    fields = rows_meta["_key"]
    if len(fields) != 1:
        raise UnsupportedHailTable(f"expected a single key field, got {fields}")
    field = fields[0]
    if codec.get("_vType") != f"Struct{{{field}:String}}" or codec.get("_eType") != f"+EBaseStruct{{{field}:EBinary}}":
        raise UnsupportedHailTable(f"unsupported row type {codec.get('_vType')} / {codec.get('_eType')}")
    if codec.get("_bufferSpec") != EXPECTED_BUFFER_SPEC:
        raise UnsupportedHailTable("unsupported buffer spec; install Hail to read this table")

    expected_counts = top["components"]["partition_counts"]["counts"]
    parts = rows_meta["_partFiles"]
    if len(parts) != len(expected_counts):
        raise UnsupportedHailTable("partition file count does not match partition_counts")

    values: list[str] = []
    digest = hashlib.sha256()  # over the partition files, in order
    for part, expected in zip(parts, expected_counts):
        raw = _get(f"{base_url}/rows/parts/{part}")
        digest.update(raw)
        rows = decode_string_rows(raw)
        if len(rows) != expected:
            raise UnsupportedHailTable(f"{part}: decoded {len(rows)} rows, metadata says {expected}")
        if any(r is None for r in rows):
            raise UnsupportedHailTable(f"{part}: missing key values")
        values.extend(rows)
    info = {
        "hail_version": top.get("hail_version", ""),
        "row_count": sum(expected_counts),
        "key_field": field,
        "parts_sha256": digest.hexdigest(),
    }
    return values, info
