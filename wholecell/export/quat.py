"""Quaternion helpers. wholecell order everywhere on disk: (w, x, y, z), scalar first, canonical sign w >= 0
(docs/reference/tile_format_v2.md section 2). Scalar-last appears only at third-party boundaries (scipy, three.js, WGSL)."""
from __future__ import annotations

import numpy as np

INV_SQRT2 = 1.0 / np.sqrt(2.0)


def to_xyzw(q):
    q = np.asarray(q)
    return q[..., [1, 2, 3, 0]]


def from_xyzw(q):
    q = np.asarray(q)
    return q[..., [3, 0, 1, 2]]


def canonical(q):
    """unit quaternions with w >= 0 (q and -q are the same rotation)"""
    q = np.asarray(q, np.float64)
    q = q / np.linalg.norm(q, axis=-1, keepdims=True)
    return np.where(q[..., :1] < 0, -q, q)


def matrix(q):
    """rotation matrices (n,3,3) of (w,x,y,z) quaternions (world = R . local + p)"""
    a, b, c, d = (np.asarray(q, np.float64)[..., i] for i in range(4))
    return np.stack([np.stack([a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], -1),
                     np.stack([2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)], -1),
                     np.stack([2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d], -1)], -2)


def encode_smallest3(q):
    """(n,4) w,x,y,z -> (qidx uint32 (n,), three int16 (n,3)): the largest |component| is dropped (made positive),
    the other three are stored as int16 = round(c * sqrt2 * 32767)"""
    q = np.asarray(q, np.float64)
    q = q / np.linalg.norm(q, axis=1, keepdims=True)
    idx = np.argmax(np.abs(q), axis=1)
    sgn = np.where(q[np.arange(len(q)), idx] < 0, -1.0, 1.0)
    q = q * sgn[:, None]
    keep = np.array([[j for j in range(4) if j != i] for i in range(4)])[idx]
    rest = np.take_along_axis(q, keep, axis=1)
    i16 = np.clip(np.rint(rest / INV_SQRT2 * 32767.0), -32767, 32767).astype(np.int16)
    return idx.astype(np.uint32), i16


def decode_smallest3(idx, i16):
    rest = i16.astype(np.float64) / 32767.0 * INV_SQRT2
    n = len(idx)
    q = np.zeros((n, 4))
    keep = np.array([[j for j in range(4) if j != i] for i in range(4)])[idx]
    np.put_along_axis(q, keep, rest, axis=1)
    big = np.sqrt(np.clip(1.0 - (rest ** 2).sum(1), 0.0, None))
    q[np.arange(n), idx] = big
    return canonical(q)


def angle_between(q1, q2):
    """rotation angle (rad) between two sets of unit quaternions"""
    d = np.abs((canonical(q1) * canonical(q2)).sum(-1))
    return 2.0 * np.arccos(np.clip(d, -1.0, 1.0))
