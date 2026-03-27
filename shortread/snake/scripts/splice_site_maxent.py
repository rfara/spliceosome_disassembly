#!/usr/bin/env python3

import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):
        def decorator(function):
            return function

        return decorator


BASE_CODE_STR = {"A": 0, "C": 1, "G": 2, "T": 3}
BASE_TO_CODE = np.full(256, -1, dtype=np.int8)
for base, code in BASE_CODE_STR.items():
    BASE_TO_CODE[ord(base)] = code
    BASE_TO_CODE[ord(base.lower())] = code

BGD5 = np.array([0.27, 0.23, 0.23, 0.27], dtype=np.float64)
CONS1_5 = np.array([0.004, 0.0032, 0.9896, 0.0032], dtype=np.float64)
CONS2_5 = np.array([0.0034, 0.0039, 0.0042, 0.9884], dtype=np.float64)

BGD3 = np.array([0.27, 0.23, 0.23, 0.27], dtype=np.float64)
CONS1_3 = np.array([0.9903, 0.0032, 0.0034, 0.0030], dtype=np.float64)
CONS2_3 = np.array([0.0027, 0.0037, 0.9905, 0.0030], dtype=np.float64)


def resource_dir():
    return Path(__file__).resolve().parents[1] / "resources" / "maxent"


def hash_seq_str(sequence):
    value = 0
    for base in sequence:
        value = value * 4 + BASE_CODE_STR[base]
    return value


@lru_cache(maxsize=1)
def load_models():
    score5_path = resource_dir() / "score5_matrix.txt"
    score3_path = resource_dir() / "score3_matrix.txt"

    matrix5 = np.full(4**7, np.nan, dtype=np.float64)
    with score5_path.open() as handle:
        for line in handle:
            key, value = line.split()
            matrix5[hash_seq_str(key)] = float(value)
    if np.isnan(matrix5).any():
        raise ValueError("Missing entries in MaxEnt 5'SS matrix")

    matrix3_dict = defaultdict(dict)
    with score3_path.open() as handle:
        for line in handle:
            index, key, value = line.split()
            matrix3_dict[int(index)][int(key)] = float(value)

    matrix3_arrays = []
    for index in range(9):
        mapping = matrix3_dict.get(index, {})
        if not mapping:
            raise ValueError(f"Empty MaxEnt 3'SS matrix entry {index}")
        max_key = max(mapping)
        array = np.full(max_key + 1, np.nan, dtype=np.float64)
        for key, value in mapping.items():
            array[key] = value
        if np.isnan(array).any():
            raise ValueError(f"Missing values in MaxEnt 3'SS matrix entry {index}")
        matrix3_arrays.append(array)

    return {
        "matrix5": matrix5,
        "matrix3": tuple(matrix3_arrays),
    }


def seq_to_codes(sequence):
    if sequence is None:
        return None
    codes = BASE_TO_CODE[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    if np.any(codes < 0):
        return None
    return codes


@njit(cache=True)
def _score5_window(codes, matrix5, cons1_5, cons2_5, bgd5):
    key0 = codes[3]
    key1 = codes[4]
    score = cons1_5[key0] * cons2_5[key1] / (bgd5[key0] * bgd5[key1])
    hashed = 0
    hashed = hashed * 4 + codes[0]
    hashed = hashed * 4 + codes[1]
    hashed = hashed * 4 + codes[2]
    hashed = hashed * 4 + codes[5]
    hashed = hashed * 4 + codes[6]
    hashed = hashed * 4 + codes[7]
    hashed = hashed * 4 + codes[8]
    return np.log2(score * matrix5[hashed])


@njit(cache=True)
def _hash_rest_window(codes, rest_start, length):
    hashed = 0
    for offset in range(length):
        rest_index = rest_start + offset
        if rest_index < 18:
            code = codes[rest_index]
        else:
            code = codes[rest_index + 2]
        hashed = hashed * 4 + code
    return hashed


@njit(cache=True)
def _score3_window(codes, m0, m1, m2, m3, m4, m5, m6, m7, m8, cons1_3, cons2_3, bgd3):
    key0 = codes[18]
    key1 = codes[19]
    score = cons1_3[key0] * cons2_3[key1] / (bgd3[key0] * bgd3[key1])
    rest_score = 1.0
    rest_score *= m0[_hash_rest_window(codes, 0, 7)]
    rest_score *= m1[_hash_rest_window(codes, 7, 7)]
    rest_score *= m2[_hash_rest_window(codes, 14, 7)]
    rest_score *= m3[_hash_rest_window(codes, 4, 7)]
    rest_score *= m4[_hash_rest_window(codes, 11, 7)]
    rest_score /= m5[_hash_rest_window(codes, 4, 3)]
    rest_score /= m6[_hash_rest_window(codes, 7, 4)]
    rest_score /= m7[_hash_rest_window(codes, 11, 3)]
    rest_score /= m8[_hash_rest_window(codes, 14, 4)]
    return np.log2(score * rest_score)


def score5(sequence):
    if sequence is None or len(sequence) != 9:
        return math.nan
    codes = seq_to_codes(sequence)
    if codes is None:
        return math.nan
    models = load_models()
    return float(_score5_window(codes, models["matrix5"], CONS1_5, CONS2_5, BGD5))


def score3(sequence):
    if sequence is None or len(sequence) != 23:
        return math.nan
    codes = seq_to_codes(sequence)
    if codes is None:
        return math.nan
    models = load_models()
    return float(_score3_window(codes, *models["matrix3"], CONS1_3, CONS2_3, BGD3))
