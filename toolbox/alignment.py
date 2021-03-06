import collections
from enum import auto, Enum
from typing import (
    Tuple,
    List,
    Union,
    Mapping,
    Iterator,
)
from itertools import chain
from collections import deque


class Operation(Enum):
    NON = auto()
    ADD = auto()
    DEL = auto()
    MOD = auto()


Result = Tuple[int, Iterator[Operation]]
Memo = List[List[Union[None, Result]]]


def align(s1: Mapping, s2: Mapping, list_all=True):
    # Rows = s1, cols = s2
    memo = [[None] * (len(s2) + 1) for _ in range(len(s1) + 1)]
    l, memo = alignh(memo, s1, s2, 0, 0)
    alignments, _ = backtrack(
        [[None] * (len(s2) + 1) for _ in range(len(s1) + 1)],
        memo,
        s1,
        s2,
        0,
        0,
        list_all,
    )
    return {"distance": l, "alignments": alignments}


def alignh(memo: Memo, s1: Mapping, s2: Mapping, i1: int, i2: int) -> Tuple[int, Memo]:
    i1, i2 = min(i1, len(s1)), min(i2, len(s2))

    if memo[i1][i2] is not None:
        return memo[i1][i2][0], memo

    if i1 >= len(s1) and i2 >= len(s2):
        memo[i1][i2] = 0, []
    elif i1 >= len(s1):
        l, memo = alignh(memo, s1, s2, i1, i2 + 1)
        memo[i1][i2] = l + 1, [Operation.ADD]
    elif i2 >= len(s2):
        l, memo = alignh(memo, s1, s2, i1 + 1, i2)
        memo[i1][i2] = l + 1, [Operation.DEL]
    elif s1[i1] == s2[i2]:
        l, memo = alignh(memo, s1, s2, i1 + 1, i2 + 1)
        memo[i1][i2] = l, [Operation.NON]
    else:
        lens = []
        for op, j1, j2 in [
            (Operation.ADD, 0, 1),
            (Operation.MOD, 1, 1),
            (Operation.DEL, 1, 0),
        ]:
            l, memo = alignh(memo, s1, s2, i1 + j1, i2 + j2)
            lens.append((l, op))

        m = min(l for l, _ in lens)
        memo[i1][i2] = 1 + m, [op for l, op in lens if l == m]

    return memo[i1][i2][0], memo


# Tried to make it faster by using deques instead of strings
# Failed
def backtrack(memo, data, s1: str, s2: str, i1: int, i2: int, list_all: bool):
    if i1 == len(s1) and i2 == len(s2):
        memo[i1][i2] = [(deque(), deque())]
    else:
        possibilities = []
        for op in data[i1][i2][1] if list_all else data[i1][i2][1][0:1]:
            if op == Operation.DEL:
                j1, j2, l, r = 1, 0, s1[i1], "-"
            elif op == Operation.ADD:
                j1, j2, l, r = 0, 1, "-", s2[i2]
            else:
                j1, j2, l, r = 1, 1, s1[i1], s2[i2]
            base, memo = backtrack(memo, data, s1, s2, i1 + j1, i2 + j2, list_all)
            base = [(a.copy(), b.copy()) for a, b in base] if list_all else base
            for a, b in base:
                a.appendleft(l)
                b.appendleft(r)
            possibilities += base
        memo[i1][i2] = possibilities

    return memo[i1][i2], memo


def prepend_op(op, result):
    l, alignments = result
    return l, [chain([op], a) for a in alignments]
