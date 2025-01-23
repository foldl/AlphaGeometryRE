# Copyright 2025 github/foldl
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

# Utilizations for graph representation.
# Mainly for listing combinations and permutations of elements.

import geometry as gm
import std/enumerate
from std/sequtils import toSeq

iterator cross_iter[T1, T2](elems1: openArray[T1]; elems2: openArray[T2]): tuple[a: T1, b: T2] {.closure.} =
    for e1 in elems1:
        for e2 in elems2:
            yield (e1, e2)

proc cross*[T1, T2](elems1: openArray[T1]; elems2: openArray[T2]): seq[tuple[a: T1, b: T2]] =
    result = cross_iter[T1, T2](elems1, elems2).toSeq()

iterator comb2_iter[T](elems: openArray[T]): tuple[a: T, b: T] {.closure.} =
    if len(elems) < 2:
        return
    for i, e1 in enumerate(elems[0 ..< len(elems) - 1]):
        for e2 in elems[i + 1 ..< len(elems)]:
            yield (e1, e2)

proc comb2*[T](elems: openArray[T]): seq[tuple[a: T, b: T]] =
    result = comb2_iter(elems).toSeq()

iterator comb3_iter[T](elems: openArray[T]): tuple[a: T, b: T, c: T] {.closure.} =
    if len(elems) < 3:
        return
    for i, e1 in enumerate(elems[0 ..< len(elems) - 2]):
        for j, e2 in enumerate(elems[i + 1 ..< len(elems) - 1]):
            for e3 in elems[i + j + 2 ..< len(elems)]:
                yield (e1, e2, e3)

proc comb3*[T](elems: openArray[T]): seq[tuple[a: T, b: T, c: T]] =
    result = comb3_iter(elems).toSeq()

iterator comb4_iter[T](elems: openArray[T]): tuple[a: T, b: T, c: T, d: T] {.closure.} =
    if len(elems) < 4:
        return
    for i, e1 in enumerate(elems[0 ..< len(elems) - 3]):
        for j, e2 in enumerate(elems[i + 1 ..< len(elems) - 2]):
            for k, e3 in enumerate(elems[i + j + 2 ..< len(elems) - 1]):
                for e4 in elems[i + j + k + 3 ..< len(elems)]:
                    yield (e1, e2, e3, e4)

proc comb4*[T](elems: openArray[T]): seq[tuple[a: T, b: T, c: T, d: T]] =
    result = comb4_iter(elems).toSeq()

iterator perm2_iter[T](elems: openArray[T]): tuple[a: T, b: T] {.closure.} =
    for (e1, e2) in comb2(elems):
        yield (e1, e2)
        yield (e2, e1)

proc perm2*[T](elems: openArray[T]): seq[tuple[a: T, b: T]] =
    result = perm2_iter(elems).toSeq()

iterator all_4points_iter[DepsT](l1, l2: Node[DepsT]): tuple[a: Point[DepsT], b: Point[DepsT], c: Point[DepsT], d: Point[DepsT]] {.closure.} =
    let p1s = l1.neighbors(Point[DepsT])
    let p2s = l2.neighbors(Point[DepsT])
    for a, b in perm2(p1s):
        for c, d in perm2(p2s):
            yield (a, b, c, d)

proc all_4points*[DepsT](l1, l2: Node[DepsT]): seq[tuple[a: Point[DepsT], b: Point[DepsT], c: Point[DepsT], d: Point[DepsT]]] =
    return all_4points_iter(l1, l2).toSeq()

iterator all_8points_iter[DepsT](l1, l2, l3, l4: Node[DepsT]): tuple[a: Point[DepsT], b: Point[DepsT], c: Point[DepsT], d: Point[DepsT], e: Point[DepsT], f: Point[DepsT], g: Point[DepsT], h: Point[DepsT]] {.closure.} =
    for a, b, c, d in all_4points(l1, l2):
        for e, f, g, h in all_4points(l3, l4):
            yield (a, b, c, d, e, f, g, h)

proc all_8points*[DepsT](l1, l2, l3, l4: Node[DepsT]): seq[tuple[a: Point[DepsT], b: Point[DepsT], c: Point[DepsT], d: Point[DepsT], e: Point[DepsT], f: Point[DepsT], g: Point[DepsT], h: Point[DepsT]]] =
    return all_8points_iter(l1, l2, l3, l4).toSeq()

iterator perm3_iter[T](elems: openArray[T]): tuple[a: T, b: T, c: T] {.closure.} =
    var i = 0
    while i < len(elems):
        var j = 0
        while j < len(elems):
            if j != i:
                var k = 0
                while k < len(elems):
                    if (k != i) and (k != j):
                        yield (elems[i], elems[j], elems[k])
                    inc(k)
            inc(j)
        inc(i)

proc perm3*[T](elems: openArray[T]): seq[tuple[a: T, b: T, c: T]] =
    return perm3_iter(elems).toSeq()

iterator perm4_iter[T](elems: openArray[T]): tuple[a: T, b: T, c: T, d: T] {.closure.} =
    var i = 0
    while i < len(elems):
        var j = 0
        while j < len(elems):
            if j != i:
                var k = 0
                while k < len(elems):
                    if (k != i) and (k != j):
                        var l = 0
                        while l < len(elems):
                            if (l != i) and (l != j) and (l != k):
                                yield (elems[i], elems[j], elems[k], elems[l])
                            inc(l)
                    inc(k)
            inc(j)
        inc(i)

proc perm4*[T](elems: openArray[T]): seq[tuple[a: T, b: T, c: T, d: T]] =
    return perm4_iter(elems).toSeq()