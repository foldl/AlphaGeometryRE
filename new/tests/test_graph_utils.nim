discard """
  joinable: false
"""
import sets
import std/[algorithm, strformat]

import ../graph_utils as gu

proc assertEqual[A, B](a: A; b: B) =
    assert a == b, fmt"{a} != {b}"

proc test_cross() =
    assert gu.cross[int, int](@[], @[1]) == @[]
    assert gu.cross[int, int](@[1], @[]) == @[]
    assert gu.cross(@[1], @[2]) == @[(1, 2)]
    assert gu.cross(@[1], @[2, 3]) == @[(1, 2), (1, 3)]

    let e1 = @[1, 2, 3]
    let e2 = @[4, 5]
    let target = @[(1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5)]
    assert gu.cross(e1, e2) == target

func test_comb2() =
    assertEqual(gu.comb2[int]([]).len(), 0)
    assertEqual(gu.comb2([1]).len(), 0)
    assertEqual(gu.comb2([1, 2]), @[(1, 2)])
    assertEqual(gu.comb2([1, 2, 3]), @[(1, 2), (1, 3), (2, 3)])

func test_comb3() =
    assertEqual(gu.comb3[int]([]).len(), 0)
    assertEqual(gu.comb3([1]).len(), 0)
    assertEqual(gu.comb3([1, 2]).len(), 0)
    assertEqual(gu.comb3([1, 2, 3]), @[(1, 2, 3)])
    assertEqual(
        gu.comb3([1, 2, 3, 4]), @[(1, 2, 3),
                                  (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    )

proc test_comb4() =
    assertEqual(gu.comb4[int]([]).len(), 0)
    assertEqual(gu.comb4([1]).len(), 0)
    assertEqual(gu.comb4([1, 2]).len(), 0)
    assertEqual(gu.comb4([1, 2, 3]).len(), 0)
    assertEqual(gu.comb4([1, 2, 3, 4]), @[(1, 2, 3, 4)])
    assertEqual(
        gu.comb4([1, 2, 3, 4, 5]),
        @[(1, 2, 3, 4), (1, 2, 3, 5), (1, 2, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5)],
    )

proc test_perm2() =
    assertEqual(gu.perm2[int]([]).len(), 0)
    assertEqual(gu.perm2([1]).len(), 0)
    assertEqual(gu.perm2([1, 2]), [(1, 2), (2, 1)])
    assertEqual(
        gu.perm2([1, 2, 3]), [(1, 2), (2, 1),
                                (1, 3), (3, 1), (2, 3), (3, 2)]
    )

proc test_perm3() =
    assertEqual(gu.perm3[int]([]).len(), 0)
    assertEqual(gu.perm3([1]).len(), 0)
    assertEqual(gu.perm3([1, 2]).len(), 0)
    assertEqual(
        gu.perm3([1, 2, 3]),
        @[(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)],
    )
    assertEqual(
        gu.perm3([1, 2, 3, 4]),
        @[
            (1, 2, 3),
            (1, 2, 4),
            (1, 3, 2),
            (1, 3, 4),
            (1, 4, 2),
            (1, 4, 3),
            (2, 1, 3),
            (2, 1, 4),
            (2, 3, 1),
            (2, 3, 4),
            (2, 4, 1),
            (2, 4, 3),
            (3, 1, 2),
            (3, 1, 4),
            (3, 2, 1),
            (3, 2, 4),
            (3, 4, 1),
            (3, 4, 2),
            (4, 1, 2),
            (4, 1, 3),
            (4, 2, 1),
            (4, 2, 3),
            (4, 3, 1),
            (4, 3, 2),
        ],
    )

proc test_perm4() =
    assertEqual(gu.perm4[int]([]).len(), 0)
    assertEqual(gu.perm4([1]).len(), 0)
    assertEqual(gu.perm4([1, 2]).len(), 0)
    assertEqual(gu.perm4([1, 2, 3]).len(), 0)
    assertEqual(
        gu.perm4([1, 2, 3, 4]),
        @[
            (1, 2, 3, 4),
            (1, 2, 4, 3),
            (1, 3, 2, 4),
            (1, 3, 4, 2),
            (1, 4, 2, 3),
            (1, 4, 3, 2),
            (2, 1, 3, 4),
            (2, 1, 4, 3),
            (2, 3, 1, 4),
            (2, 3, 4, 1),
            (2, 4, 1, 3),
            (2, 4, 3, 1),
            (3, 1, 2, 4),
            (3, 1, 4, 2),
            (3, 2, 1, 4),
            (3, 2, 4, 1),
            (3, 4, 1, 2),
            (3, 4, 2, 1),
            (4, 1, 2, 3),
            (4, 1, 3, 2),
            (4, 2, 1, 3),
            (4, 2, 3, 1),
            (4, 3, 1, 2),
            (4, 3, 2, 1),
        ],
    )

proc main() =
    test_cross()
    test_comb2()
    test_comb3()
    test_comb4()
    test_perm2()
    test_perm3()
    test_perm4()

main()