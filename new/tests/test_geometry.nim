discard """
  joinable: false
"""
import sets
import std/[algorithm, strformat]

import ../geometry as gm

type DepsT = string

proc setup_equality_example(): tuple[a: gm.Node[string], b: Node[string], c: Node[string], d: Node[string], e: Node[string], f: Node[string], g: Node[string], h: Node[string]] =
    # Create 4 nodes a, b, c, d
    # and their lengths
    var a = gm.newNodeType[gm.Segment[string]]("a")
    var la = gm.newNodeType[gm.Length[string]]("l(a)")
    a.connect_to(la)
    la.connect_to(a)

    var b = gm.newNodeType[Segment[string]]("b")
    var lb = gm.newNodeType[Length[string]]("l(b)")
    b.connect_to(lb)
    lb.connect_to(b)

    var c = gm.newNodeType[Segment[string]]("c")
    var lc = gm.newNodeType[Length[string]]("l(c)")
    c.connect_to(lc)
    lc.connect_to(c)

    var d = gm.newNodeType[Segment[string]]("d")
    var ld = gm.newNodeType[Length[string]]("l(d)")
    d.connect_to(ld)
    ld.connect_to(d)

    # Now let a=b, b=c, a=c, c=d
    la.merge([lb], "fact1")
    lb.merge([lc], "fact2")
    la.merge([lc], "fact3")
    lc.merge([ld], "fact4")
    return (a, b, c, d, la, lb, lc, ld)

proc test_merged_node_representative() =
    let (a, b, c, d, la, lb, lc, ld) = setup_equality_example()

    # all nodes are now represented by la.
    assert la.rep() == la
    assert lb.rep() == la
    assert lc.rep() == la
    assert ld.rep() == la

proc test_merged_node_equivalence() =
    let (a, b, c, d, la, lb, lc, ld) = setup_equality_example()

    # all la, lb, lc, ld are equivalent
    assert la.equivs() == [la, lb, lc, ld].toHashSet()
    assert lb.equivs() == [la, lb, lc, ld].toHashSet()
    assert lc.equivs() == [la, lb, lc, ld].toHashSet()
    assert ld.equivs() == [la, lb, lc, ld].toHashSet()

proc test_bfs_for_equality_transitivity() =
    let (a, b, c, d, la, lb, lc, ld) = setup_equality_example()

    # check that a==d because fact3 & fact4, not fact1 & fact2
    let why = gm.why_equal(a, d)
    assert why.sorted() == @["fact3", "fact4"].sorted(), fmt"{why} != [fact3, fact4]"

proc main() =
    test_merged_node_representative()
    test_merged_node_equivalence()
    test_bfs_for_equality_transitivity()


main()



