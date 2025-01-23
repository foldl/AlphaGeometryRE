discard """
  joinable: false
"""
import sets
import std/[algorithm, strformat, math, random, options]

import ../numericals as nm
import ../n_utils

proc assertEqual[A, B](a: A; b: B) =
    assert a == b, fmt"{a} != {b}"

proc assertTrue(a: bool) =
    assert a, "not true"

proc assertTrue(a: Option[bool]) =
    assertTrue a.get(false)

proc assertFalse(a: bool) =
    assert not a, "not false"

proc assertFalse(a: Option[bool]) =
    assertFalse a.get(true)

proc assertAlmostEqual(a, b: float, tol: float = 1e-10) =
    assert abs(a - b) < tol, fmt"{a} != {b}"

proc test_sketch_ieq_triangle() =
    let (a, b, c) = nm.sketch_ieq_triangle([])
    assertAlmostEqual(a.distance(b), b.distance(c))
    assertAlmostEqual(c.distance(a), b.distance(c))

proc test_sketch_2l1c() =
    let p = nm.newPoint(0.0, 0.0)
    let pi = math.PI
    let anga = rand(-0.4 * pi .. 0.4 * pi)
    let a = newPoint(math.cos(anga), math.sin(anga))
    let angb = rand(0.6 * pi .. 1.4 * pi)
    let b = newPoint(math.cos(angb), math.sin(angb))

    let angc = rand(anga + 0.05 * pi .. angb - 0.05 * pi)
    let c = newPoint(math.cos(angc), math.sin(angc)) * rand(0.2 .. 0.8)

    let (x, y, z, i) = nm.sketch_2l1c([a, b, c, p])
    assert(check_coll([x, c, a]))
    assert(check_coll([y, c, b]))
    assertAlmostEqual(z.distance(p), 1.0)
    assert(check_coll([p, i, z]))
    assert(newLine(i, x).is_perp(newLine(c, a)))
    assert(newLine(i, y).is_perp(newLine(c, b)))
    assertAlmostEqual(i.distance(x), i.distance(y))
    assertAlmostEqual(i.distance(x), i.distance(z))

proc test_sketch_3peq() =
    let (a, b, c) = seq2tuple3 random_points(3)
    let (x, y, z) = nm.sketch_3peq([a, b, c])

    assert(check_coll([a, b, x]))
    assert(check_coll([a, c, y]))
    assert(check_coll([b, c, z]))
    assert(check_coll([x, y, z]))
    assertAlmostEqual(z.distance(x), z.distance(y))

proc test_sketch_aline() =
    let (a, b, c, d, e) = seq2tuple5(random_points(5))
    let ex = nm.sketch_aline([a, b, c, d, e])
    assert(ex of HalfLine)
    assertEqual(HalfLine(ex).tail, e)
    let x = HalfLine(ex).head
    assertAlmostEqual(ang_between(b, a, c), ang_between(e, d, x))

proc test_sketch_amirror() =
    let (a, b, c) = seq2tuple3(random_points(3))
    let bx0 = nm.sketch_amirror([a, b, c])
    assert(bx0 of HalfLine)

    let bx = HalfLine(bx0)
    assert bx.tail == b
    let x = bx.head

    let ang1 = ang_between(b, a, c)
    let ang2 = ang_between(b, c, x)
    assertAlmostEqual(ang1, ang2)

proc test_sketch_bisectr() =
    let (a, b, c) = seq2tuple3(random_points(3))
    let line = nm.sketch_bisect([a, b, c])
    assertAlmostEqual(b.distance(line), 0.0)

    let l = a.perpendicular_line(line)
    let x = line_line_intersection(l, newLine(b, c))
    assertAlmostEqual(a.distance(line), x.distance(line))

    let (d, _) = line_circle_intersection(line, newCircle(b, radius=1))
    let ang1 = ang_between(b, a, d)
    let ang2 = ang_between(b, d, c)
    assertAlmostEqual(ang1, ang2)

proc test_sketch_bline() =
    let (a, b) = seq2tuple2(random_points(2))
    let l = nm.sketch_bline([a, b])
    assert(newLine(a, b).is_perp(l))
    assertAlmostEqual(a.distance(l), b.distance(l))

proc test_sketch_cc_tangent() =
    let o = newPoint(0.0, 0.0)
    let w = newPoint(1.0, 0.0)

    let ra = rand(0.0 .. math.PI)
    let rb = rand(0.0 .. math.PI)

    let a = o + rand(0.0 .. 0.6) * newPoint(math.cos(ra), math.sin(ra))
    let b = w + rand(0.4 .. 1.0) * newPoint(math.sin(rb), math.cos(rb))

    let (x, y, z, t) = nm.sketch_cc_tangent([o, a, w, b])
    let xy = newLine(x, y)
    let zt = newLine(z, t)
    assertAlmostEqual(o.distance(xy), o.distance(a))
    assertAlmostEqual(o.distance(zt), o.distance(a))
    assertAlmostEqual(w.distance(xy), w.distance(b))
    assertAlmostEqual(w.distance(zt), w.distance(b))

proc test_sketch_circle() =
    let (a, b, c) = seq2tuple3(random_points(3))
    let circle = nm.sketch_circle([a, b, c])
    assertAlmostEqual(circle.center.distance(a), 0.0)
    assertAlmostEqual(circle.radius, b.distance(c))

proc test_sketch_e5128() =
    let b = newPoint(0.0, 0.0)
    let c = newPoint(0.0, 1.0)
    let ang = rand(-math.PI / 2 .. 3 * math.PI / 2)
    let d = head_from(c, ang, 1.0)
    let a = newPoint(rand(0.5 .. 2.0), 0.0)

    let (e, g) = nm.sketch_e5128([a, b, c, d])
    let ang1 = ang_between(a, b, d)
    let ang2 = ang_between(e, a, g)
    assertAlmostEqual(ang1, ang2)

proc test_sketch_eq_quadrangle() =
    let (a, b, c, d) = nm.sketch_eq_quadrangle([])
    assertAlmostEqual(a.distance(d), c.distance(b))
    let ac = newLine(a, c)
    assert ac.diff_side(b, d).get(false), fmt"{(ac.at(b), ac.at(d))}"
    let bd = newLine(b, d)
    assert bd.diff_side(a, c).get(false), fmt"{(bd.at(a), bd.at(c))}"

proc test_sketch_eq_trapezoid() =
    let (a, b, c, d) = nm.sketch_eq_trapezoid([])
    assert newLine(a, b).is_parallel(newLine(c, d))
    assertAlmostEqual(a.distance(d), b.distance(c))

proc test_sketch_eqangle3() =
    let points = random_points(5)
    let x = nm.sketch_eqangle3(points).sample_within(points)[0]
    let (a, b, d, e, f) = seq2tuple5(points)
    assertTrue(check_eqangle([x, a, x, b, d, e, d, f]))

proc test_sketch_eqangle2() =
    let (a, b, c) = seq2tuple3(random_points(3))
    let x = nm.sketch_eqangle2([a, b, c])
    let ang1 = ang_between(a, b, x)
    let ang2 = ang_between(c, x, b)
    assertAlmostEqual(ang1, ang2)

proc test_sketch_edia_quadrangle() =
    let (a, b, c, d) = nm.sketch_eqdia_quadrangle([])
    assert newLine(a, c).diff_side(b, d).get(false)
    assert newLine(b, d).diff_side(a, c).get(false)
    assertAlmostEqual(a.distance(c), b.distance(d))

proc test_sketch_isos() =
    let (a, b, c) = nm.sketch_isos([])
    assertAlmostEqual(a.distance(b), a.distance(c))
    assertAlmostEqual(ang_between(b, a, c), ang_between(c, b, a))

proc test_sketch_quadrange() =
    let (a, b, c, d) = nm.sketch_quadrangle([])
    assertTrue(newLine(a, c).diff_side(b, d))
    assertTrue(newLine(b, d).diff_side(a, c))

proc test_sketch_r_trapezoid() =
    let (a, b, c, d) = nm.sketch_r_trapezoid([])
    assertTrue(newLine(a, b).is_perp(newLine(a, d)))
    assertTrue(newLine(a, b).is_parallel(newLine(c, d)))
    assertTrue(newLine(a, c).diff_side(b, d))
    assertTrue(newLine(b, d).diff_side(a, c))

proc test_sketch_r_triangle() =
    let (a, b, c) = nm.sketch_r_triangle([])
    assertTrue(newLine(a, b).is_perp(newLine(a, c)))

proc test_sketch_rectangle() =
    let (a, b, c, d) = nm.sketch_rectangle([])
    assertTrue(newLine(a, b).is_perp(newLine(b, c)))
    assertTrue(newLine(b, c).is_perp(newLine(c, d)))
    assertTrue(newLine(c, d).is_perp(newLine(d, a)))

proc test_sketch_reflect() =
    let (a, b, c) = seq2tuple3(random_points(3))
    let x = nm.sketch_reflect([a, b, c])
    assertTrue(newLine(a, x).is_perp(newLine(b, c)))
    assertAlmostEqual(x.distance(newLine(b, c)), a.distance(newLine(b, c)))

proc test_sketch_risos() =
    let (a, b, c) = nm.sketch_risos([])
    assertAlmostEqual(a.distance(b), a.distance(c))
    assertTrue(newLine(a, b).is_perp(newLine(a, c)))

proc test_sketch_rotaten90() =
    let (a, b) = seq2tuple2(random_points(2))
    let x = nm.sketch_rotaten90([a, b])
    assertAlmostEqual(a.distance(x), a.distance(b))
    assertTrue(newLine(a, x).is_perp(newLine(a, b)))
    let d = newPoint(0.0, 0.0)
    let e = newPoint(0.0, 1.0)
    let f = newPoint(1.0, 0.0)
    assertAlmostEqual(ang_between(d, e, f), ang_between(a, b, x))

proc test_sketch_rotatep90() =
    let (a, b) = seq2tuple2(random_points(2))
    let x = nm.sketch_rotatep90([a, b])
    assertAlmostEqual(a.distance(x), a.distance(b))
    assertTrue(newLine(a, x).is_perp(newLine(a, b)))
    let d = newPoint(0.0, 0.0)
    let e = newPoint(0.0, 1.0)
    let f = newPoint(1.0, 0.0)
    assertAlmostEqual(ang_between(d, f, e), ang_between(a, b, x))

proc test_sketch_s_angle() =
    let (a, b) = seq2tuple2(random_points(2))
    let y = rand(0.0 .. math.PI)
    let bx0 = nm.sketch_s_angle([a, b], [y / math.PI * 180])
    assert(bx0 of HalfLine)

    let bx = HalfLine(bx0)
    assertEqual(bx.tail, b)
    let x = bx.head

    let d = newPoint(1.0, 0.0)
    let e = newPoint(0.0, 0.0)
    let f = newPoint(math.cos(y), math.sin(y))
    assertAlmostEqual(ang_between(e, d, f), ang_between(b, a, x))

proc test_sketch_shift() =
    let (a, b, c) = seq2tuple3(random_points(3))
    let x = nm.sketch_shift([a, b, c])
    assertTrue((b - a).close(x - c))

proc test_sketch_square() =
    let (a, b) = seq2tuple2(random_points(2))
    let (c, d) = nm.sketch_square([a, b])
    assertTrue(newLine(a, b).is_perp(newLine(b, c)))
    assertTrue(newLine(b, c).is_perp(newLine(c, d)))
    assertTrue(newLine(c, d).is_perp(newLine(d, a)))
    assertAlmostEqual(a.distance(b), b.distance(c))

proc test_sketch_isquare() =
    let (a, b, c, d) = nm.sketch_isquare([])
    assertTrue(newLine(a, b).is_perp(newLine(b, c)))
    assertTrue(newLine(b, c).is_perp(newLine(c, d)))
    assertTrue(newLine(c, d).is_perp(newLine(d, a)))
    assertAlmostEqual(a.distance(b), b.distance(c))

proc test_sketch_trapezoid() =
    let (a, b, c, d) = nm.sketch_trapezoid([])
    assertTrue(newLine(a, b).is_parallel(newLine(c, d)))
    assertTrue(newLine(a, c).diff_side(b, d))
    assertTrue(newLine(b, d).diff_side(a, c))

proc test_sketch_triangle() =
    let (a, b, c) = nm.sketch_triangle([])
    assertFalse(check_coll([a, b, c]))

proc test_sketch_triangle12() =
    let (a, b, c) = nm.sketch_triangle12([])
    assertAlmostEqual(a.distance(b) * 2, a.distance(c))

proc test_sketch_trisect() =
    let (a, b, c) = seq2tuple3(random_points(3))
    let (x, y) = nm.sketch_trisect([a, b, c])
    assertAlmostEqual(ang_between(b, a, x), ang_between(b, x, y))
    assertAlmostEqual(ang_between(b, x, y), ang_between(b, y, c))
    assertAlmostEqual(ang_between(b, a, x) * 3, ang_between(b, a, c))

proc test_sketch_trisegment() =
    let (a, b) = seq2tuple2(random_points(2))
    let (x, y) = nm.sketch_trisegment([a, b])
    assertAlmostEqual(a.distance(x) + x.distance(y) + y.distance(b), a.distance(b))
    assertAlmostEqual(a.distance(x), x.distance(y))
    assertAlmostEqual(x.distance(y), y.distance(b))

proc main() =
    test_sketch_ieq_triangle()
    test_sketch_2l1c()
    test_sketch_3peq()
    test_sketch_aline()
    test_sketch_amirror()
    test_sketch_bisectr()
    test_sketch_bline()
    test_sketch_cc_tangent()
    test_sketch_circle()
    test_sketch_e5128()
    test_sketch_eq_quadrangle()
    test_sketch_eq_trapezoid()
    test_sketch_eqangle3()
    test_sketch_eqangle2()
    test_sketch_edia_quadrangle()
    test_sketch_isos()
    test_sketch_quadrange()
    test_sketch_r_trapezoid()
    test_sketch_r_triangle()
    test_sketch_rectangle()
    test_sketch_reflect()
    test_sketch_risos()
    test_sketch_rotaten90()
    test_sketch_rotatep90()
    test_sketch_s_angle()
    test_sketch_shift()
    test_sketch_square()
    test_sketch_isquare()
    test_sketch_trapezoid()
    test_sketch_triangle()
    test_sketch_triangle12()
    test_sketch_trisect()
    test_sketch_trisegment()

main()
