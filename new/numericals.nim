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

import std/math
from std/random import rand
import std/[strformat, options, random, sequtils, sugar]

import geometry as gm
import n_utils

const ATOM = 1e-12

# Naming in geometry is a little different
# we stick to geometry naming to better read the code.

type
    BoxedFloat * = ref object of RootObj
        v: float

    Point* = ref object of RootObj
        x, y: float

    Line* = ref object of RootObj
        a, b, c: float # a x + b y + c = 0

    HalfLine* = ref object of Line
        tail*: Point
        head*: Point

    Circle* = ref object of RootObj
        center*: Point
        radius*: float
        r2: float       # sqr(radius)

    HoleCircle* = ref object of Circle
        hole: Point

proc line_line_intersection*(l1, l2: Line): Point
proc parallel_line*(self: Line, p: Point): Line
proc sign*(self: Line, point: Point): int
proc close_enough*(a, b: float; tol: float = 1e-12): bool
proc line_circle_intersection*(line: Line, circle: Circle): (Point, Point)
proc newCircle*(center: Point, radius: float): Circle
proc circle_circle_intersection*(c1, c2: Circle): (Point, Point)
proc ang_between*(tail: Point; head1, head2: Point): float

proc check_eqangle*(points: openArray[Point]): bool
proc check_eqratio*(points: openArray[Point]): bool
proc check_simtri*(points: openArray[Point]): bool
proc check_contri*(points: openArray[Point]): bool

proc sketch_bisect*(args: openArray[Point]): Line
proc sketch_exbisect*(args: openArray[Point]): Line
proc sketch_aline*(args: openArray[Point]): HalfLine

proc newPoint*(x, y: float): Point = result = Point(x: x, y: y)

proc `<`*(p1, p2: Point): bool = result = (p1.x, p1.y) < (p2.x, p2.y)

proc `>`*(p1, p2: Point): bool = result = (p1.x, p1.y) > (p2.x, p2.y)

proc `+`*(p1, p2: Point): Point = newPoint(p1.x + p2.x, p1.y + p2.y)

proc `-`*(p1, p2: Point): Point = newPoint(p1.x - p2.x, p1.y - p2.y)

proc `*`*(p: Point; f: float): Point = newPoint(p.x * f, p.y * f)

proc `*`*(f: float; p: Point): Point = p * f

proc `/`*(p: Point; f: float): Point = newPoint(p.x / f, p.y / f)

proc `div`*(p: Point; f: float): Point = newPoint(floor(p.x / f), floor(p.y / f))

proc `$`*(p: Point): string = fmt"P({p.x}, {p.y})"

proc close*(self, point: Point; tol: float = 1e-12): bool=
    return abs(self.x - point.x) < tol and abs(self.y - point.y) < tol

proc midpoint*(p1, p2: Point): Point = newPoint(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y))

proc at*(self: Line; x, y: float): float

proc distance*(self: Line; p: Point): float =
    let (a, b) = (self.a, self.b)
    return abs(self.at(p.x, p.y)) / sqrt(a * a + b * b)

proc distance*(self: Point, p: Point): float = hypot(self.x - p.x, self.y - p.y)

proc distance*(self: Point, p: Line): float = p.distance(self)

proc distance*(self: Point, p: Circle): float = abs(p.radius - self.distance(p.center))

proc rotate*(self: Point; sinb, cosb: float): Point =
    return newPoint(self.x * cosb - self.y * sinb, self.x * sinb + self.y * cosb)

proc rotate*(self: Point; ang: float): Point = self.rotate(sin(ang), cos(ang))

proc flip*(self: Point): Point = newPoint(-self.x, self.y)

proc perpendicular_line*(self: Line, p: Point): Line

proc perpendicular_line*(self: Point, line: Line): Line = line.perpendicular_line(self)

proc foot*(self: Point, line: Line): Point =
    let l = line.perpendicular_line(self)
    result = line_line_intersection(l, line)

proc foot*(self: Point, circle: Circle): Point =
    let c = circle.center
    let r = circle.radius
    result = c + (self - c) * r / self.distance(c)

proc parallel_line*(self: Point, line: Line): Line = line.parallel_line(self)

proc norm*(self: Point): float = hypot(self.x, self.y)

proc cos*(self, other: Point): float =
    let (x, y) = (self.x, self.y)
    let (a, b) = (other.x, other.y)
    return (x * a + y * b) / self.norm() / other.norm()

proc dot*(self, other: Point): float =
    return self.x * other.x + self.y * other.y

proc sign*(self: Point, line: Line): int = line.sign(self)

proc is_same(self, other: Point): bool = self.distance(other) <= ATOM

proc newLine*(a0, b0, c0: float): Line =
    var (a, b, c) = (a0, b0, c0)
    if a < 0.0 or a == 0.0 and b > 0.0:
        (a, b, c) = (-a, -b, -c)
    result = Line(a: a, b: b, c: c)

proc newLine*(p1, p2: Point): Line =
    result = newLine(p1.y - p2.y, p2.x - p1.x, p1.x * p2.y - p2.x * p1.y)

proc parallel_line*(self: Line, p: Point): Line =
    return newLine(self.a, self.b, -self.a * p.x - self.b * p.y)

proc perpendicular_line*(self: Line, p: Point): Line = newLine(p, p + newPoint(self.a, self.b))

proc greater_than(self, other: Line): bool =
        let (a, b) = (self.a, self.b)
        let (x, y) = (other.a, other.b)
        # b/a > y/x
        return b * x > a * y

proc `>`*(self, other: Line): bool = self.greater_than(other)

proc `<`*(self, other: Line): bool = other.greater_than(self)

proc same*(self, other: Line): bool =
    let (a, b, c) = (self.a, self.b, self.c)
    let (x, y, z) = (other.a, other.b, other.c)
    return close_enough(a * y, b * x) and close_enough(b * z, c * y)

proc equal*(self, other: Line): bool =
    let (a, b) = (self.a, self.b)
    let (x, y) = (other.a, other.b)
    # b/a == y/x
    return b * x == a * y

proc less_than*(self, other: Line): bool =
    let (a, b) = (self.a, self.b)
    let (x, y) = (other.a, other.b)
    # b/a > y/x
    return b * x < a * y

proc intersect*(self: Line; obj: Line): (Point, ) = (line_line_intersection(self, obj),)

proc intersect*(self: Line; obj: Circle): (Point, Point) = line_circle_intersection(self, obj)

proc at*(self: Line; x: Point): float = self.at(x.x, x.y)

proc at*(self: Line; x, y: float): float = x * self.a + y * self.b + self.c

proc is_parallel*(self, other: Line): bool = abs(self.a * other.b - self.b * other.a) < ATOM

proc is_perp*(self, other: Line): bool = abs(self.a * other.a + self.b * other.b) < ATOM

proc cross*(self, other: Line): float = self.a * other.b - self.b * other.a

proc dot*(self, other: Line): float = self.a * other.a + self.b * other.b

# Get a point on line closest to (x, y).
# ax + by + c = 0
proc point_at*(self: Line; x: Option[float], y: Option[float]): Point =
    let (a, b, c) = (self.a, self.b, self.c)
    if x.isNone and y.isSome():
        if a != 0:
            return newPoint((-c - b * y.get()) / a, y.get())
        else:
            return nil
    elif x.isSome() and y.isNone():
        if b != 0:
            return newPoint(x.get(), (-c - a * x.get()) / b)  # pylint: disable=invalid-unary-operand-type
        else:
            return nil
    else:
        if a * x.get() + b * y.get() + c == 0.0:
            return newPoint(x.get(), y.get())
        return nil

proc diff_side*(self: Line; p1, p2: Point): Option[bool] =
    let d1 = self.at(p1.x, p1.y)
    let d2 = self.at(p2.x, p2.y)
    if d1 == 0 or d2 == 0:
        return none(bool)
    return some(d1 * d2 < 0)

proc same_side*(self: Line; p1, p2: Point): Option[bool] =
    let d1 = self.at(p1.x, p1.y)
    let d2 = self.at(p2.x, p2.y)
    if d1 == 0 or d2 == 0:
        return none(bool)
    return some(d1 * d2 > 0)

proc sign*(self: Line, point: Point): int =
    let s = self.at(point.x, point.y)
    if s > 0:
        return 1
    elif s < 0:
        return -1
    return 0

proc is_same*(self, other: Line): bool =
    let (a, b, c) = (self.a, self.b, self.c)
    let (x, y, z) = (other.a, other.b, other.c)
    return abs(a * y - b * x) <= ATOM and abs(b * z - c * y) <= ATOM

proc close_enough*(a, b: float; tol: float = 1e-12): bool =
    return abs(a - b) < tol

# Sample a point within the boundary of points.
proc sample_within*(self: Line; points: openArray[Point]; n: int = 5): seq[Point] =
    var center = sum(points) * (1 / len(points))
    let radius = max(points.mapIt(it.distance(center)))
    if close_enough(center.distance(self), radius):
        center = center.foot(self)
    var (a, b) = line_circle_intersection(self, newCircle(center.foot(self), radius))

    if self of HalfLine:
        let h = HalfLine(self)
        if (a - h.tail).dot(h.head - h.tail) > 0:
            (a, b) = (h.tail, a)
        else:
            (a, b) = (h.tail, b)

    var r: Point = nil
    var best = -1.0
    var nn = n
    while nn > 0:
        dec(nn)
        let x = a + (b - a) * rand(1.0)
        let mind = min(points.mapIt(x.distance(it)))
        if mind > best:
            best = mind
            r = x

    return @[r]

type
    InvalidLineIntersectError* = object of ValueError

proc newHalfLine*(tail, head: Point): HalfLine =
    result = HalfLine(tail: tail, head: head)
    let r = newLine(tail, head)
    (result.a, result.b, result.c) = (r.a, r.b, r.c)

proc intersect*(self: HalfLine, obj: Line): Point = line_line_intersection(Line(self), obj)

proc intersect*(self: HalfLine, obj: Circle): Point =
    var exclude = @[self.tail]
    if obj of HoleCircle:
        exclude.add(HoleCircle(obj).hole)

    let (a, b) = line_circle_intersection(Line(self), obj)
    if exclude.mapIt(a.close(it)).anyIt(it):
        return b
    if exclude.mapIt(b.close(it)).anyIt(it):
        return a

    let v = self.head - self.tail
    let va = a - self.tail
    let vb = b - self.tail
    if v.dot(va) > 0:
        return a
    if v.dot(vb) > 0:
        return b
    raise newException(InvalidLineIntersectError, "")

proc perpendicular_bisector(p1, p2: Point): Line =
    let midpoint = (p1 + p2) * 0.5
    return newLine(midpoint, midpoint + newPoint(p2.y - p1.y, p1.x - p2.x))

proc same_sign*(a, b, c, d, e, f: Point): bool =
    let (ab, cb) = (a - b, c - b)
    let (de, fe) = (d - e, f - e)
    return (ab.x * cb.y - ab.y * cb.x) * (de.x * fe.y - de.y * fe.x) > 0

proc newCircle*(center: Point, radius: float): Circle =
    result = Circle(center: center, radius: radius)
    result.r2 = radius * radius

proc newCircle*(center: Point, p: Point): Circle =
    let radius = hypot(center.x - p.x, center.y - p.y)
    result = newCircle(center, radius)

proc newCircle*(p1, p2, p3: Point): Circle =
    let l12 = perpendicular_bisector(p1, p2)
    let l23 = perpendicular_bisector(p2, p3)
    let center = line_line_intersection(l12, l23)
    result = newCircle(center, p1)

proc intersect*(self: Circle, obj: Line): (Point, Point) = obj.intersect(self)

proc intersect*(self: Circle, obj: Circle): (Point, Point) = circle_circle_intersection(self, obj)

# Sample a point within the boundary of points.
proc sample_within*(self: Circle, points: openArray[Point], n: int = 5): seq[Point] =
    var r: Point = nil
    var best = -1.0
    var nn = n
    while nn > 0:
        dec(nn)
        let ang = rand(2.0) * math.PI
        let x = self.center + newPoint(cos(ang), sin(ang)) * self.radius
        let mind = min(points.mapIt(x.distance(it)))
        if mind > best:
            best = mind
            r = x

    return @[r]

proc newHoleCircle(center: Point, radius: float, hole: Point): HoleCircle =
    result = HoleCircle(center: center, radius: radius, hole: hole)
    result.r2 = radius * radius

method intersect*(self: HoleCircle; obj: Line): Point =
    let (a, b) = line_circle_intersection(obj, self)
    if a.close(self.hole):
        return b
    return a

method intersect*(self: HoleCircle; obj: HalfLine): Point =
    return obj.intersect(self)

method intersect*(self: HoleCircle; obj: Circle): Point =
    let (a, b) = circle_circle_intersection(obj, self)
    if a.close(self.hole):
        return b
    return a

method intersect*(self: HoleCircle; obj: HoleCircle): Point =
    let (a, b) = circle_circle_intersection(obj, self)
    if a.close(self.hole) or a.close(obj.hole):
        return b
    return a

# Solve a x^2 + bx + c = 0.
proc solve_quad*(a, b, c: float): Option[(float, float)] =
    var a2 = 2 * a
    let d = b * b - 2 * a2 * c
    if d < 0:
        return none((float, float))

    let y = math.sqrt(d)
    return some(((-b - y) / a2, (-b + y) / a2))

type
    InvalidQuadSolveError* = object of ValueError

# Returns a pair of Points as intersections of c1 and c2.
proc circle_circle_intersection*(c1, c2: Circle): (Point, Point) =
    # circle 1: (x0, y0), radius r0
    # circle 2: (x1, y1), radius r1
    let (x0, y0, r0) = (c1.center.x, c1.center.y, c1.radius)
    let (x1, y1, r1) = (c2.center.x, c2.center.y, c2.radius)

    let d = hypot(x1 - x0, y1 - y0)
    if d == 0:
        raise newException(InvalidQuadSolveError, "")

    let a = (r0 * r0 - r1 * r1 + d * d) / (2 * d)
    var h = r0 * r0 - a * a
    if h < 0:
        raise newException(InvalidQuadSolveError, "")

    h = sqrt(h)
    let x2 = x0 + a * (x1 - x0) / d
    let y2 = y0 + a * (y1 - y0) / d
    let x3 = x2 + h * (y1 - y0) / d
    let y3 = y2 - h * (x1 - x0) / d
    let x4 = x2 - h * (y1 - y0) / d
    let y4 = y2 + h * (x1 - x0) / d

    return (newPoint(x3, y3), newPoint(x4, y4))

# Returns a pair of points as intersections of line and circle.
proc line_circle_intersection*(line: Line, circle: Circle): (Point, Point) =
    let (a, b, c) = (line.a, line.b, line.c)
    let r = circle.radius
    let center = circle.center
    let (p, q) = (center.x, center.y)

    if b == 0:
        let x = -c / a
        let x_p = x - p
        let x_p2 = x_p * x_p
        let y = solve_quad(1, -2 * q, q * q + x_p2 - r * r)
        if y.isNone():
            raise newException(InvalidQuadSolveError, "")

        let (y1, y2) = y.get()
        return (newPoint(x, y1), newPoint(x, y2))

    if a == 0:
        let y = -c / b
        let y_q = y - q
        let y_q2 = y_q * y_q
        let x = solve_quad(1, -2 * p, p * p + y_q2 - r * r)
        if x.isNone():
            raise newException(InvalidQuadSolveError, "")
        let (x1, x2) = x.get()
        return (newPoint(x1, y), newPoint(x2, y))

    let c_ap = c + a * p
    let a2 = a * a
    let y = solve_quad(a2 + b * b, 2 * (b * c_ap - a2 * q), c_ap * c_ap + a2 * (q * q - r * r))
    if y.isNone():
        raise newException(InvalidQuadSolveError, "")
    let (y1, y2) = y.get()

    return (newPoint(-(b * y1 + c) / a, y1), newPoint(-(b * y2 + c) / a, y2))

# Whether a is between b & c.
proc check_between(a, b, c: Point): bool =
    return (a - b).dot(c - b) > 0 and (a - c).dot(b - c) > 0

proc circle_segment_intersect*(circle: Circle, p1, p2: Point): seq[Point] =
    let l = newLine(p1, p2)
    let (px, py) = line_circle_intersection(l, circle)

    result = @[]
    if check_between(px, p1, p2):
        result.add(px)
    if check_between(py, p1, p2):
        result.add(py)

proc line_segment_intersection*(l: Line, A, B: Point): Point =
    let (a, b, c) = (l.a, l.b, l.c)
    let (x1, y1, x2, y2) = (A.x, A.y, B.x, B.y)
    let (dx, dy) = (x2 - x1, y2 - y1)
    let d = a * dx + b * dy
    if d == 0:
        raise newException(InvalidLineIntersectError, "")
    let alpha = (-c - a * x1 - b * y1) / d
    return newPoint(x1 + alpha * dx, y1 + alpha * dy)

proc line_line_intersection*(l1, l2: Line): Point =
    let (a1, b1, c1) = (l1.a, l1.b, l1.c)
    let (a2, b2, c2) = (l2.a, l2.b, l2.c)
    # a1x + b1y + c1 = 0
    # a2x + b2y + c2 = 0
    let d = a1 * b2 - a2 * b1
    if d == 0:
        raise newException(InvalidLineIntersectError, "")
    return newPoint((c2 * b1 - c1 * b2) / d, (c1 * a2 - c2 * a1) / d)

proc check_too_close*(newpoints, points: openArray[Point], tol: float = 0.1): bool =
    if len(points) < 1:
        return false
    let avg = sum(points) * 1.0 / float(len(points))
    let mindist = min(points.mapIt(it.distance(avg)))
    for p0 in newpoints:
        for p1 in points:
            if p0.distance(p1) < tol * mindist:
                return true
    return false

proc check_too_far*(newpoints, points: openArray[Point]; tol: float = 4): bool =
    if len(points) < 2:
        return false
    let avg = sum(points) * 1.0 / float(len(points))
    let maxdist = max(points.mapIt(it.distance(avg)))
    for p in newpoints:
        if p.distance(avg) > maxdist * tol:
            return true
    return false

proc check_aconst*(args: seq[RootRef]): bool =
    var (ra, rb, rc, rd, rnum, rden) = seq2tuple6(args)
    var (a, b, c, d) = (Point(ra), Point(rb), Point(rc), Point(rd))
    let num = BoxedFloat(rnum).v
    let den = BoxedFloat(rden).v
    d = d + a - c
    var ang = ang_between(a, b, d)
    if ang < 0:
        ang += math.PI
    return close_enough(ang, num * math.PI / den)

proc check_circle*(points: openArray[Point]): bool =
    if len(points) != 4:
        return false
    let (o, a, b, c) = seq2tuple4(points)
    let (oa, ob, oc) = (o.distance(a), o.distance(b), o.distance(c))
    return close_enough(oa, ob) and close_enough(ob, oc)

proc check_coll*(points: openArray[Point]): bool =
    let (a, b) = (points[0], points[1])
    let l = newLine(a, b)
    for p in points[2..<len(points)]:
        if abs(l.at(p.x, p.y)) > ATOM:
            return false
    return true

proc check_ncoll*(points: openArray[Point]): bool =
    return not check_coll(points)

proc check_para*(points: openArray[Point]): bool =
    let (a, b, c, d) = seq2tuple4(points)
    let ab = newLine(a, b)
    let cd = newLine(c, d)
    if ab.same(cd):
        return false
    return ab.is_parallel(cd)

proc check_para_or_coll*(points: openArray[Point]): bool =
    return check_para(points) or check_coll(points)

proc check_perp*(points: openArray[Point]): bool =
    let (a, b, c, d) = seq2tuple4(points)
    let ab = newLine(a, b)
    let cd = newLine(c, d)
    return ab.is_perp(cd)

proc check_cyclic*(points: openArray[Point]): bool =
    let points_set = deduplicate(points)
    let (a , b, c) = seq2tuple3(points_set)
    let ps = points[3..<len(points_set)]
    if check_coll([a, b, c]):
        return false
    let circle = newCircle(a, b, c)
    for d in ps:
        if not close_enough(d.distance(circle.center), circle.radius):
            return false
    return true

proc check*(name: string; args: seq[Point]): Option[bool] =
    type check_fun = proc(args: openArray[Point]): bool
    var fun: check_fun = nil

    if name == "eqangle6":
        fun = check_eqangle
    elif name == "eqratio6":
        fun = check_eqratio
    elif name in ["simtri2", "simtri*"]:
        fun = check_simtri
    elif name in ["contri2", "contri*"]:
        fun = check_contri
    elif name == "para":
        fun = check_para_or_coll
    elif name == "on_line":
        fun = check_coll
    elif name in ["rcompute", "acompute"]:
        return some(true)
    elif name in ["fixl", "fixc", "fixb", "fixt", "fixp"]:
        return some(true)
    else:
        discard

    if fun == nil: return none(bool)

    return some(fun(args))

# Numerical check.
proc check*(name: string; args: seq[gm.Point[string]]): Option[bool] =
    result = check(name, args.mapIt(if it.num of Point: Point(it.num) else: raise newException(ObjectConversionDefect, "")))

proc check_sameside*(points: openArray[Point]): bool =
    let (b, a, c, y, x, z) = seq2tuple6(points)
    # whether b is to the same side of a & c as y is to x & z
    let ba = b - a
    let bc = b - c
    let yx = y - x
    let yz = y - z
    return ba.dot(bc) * yx.dot(yz) > 0

proc bring_together*(a, b, c, d: Point): (Point, Point, Point, Point) =
    let ab = newLine(a, b)
    let cd = newLine(c, d)
    let x = line_line_intersection(ab, cd)
    let unit = newCircle(x, 1.0)
    let (y, _) = line_circle_intersection(ab, unit)
    let (z, _) = line_circle_intersection(cd, unit)
    return (x, y, x, z)

proc same_clock*(a, b, c, d, e, f: Point): bool =
    let ba = b - a
    let cb = c - b
    let ed = e - d
    let fe = f - e
    return (ba.x * cb.y - ba.y * cb.x) * (ed.x * fe.y - ed.y * fe.x) > 0

# Check if the angle is equal to the given constant.
proc check_const_angle*(points: openArray[Point], values: openArray[float]): bool =
    var (a, b, c, d) = seq2tuple4(points)
    let (m, n) = seq2tuple2(values)
    (a, b, c, d) = bring_together(a, b, c, d)
    let ba = b - a
    let dc = d - c

    let a3 = math.arctan2(ba.y, ba.x)
    let a4 = math.arctan2(dc.y, dc.x)
    let y = a3 - a4

    return close_enough((m / n) mod 1, (y / math.PI) mod 1)

# Check if 8 points make 2 equal angles.
proc check_eqangle*(points: openArray[Point]): bool =
    var (a, b, c, d, e, f, g, h) = seq2tuple8(points)

    var ab = newLine(a, b)
    var cd = newLine(c, d)
    var ef = newLine(e, f)
    var gh = newLine(g, h)

    if ab.is_parallel(cd):
        return ef.is_parallel(gh)
    if ef.is_parallel(gh):
        return ab.is_parallel(cd)

    (a, b, c, d) = bring_together(a, b, c, d)
    (e, f, g, h) = bring_together(e, f, g, h)

    var ba = b - a
    let dc = d - c
    let fe = f - e
    let hg = h - g

    let sameclock = (ba.x * dc.y - ba.y * dc.x) * (fe.x * hg.y - fe.y * hg.x) > 0
    if not sameclock:
        ba = ba * -1.0

    let a1 = math.arctan2(fe.y, fe.x)
    let a2 = math.arctan2(hg.y, hg.x)
    let x = a1 - a2

    let a3 = math.arctan2(ba.y, ba.x)
    let a4 = math.arctan2(dc.y, dc.x)
    let y = a3 - a4

    let xy = (x - y) mod (2 * math.PI)
    return close_enough(xy, 0, tol=1e-11) or close_enough(xy, 2 * math.PI, tol=1e-11)

proc check_eqratio*(points: openArray[Point]): bool =
    let (a, b, c, d, e, f, g, h) = seq2tuple8(points)
    let ab = a.distance(b)
    let cd = c.distance(d)
    let ef = e.distance(f)
    let gh = g.distance(h)
    return close_enough(ab * gh, cd * ef)

proc check_cong*(points: openArray[Point]): bool =
    let (a, b, c, d) = seq2tuple4(points)
    return close_enough(a.distance(b), c.distance(d))

proc check_midp*(points: openArray[Point]): bool =
    let (a, b, c) = seq2tuple3(points)
    return check_coll(points) and close_enough(a.distance(b), a.distance(c))

# Check if 6 points make a pair of similar triangles.
proc check_simtri*(points: openArray[Point]): bool =
    let (a, b, c, x, y, z) = seq2tuple6(points)
    let ab = a.distance(b)
    let bc = b.distance(c)
    let ca = c.distance(a)
    let xy = x.distance(y)
    let yz = y.distance(z)
    let zx = z.distance(x)
    let tol = 1e-9
    return close_enough(ab * yz, bc * xy, tol) and close_enough(bc * zx, ca * yz, tol)

proc check_contri*(points: openArray[Point]): bool =
    let (a, b, c, x, y, z) = seq2tuple6(points)
    let ab = a.distance(b)
    let bc = b.distance(c)
    let ca = c.distance(a)
    let xy = x.distance(y)
    let yz = y.distance(z)
    let zx = z.distance(x)
    let tol = 1e-9
    return close_enough(ab, xy, tol) and close_enough(bc, yz, tol) and close_enough(ca, zx, tol)

proc check_ratio*(points: openArray[Point], values: openArray[float]): bool =
    let (a, b, c, d) = seq2tuple4(points)
    let (m, n) = seq2tuple2(values)
    let ab = a.distance(b)
    let cd = c.distance(d)
    return close_enough(ab * n, cd * m)

proc assert_close_enough*(a, b: float; tol: float = 1e-12) =
    assert(close_enough(a, b, tol), fmt"|{a}-{b}| = {abs(a-b)} >= {tol}")

proc ang_of*(tail, head: Point): float =
    let vector = head - tail
    return arctan2(vector.y, vector.x) mod (2 * math.PI)

proc ang_between*(tail: Point; head1, head2: Point): float =
    let ang1 = ang_of(tail, head1)
    let ang2 = ang_of(tail, head2)
    let diff = ang1 - ang2
    # return diff % (2*math.PI)
    if diff > math.PI:
        return diff - 2 * math.PI
    if diff < -math.PI:
        return 2 * math.PI + diff
    return diff

proc head_from*(tail: Point; ang: float; length: float = 1): Point =
    let vector = newPoint(math.cos(ang) * length, math.sin(ang) * length)
    return tail + vector

proc random_points*(n: int = 3): seq[Point] =
    return (1..n).mapIt(newPoint(rand(-1.0..1.0), rand(-1.0..1.0)))

func sum(points: openArray[Point]): Point =
    result = points[0]
    for i in 1 ..< len(points):
        result = result + points[i]

# Random rotate-flip-scale-shift a point cloud.
proc random_rfss*(points0: openArray[Point]): seq[Point] =
    # center point cloud.
    let average = sum(points0) * (1.0 / float(len(points0)))
    let points = points0.mapIt(it - average)

    # rotate
    let ang = rand(0.0..2 * math.PI)
    let (sin, cos) = (math.sin(ang), math.cos(ang))
    # scale and shift
    let scale = rand(0.5..2.0)
    let shift = newPoint(rand(-1.0..1.0), rand(-1.0..1.0))
    result = points.mapIt(it.rotate(sin, cos) * scale + shift)

    # randomly flip
    if rand(1.0) < 0.5:
        result = points.mapIt(it.flip())

proc random_rfss(a, b, c: Point): (Point, Point, Point) =
    let l = random_rfss([a, b, c])
    return (l[0], l[1], l[2])

proc random_rfss(a, b, c, d: Point): (Point, Point, Point, Point) =
    let l = random_rfss([a, b, c, d])
    return (l[0], l[1], l[2], l[3])

proc random_rfss(a, b, c, d, e: Point): (Point, Point, Point, Point, Point) =
    let l = random_rfss([a, b, d, e])
    return (l[0], l[1], l[2], l[3], l[4])

# Reduce intersecting objects into one point of intersections.
proc reduce*[T](objs: openArray[T], existing_points: openArray[Point]): seq[Point] =
    if len(objs) == 1:
        return objs[0].sample_within(existing_points)
    elif len(objs) == 2:
        let (a, b) = objs
        let r = a.intersect(b)
        if r of Point:
            return @[r]
        let a, b = r
        let a_close = any(existing_points.mapIt(a.close(it)))
        if a_close:
            return @[b]
        let b_close = any(existing_points.mapIt(b.close(it)))
        if b_close:
            return @[a]
        return @[random.sample([a, b])]
    else:
        raise ValueError(fmt"Cannot reduce {objs}")

proc sketch_on_opline*(args: openArray[Point]): HalfLine =
    let (a, b) = seq2tuple2(args)
    return newHalfLine(a, a + a - b)

proc sketch_on_hline*(args: openArray[Point]): HalfLine =
    let (a, b) = seq2tuple2(args)
    return newHalfLine(a, b)

proc sketch_ieq_triangle*(args: openArray[Point]): (Point, Point, Point) =
    let a = newPoint(0.0, 0.0)
    let b = newPoint(1.0, 0.0)
    let (c, _) = newCircle(a, b).intersect(newCircle(b, a))
    return (a, b, c)

proc sketch_incenter2*(args: openArray[Point]): (Point, Point, Point, Point) =
    let (a, b, c) = seq2tuple3(args)
    let l1 = sketch_bisect([b, a, c])
    let l2 = sketch_bisect([a, b, c])
    let i = line_line_intersection(l1, l2)
    let x = i.foot(newLine(b, c))
    let y = i.foot(newLine(c, a))
    let z = i.foot(newLine(a, b))
    return (x, y, z, i)

proc sketch_excenter2*(args: openArray[Point]): (Point, Point, Point, Point) =
    let (a, b, c) = seq2tuple3(args)
    let l1 = sketch_bisect([b, a, c])
    let l2 = sketch_exbisect([a, b, c])
    let i = line_line_intersection(l1, l2)
    let x = i.foot(newLine(b, c))
    let y = i.foot(newLine(c, a))
    let z = i.foot(newLine(a, b))
    return (x, y, z, i)

proc sketch_centroid*(args: openArray[Point]): (Point, Point, Point, Point) =
    let (a, b, c) = seq2tuple3(args)
    let x = (b + c) * 0.5
    let y = (c + a) * 0.5
    let z = (a + b) * 0.5
    let i = line_line_intersection(newLine(a, x), newLine(b, y))
    return (x, y, z, i)

proc sketch_ninepoints*(args: openArray[Point]): (Point, Point, Point, Point) =
    let (a, b, c) = seq2tuple3(args)
    let x = (b + c) * 0.5
    let y = (c + a) * 0.5
    let z = (a + b) * 0.5
    let cc = newCircle(x, y, z)
    return (x, y, z, cc.center)

# Sketch a circle touching two lines and another circle.
proc sketch_2l1c*(args: openArray[Point]): (Point, Point, Point, Point) =
    let (a, b, c, p) = seq2tuple4(args)
    let (bc, ac) = (newLine(b, c), newLine(a, c))
    let circle = newCircle(p, a)

    var (d, d2) = line_circle_intersection(p.perpendicular_line(bc), circle)
    if bc.diff_side(d2, a).get(false):
        d = d2

    var (e, e2) = line_circle_intersection(p.perpendicular_line(ac), circle)
    if ac.diff_side(e2, b).get(false):
        e = e2

    let df = d.perpendicular_line(newLine(p, d))
    let ef = e.perpendicular_line(newLine(p, e))
    let f = line_line_intersection(df, ef)

    var (g, g2) = line_circle_intersection(newLine(c, f), circle)
    if bc.same_side(g2, a).get(false):
        g = g2

    let b2 = c + (b - c) / b.distance(c)
    let a2 = c + (a - c) / a.distance(c)
    let m = (a2 + b2) * 0.5
    let x = line_line_intersection(newLine(c, m), newLine(p, g))
    return (x.foot(ac), x.foot(bc), g, x)

proc sketch_3peq*(args: openArray[Point]): (Point, Point, Point) =
    let (a, b, c) = seq2tuple3(args)
    let (ab, bc, ca) = (newLine(a, b), newLine(b, c), newLine(c, a))

    let z = b + (c - b) * rand(-0.5..1.5)

    let z2 = z * 2 - c
    let l = z2.parallel_line(ca)
    let x = line_line_intersection(l, ab)
    let y = z * 2 - x
    return (x, y, z)

proc sketch*[T, DepsT](name: string; args: seq[gm.Point[DepsT]]): seq[T] =
    return sketch[T](name, args.mapIt(it.num))

proc sketch*[T](name: string; args: openArray[Point]): seq[T] =
    # TODO
    return @[]

# Try to sketch an intersection between two objects.
proc try_to_sketch_intersect*(name1: string, args1: seq[gm.Point or Point],
    name2: string, args2: seq[gm.Point or Point], existing_points: seq[Point]): Option[Point] =
    var obj1 = sketch(name1, args1)[0]
    var obj2 = sketch(name2, args2)[0]

    func pp(x1, x2: Point): Option[Point] =
        let close1 = check_too_close([x1], existing_points)
        let far1 = check_too_far([x1], existing_points)
        if not close1 and not far1:
            return x1
        let close2 = check_too_close([x2], existing_points)
        let far2 = check_too_far([x2], existing_points)
        if not close2 and not far2:
            return x2
        return none(Point)

    if (obj1 of Line) and (obj2 of Line):
        try:
            return line_line_intersection(Line(obj1), Line(obj2))
        except:
            return none(Point)
    elif (obj1 of Circle) and (obj2 of Circle):
        try:
            let (x1, x2) = circle_circle_intersection(Circle(obj1), Circle(obj2))
            return pp(x1, x2)
        except:
            return none(Point)
    else:
        if (obj2 of Line) and (obj1 of Circle):
            (obj1, obj2) = (obj2, obj1)

        try:
            let (x1, x2) = line_circle_intersection(Line(obj1), Circle(obj2))
            return pp(x1, x2)
        except:
            return none(Point)

proc sketch_acircle*(args: openArray[Point]): Circle =
    let (a, b, c, d, f) = seq2tuple5(args)
    let de = sketch_aline([c, a, b, f, d])
    let fe = sketch_aline([a, c, b, d, f])
    let e = line_line_intersection(de, fe)
    return newCircle(d, e, f)

# Sketch the construction aline.
proc sketch_aline*(args: openArray[Point]): HalfLine =
    let (A, B, C, D, E) = seq2tuple5(args)
    let ab = A - B
    let cb = C - B
    let de = D - E

    let dab = A.distance(B)
    let ang_ab = math.arctan2(ab.y / dab, ab.x / dab)

    let dcb = C.distance(B)
    let ang_bc = math.arctan2(cb.y / dcb, cb.x / dcb)

    let dde = D.distance(E)
    let ang_de = math.arctan2(de.y / dde, de.x / dde)

    let ang_ex = ang_de + ang_bc - ang_ab
    let X = E + newPoint(math.cos(ang_ex), math.sin(ang_ex))
    return newHalfLine(E, X)

# Sketch the angle mirror.
proc sketch_amirror*(args: openArray[Point]): HalfLine =
    let (A, B, C) = seq2tuple3(args)
    let ab = A - B
    let cb = C - B

    let dab = A.distance(B)
    let ang_ab = math.arctan2(ab.y / dab, ab.x / dab)
    let dcb = C.distance(B)
    let ang_bc = math.arctan2(cb.y / dcb, cb.x / dcb)

    let ang_bx = 2 * ang_bc - ang_ab
    let X = B + newPoint(math.cos(ang_bx), math.sin(ang_bx))
    return newHalfLine(B, X)

proc sketch_bisect*(args: openArray[Point]): Line =
    let (a, b, c) = seq2tuple3(args)
    let ab = a.distance(b)
    let bc = b.distance(c)
    let x = b + (c - b) * (ab / bc)
    let m = (a + x) * 0.5
    return newLine(b, m)

proc sketch_exbisect*(args: openArray[Point]): Line =
    let (a, b, c) = seq2tuple3(args)
    return sketch_bisect(args).perpendicular_line(b)

proc sketch_bline*(args: openArray[Point]): Line =
    let (a, b) = seq2tuple2(args)
    let m = (a + b) * 0.5
    return m.perpendicular_line(newLine(a, b))

proc sketch_dia*(args: openArray[Point]): Circle =
    let (a, b) = seq2tuple2(args)
    return newCircle((a + b) * 0.5, a)

proc sketch_tangent*(args: openArray[Point]): (Point, Point) =
    let (a, o, b) = seq2tuple3(args)
    let dia = sketch_dia([a, o])
    return circle_circle_intersection(newCircle(o, b), dia)

proc sketch_circle*(args: openArray[Point]): Circle =
    let (a, b, c) = seq2tuple3(args)
    return newCircle(a, b.distance(c))

# Sketch tangents to two circles.
proc sketch_cc_tangent*(args: openArray[Point]): (Point, Point, Point, Point) =
    var (o, a, w, b) = seq2tuple4(args)
    var (ra, rb) = (o.distance(a), w.distance(b))

    var ow = newLine(o, w)
    if close_enough(ra, rb):
        let oo = ow.perpendicular_line(o)
        let oa = newCircle(o, ra)
        let (x, z) = line_circle_intersection(oo, oa)
        let y = x + w - o
        let t = z + w - o
        return (x, y, z, t)

    let swap = rb > ra
    if swap:
        (o, a, w, b) = (w, b, o, a)
        (ra, rb) = (rb, ra)

    let oa = newCircle(o, ra)
    let q = o + (w - o) * ra / (ra - rb)

    var (x, z) = circle_circle_intersection(sketch_dia([o, q]), oa)
    var y = w.foot(newLine(x, q))
    var t = w.foot(newLine(z, q))

    if swap:
        (x, y, z, t) = (y, x, t, z)

    return (x, y, z, t)

proc sketch_hcircle*(args: openArray[Point]): HoleCircle =
    let (a, b) = seq2tuple2(args)
    return newHoleCircle(center=a, radius=a.distance(b), hole=b)

proc sketch_e5128*(args: openArray[Point]): (Point, Point) =
    let (a, b, c, d) = seq2tuple4(args)
    let ad = newLine(a, d)

    let g = (a + b) * 0.5
    let de = newLine(d, g)

    var (e, f) = line_circle_intersection(de, newCircle(c, b))

    if e.distance(d) < f.distance(d):
        e = f
    return (e, g)

# Sketch quadrangle with two equal opposite sides.
proc sketch_eq_quadrangle*(args: openArray[Point]): (Point, Point, Point, Point) =
    var a = newPoint(0.0, 0.0)
    var b = newPoint(1.0, 0.0)

    let length = rand(0.5..2.0)
    var ang = rand(math.PI / 3 .. math.PI * 2 / 3)
    var d = head_from(a, ang, length)

    ang = ang_of(b, d)
    ang = rand(ang / 10 .. ang / 9)
    var c = head_from(b, ang, length)
    (a, b, c, d) = seq2tuple4(random_rfss([a, b, c, d]))
    return (a, b, c, d)

proc sketch_eq_trapezoid*(args: openArray[Point]): (Point, Point, Point, Point) =
    var a = newPoint(0.0, 0.0)
    var b = newPoint(1.0, 0.0)
    let l = rand(0.5 .. 2.0)

    let height = rand(0.5 .. 2.0)
    var c = newPoint(0.5 + l / 2.0, height)
    var d = newPoint(0.5 - l / 2.0, height)

    (a, b, c, d) = random_rfss(a, b, c, d)
    return (a, b, c, d)

# Sketch the def eqangle2.
proc sketch_eqangle2*(args: openArray[Point]): Point =
    let (a, b, c) = seq2tuple3(args)

    let d = c * 2 - b

    let ba = b.distance(a)
    let bc = b.distance(c)
    let l = ba * ba / bc
    var be = 0.0

    if rand(0.0 .. 1.0) < 0.5:
        be = min(l, bc)
        be = rand(be * 0.1 .. be * 0.9)
    else:
        be = max(l, bc)
        be = rand(be * 1.1 .. be * 1.5)

    let e = b + (c - b) * (be / bc)
    let y = b + (a - b) * (be / l)
    return line_line_intersection(newLine(c, y), newLine(a, e))

proc sketch_eqangle3*(args: openArray[Point]): Circle =
    let (a, b, d, e, f) = seq2tuple5(args)
    let de = d.distance(e)
    let ef = e.distance(f)
    let ab = b.distance(a)
    let ang_ax = ang_of(a, b) + ang_between(e, d, f)
    let x = head_from(a, ang_ax, length=de / ef * ab)
    return newCircle(a, b, x)

# Sketch quadrangle with two equal diagonals.
proc sketch_eqdia_quadrangle*(args: openArray[Point]): (Point, Point, Point, Point) =
    let m = rand(0.3 .. 0.7)
    let n = rand(0.3 .. 0.7)
    var a = newPoint(-m, 0.0)
    var c = newPoint(1 - m, 0.0)
    var b = newPoint(0.0, -n)
    var d = newPoint(0.0, 1 - n)

    let ang = rand(-0.25 * math.PI .. 0.25 * math.PI)
    let (sin, cos) = (math.sin(ang), math.cos(ang))
    b = b.rotate(sin, cos)
    d = d.rotate(sin, cos)
    (a, b, c, d) = random_rfss(a, b, c, d)
    return (a, b, c, d)

proc sketch_free*(args: openArray[Point]): Point =
    return random_points(1)[0]

proc sketch_isos*(args: openArray[Point]): (Point, Point, Point) =
    let base = rand(0.5 .. 1.5)
    let height = rand(0.5 .. 1.5)

    var b = newPoint(-base / 2, 0.0)
    var c = newPoint(base / 2, 0.0)
    var a = newPoint(0.0, height)
    (a, b, c) = random_rfss(a, b, c)
    return (a, b, c)

proc sketch_line*(args: openArray[Point]): Line =
    let (a, b) = seq2tuple2(args)
    return newLine(a, b)

proc sketch_cyclic*(args: openArray[Point]): Circle =
    let (a, b, c) = seq2tuple3(args)
    return newCircle(a, b, c)

proc sketch_hline*(args: openArray[Point]): HalfLine =
    let (a, b) = seq2tuple2(args)
    return newHalfLine(a, b)

proc sketch_midp*(args: openArray[Point]): Point =
    let (a, b) = seq2tuple2(args)
    return (a + b) * 0.5

proc sketch_pentagon*(args: openArray[Point]): (Point, Point, Point, Point, Point) =
    var points = @[newPoint(1.0, 0.0)]
    var ang = 0.0

    for i in 1..4:
        ang += (2 * math.PI - ang) / (5.0 - float(i)) * rand(0.5 .. 1.5)
        let point = newPoint(math.cos(ang), math.sin(ang))
        points.add(point)

    var (a, b, c, d, e) = seq2tuple5(points)
    (a, b, c, d, e) = random_rfss(a, b, c, d, e)
    return (a, b, c, d, e)

proc sketch_pline*(args: openArray[Point]): Line =
    let (a, b, c) = seq2tuple3(args)
    return a.parallel_line(newLine(b, c))

proc sketch_pmirror*(args: openArray[Point]): Point =
    let (a, b) = seq2tuple2(args)
    return b * 2 - a

# Sketch a random quadrangle.
proc sketch_quadrangle*(args: openArray[Point]): (Point, Point, Point, Point) =
    let m = rand(0.3 .. 0.7)
    let n = rand(0.3 .. 0.7)

    var a = newPoint(-m, 0.0)
    var c = newPoint(1 - m, 0.0)
    var b = newPoint(0.0, -rand(0.25 .. 0.75))
    var d = newPoint(0.0, rand(0.25 .. 0.75))

    let ang = rand(-0.25 * math.PI .. 0.25 * math.PI)
    let (sin, cos) = (math.sin(ang), math.cos(ang))
    b = b.rotate(sin, cos)
    d = d.rotate(sin, cos)
    (a, b, c, d) = random_rfss(a, b, c, d)
    return (a, b, c, d)

proc sketch_r_trapezoid*(args: openArray[Point]): (Point, Point, Point, Point) =
    var a = newPoint(0.0, 1.0)
    var d = newPoint(0.0, 0.0)
    var b = newPoint(rand(0.5 .. 1.5), 1.0)
    var c = newPoint(rand(0.5 .. 1.5), 0.0)
    (a, b, c, d) = random_rfss(a, b, c, d)
    return (a, b, c, d)

proc sketch_r_triangle*(args: openArray[Point]): (Point, Point, Point) =
    var a = newPoint(0.0, 0.0)
    var b = newPoint(0.0, rand(0.5 .. 2.0))
    var c = newPoint(rand(0.5 .. 2.0), 0.0)
    (a, b, c) = random_rfss(a, b, c)
    return (a, b, c)

proc sketch_rectangle*(args: openArray[Point]): (Point, Point, Point, Point) =
    var a = newPoint(0.0, 0.0)
    var b = newPoint(0.0, 1.0)
    let l = rand(0.5 .. 2.0)
    var c = newPoint(l, 1.0)
    var d = newPoint(l, 0.0)
    (a, b, c, d) = random_rfss(a, b, c, d)
    return (a, b, c, d)

proc sketch_reflect*(args: openArray[Point]): Point =
    let (a, b, c) = seq2tuple3(args)
    let m = a.foot(newLine(b, c))
    return m * 2 - a

proc sketch_risos*(args: openArray[Point]): (Point, Point, Point) =
    var a = newPoint(0.0, 0.0)
    var b = newPoint(0.0, 1.0)
    var c = newPoint(1.0, 0.0)
    (a, b, c) = random_rfss(a, b, c)
    return (a, b, c)

proc sketch_rotaten90*(args: openArray[Point]): Point =
    let (a, b) = seq2tuple2(args)
    let ang = -math.PI / 2
    return a + (b - a).rotate(math.sin(ang), math.cos(ang))

proc sketch_rotatep90*(args: openArray[Point]): Point =
    let (a, b) = seq2tuple2(args)
    let ang = math.PI / 2
    return a + (b - a).rotate(math.sin(ang), math.cos(ang))

proc sketch_s_angle*(args: openArray[Point], values: openArray[float]): HalfLine =
    let (a, b) = seq2tuple2(args)
    let y = values[0]
    let ang = y / 180 * math.PI
    let x = b + (a - b).rotate(ang)
    return newHalfLine(b, x)

proc sketch_segment*(args: openArray[Point]): (Point, Point) =
    let (a, b) = seq2tuple2(random_points(2))
    return (a, b)

proc sketch_shift*(args: openArray[Point]): Point =
    let (a, b, c) = seq2tuple3(args)
    return c + (b - a)

proc sketch_square*(args: openArray[Point]): (Point, Point) =
    let (a, b) = seq2tuple2(args)
    let c = b + (a - b).rotate(-math.PI / 2)
    let d = a + (b - a).rotate(math.PI / 2)
    return (c, d)

proc sketch_isquare*(args: openArray[Point]): (Point, Point, Point, Point) =
    var a = newPoint(0.0, 0.0)
    var b = newPoint(1.0, 0.0)
    var c = newPoint(1.0, 1.0)
    var d = newPoint(0.0, 1.0)
    (a, b, c, d) = random_rfss(a, b, c, d)
    return (a, b, c, d)

proc sketch_tline*(args: openArray[Point]): Line =
    let (a, b, c) = seq2tuple3(args)
    return a.perpendicular_line(newLine(b, c))

proc sketch_trapezoid*(args: openArray[Point]): (Point, Point, Point, Point) =
    var d = newPoint(0.0, 0.0)
    var c = newPoint(1.0, 0.0)

    let base = rand(0.5 .. 2.0)
    let height = rand(0.5 .. 2.0)
    var a = newPoint(rand(0.2 .. 0.5), height)
    var b = newPoint(a.x + base, height)
    (a, b, c, d) = random_rfss(a, b, c, d)
    return (a, b, c, d)

proc sketch_triangle*(args: openArray[Point]): (Point, Point, Point) =
    let a = newPoint(0.0, 0.0)
    let b = newPoint(1.0, 0.0)
    let ac = rand(0.5 .. 2.0)
    let ang = rand(0.2 .. 0.8) * math.PI
    let c = head_from(a, ang, ac)
    return (a, b, c)

proc sketch_triangle12*(args: openArray[Point]): (Point, Point, Point) =
    var b = newPoint(0.0, 0.0)
    var c = newPoint(rand(1.5 .. 2.5), 0.0)
    var (a, _) = circle_circle_intersection(newCircle(b, 1.0), newCircle(c, 2.0))
    (a, b, c) = random_rfss(a, b, c)
    return (a, b, c)

# Sketch two trisectors of an angle.
proc sketch_trisect*(args: openArray[Point]): (Point, Point) =
    let (a, b, c) = seq2tuple3(args)
    var ang1 = ang_of(b, a)
    var ang2 = ang_of(b, c)

    var swap = 0
    if ang1 > ang2:
        (ang1, ang2) = (ang2, ang1)
        swap += 1

    if ang2 - ang1 > math.PI:
        (ang1, ang2) = (ang2, ang1 + 2 * math.PI)
        swap += 1

    let angx = ang1 + (ang2 - ang1) / 3
    let angy = ang2 - (ang2 - ang1) / 3

    var x = b + newPoint(math.cos(angx), math.sin(angx))
    var y = b + newPoint(math.cos(angy), math.sin(angy))

    let ac = newLine(a, c)
    x = line_line_intersection(newLine(b, x), ac)
    y = line_line_intersection(newLine(b, y), ac)

    if swap == 1:
        return (y, x)
    return (x, y)

proc sketch_trisegment*(args: openArray[Point]): (Point, Point) =
    let (a, b) = seq2tuple2(args)
    let (x, y) = (a + (b - a) * (1.0 / 3), a + (b - a) * (2.0 / 3))
    return (x, y)
