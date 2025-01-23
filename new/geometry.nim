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

import std/tables
from std/sequtils import zip
import sets
import std/enumerate
import std/[sequtils, sugar, strformat, hashes]

type
    Node*[DepsT] = ref object of RootObj
        name: string
        edge_graph: Table[Node[DepsT], Table[Node[DepsT], seq[DepsT]]]
        merge_graph: Table[Node[DepsT], DepsT]
        rep_by: Node[DepsT]
        members: HashSet[Node[DepsT]]

        val: Node[DepsT]
        obj: Node[DepsT]

        num*: RootRef
        change: HashSet[Node[DepsT]]

        new_name: string = ""

    Point*[DepsT] = ref object of Node[DepsT]

    Line*[DepsT] = ref object of Node[DepsT]

    Segment*[DepsT] = ref object of Node[DepsT]

    Circle*[DepsT] = ref object of Node[DepsT]

    Direction*[DepsT] = ref object of Node[DepsT]

    Angle*[DepsT] = ref object of Node[DepsT]
        d: tuple[n1: Direction[DepsT], n2: Direction[DepsT]]

    Measure*[DepsT] = ref object of Node[DepsT]

    Length*[DepsT] = ref object of Node[DepsT]
        d: tuple[n1: Node[DepsT], n2: Node[DepsT]]

    Ratio*[DepsT] = ref object of Node[DepsT]
        l: (Length[DepsT], Length[DepsT])

    Value*[DepsT] = ref object of Node[DepsT]

proc hash*[DepsT](node: Node[DepsT]): Hash =
    return hash(node.name)

proc `$`*[DepsT](node: Node[DepsT]): string =
    return "Node(" & ($ cast[uint](node)) & ": " & node.name & ")"

proc initNode[DepsT](node: var Node[DepsT]; name: string) =
    node.name = name
    node.edge_graph = initTable[Node[DepsT], Table[Node[DepsT], seq[DepsT]]]()
    node.merge_graph = initTable[Node[DepsT], DepsT]()
    node.rep_by = nil
    node.members = initHashSet[Node[DepsT]]()
    node.val = nil
    node.obj = nil
    node.num = nil
    node.change = initHashSet[Node[DepsT]]()
    node.members.incl(node)

proc newNodeType*[NodeType](name: string = ""): NodeType =
    result = NodeType()
    initNode(result, name)

proc newNode*[DepsT](name: string = ""): Node[DepsT] =
    result = newNodeType[Node[DepsT]](name)

proc update*[A, B](self: var Table[A, B], target: Table[A, B]) =
    for k, v in target.pairs:
        self[k] = v

proc merge_edge_graph*[DepsT](self: Node[DepsT], new_edge_graph: Table[Node[DepsT], Table[Node[DepsT], seq[DepsT]]])

proc set_rep*[DepsT](self: Node[DepsT], node: Node[DepsT]) =
    if self == node:
        return
    self.rep_by = node
    node.merge_edge_graph(self.edge_graph)
    for x in self.members.items():
        node.members.incl(x)

proc rep*[DepsT](self: Node[DepsT]): Node[DepsT] =
    var x = self
    while x.rep_by != nil:
        x = x.rep_by
    return x

proc why_rep*[DepsT](self: Node[DepsT]): seq[Node[DepsT]] =
    return self.why_equal([self.rep()], nil)

type RepWhy[DepsT] = tuple[rep: Node[DepsT], why: seq[Node]]

proc rep_and_why*[DepsT](self: Node[DepsT]): RepWhy =
    let a_rep = self.rep()
    return (rep: a_rep, why: self.why_equal([rep], nil))

# Neighbors of this node in the proof state graph.
proc neighbors_set*[DepsT](self: Node[DepsT]; oftype: typedesc; do_rep: bool = true): HashSet[Node[DepsT]] =
    let rep = if do_rep: self.rep() else: self

    result = initHashSet[Node]()

    for n in rep.edge_graph.keys():
        if oftype is nil or (n of oftype):
            if do_rep:
                result.add(n.rep())
            else:
                result.add(n)

    return result

proc neighbors*[DepsT](self: Node[DepsT]; oftype: typedesc; do_rep: bool = true): seq[Node[DepsT]] =
    let s = self.neighbors_set(oftype, do_rep)
    return s.toSeq()

proc merge_edge_graph*[DepsT](self: Node[DepsT], new_edge_graph: Table[Node[DepsT], Table[Node[DepsT], seq[DepsT]]]) =
    for x, xdict in new_edge_graph.pairs:
        if x in self.edge_graph:
            self.edge_graph[x].update(xdict)
        else:
            self.edge_graph[x] = xdict

proc merge_one*[DepsT, NodeT](self, node: NodeT, deps: DepsT) =
    node.rep().set_rep(self.rep())

    if node in self.merge_graph:
        return

    self.merge_graph[node] = deps
    node.merge_graph[self] = deps

proc merge*[DepsT, NodeT](self: NodeT, nodes: openArray[NodeT], deps: DepsT) =
    for i in 0..<len(nodes):
        self.merge_one(nodes[i], deps)

proc is_val*[DepsT](self, node: Node[DepsT]): bool =
    return (self of Line[DepsT]) and (node of Direction[DepsT]) or
           (self of Segment[DepsT]) and (node of Length[DepsT]) or
           (self of Angle[DepsT]) and (node of Measure[DepsT]) or
           (self of Ratio[DepsT]) and (node of Value[DepsT])

proc set_val*[DepsT](self, node: Node[DepsT]) =
    self.val = node

proc set_obj*[DepsT](self, node: Node[DepsT]) =
    self.obj = node

proc val*[DepsT](self: Node[DepsT]): Node[DepsT] =
    return if self.val == nil: nil else: self.val.rep()

proc obj*[DepsT](self: Node[DepsT]): Node[DepsT] =
    return if self.obj == nil: nil else: self.obj.rep()

proc equivs*(self: Node): HashSet[Node] =
    return self.rep().members

proc connect_to*[DepsT](self: var Node[DepsT]; node: Node[DepsT]; deps: seq[DepsT] = @[]) =
    var rep = self.rep()

    if node in rep.edge_graph:
        rep.edge_graph[node].update({self: deps}.toTable())
    else:
        rep.edge_graph[node] = {self: deps}.toTable()

    if self.is_val(node):
        self.set_val(node)
        node.set_obj(self)

# What are the equivalent nodes up to a certain level.
proc equivs_upto*[DepsT](self: Node[DepsT], level: int = -1): Table[Node[DepsT], Node[DepsT]] =
    result = initTable[Node[DepsT], Node[DepsT]]()
    result[self] = nil

    var visited = initHashSet[Node[DepsT]]()
    var queue: seq[Node[DepsT]] = @[self]
    var i = 0

    while i < len(queue):
        let current = queue[i]
        i += 1
        visited.incl(current)

        for neighbor in current.merge_graph.keys():
            if (level >= 0) and (current.merge_graph[neighbor].level >= level):
                continue
            if not (neighbor in visited):
                queue.append(neighbor)
                result[neighbor] = current

    return result

proc bfs_backtrack[DepsT](root: Node[DepsT], leafs: seq[Node[DepsT]], parent: Table[Node[DepsT], Node[DepsT]]): seq[DepsT]

# BFS why this node is equal to other nodes.

proc why_equal*[DepsT](self: Node[DepsT]; others: seq[Node[DepsT]]): seq[DepsT] =
    var others0 = others.toHashSet()
    var found = 0

    var parent = initTable[Node[DepsT], Node[DepsT]]()
    var queue: seq[Node[DepsT]] = @[self]
    var i = 0

    while i < len(queue):
        let current = queue[i]
        if current in others0:
            found += 1
        if found == len(others0):
            break

        i += 1

        for neighbor in current.merge_graph.keys():
            if not (neighbor in parent):
                queue.add(neighbor)
                parent[neighbor] = current

    return bfs_backtrack(self, others, parent)

proc why_equal*[DepsT](self: Node[DepsT]; others: seq[Node[DepsT]]; level: int): seq[DepsT] =
    var others0 = others.toHashSet()
    var found = 0

    var parent = initTable[Node[DepsT], Node[DepsT]]()
    var queue: seq[Node[DepsT]] = @[self]
    var i = 0

    while i < len(queue):
        let current = queue[i]
        if current in others0:
            found += 1
        if found == len(others0):
            break

        i += 1

        for neighbor in current.merge_graph.keys():
            if (level > 0) and (current.merge_graph[neighbor].level >= level):
                continue
            if not (neighbor in parent):
                queue.append(neighbor)
                parent[neighbor] = current

    return bfs_backtrack(self, others, parent)

# BFS for why self is equal to at least one member of each group.
proc why_equal_groups*[DepsT](self: Node[DepsT], groups: seq[seq[Node[DepsT]]], level: int): tuple[a: seq[Node[DepsT]], b: seq[Node[DepsT]]] =
    var others = newSeq[Node[DepsT]](groups.len)
    var found = 0

    var parent = initTable[Node[DepsT], Node[DepsT]]()
    var queue: seq[Node[DepsT]] = @[self]
    var i = 0

    while i < len(queue):
        let current = queue[i]

        for (j, grp) in enumerate(groups):
            if (others[j] is nil) and (current in grp):
                others[j] = current
                found += 1

        if found == len(others):
            break

        i += 1

        for neighbor in current.merge_graph.keys():
            if (level >= 0) and (current.merge_graph[neighbor].level >= level):
                continue
            if not (neighbor in parent):
                queue.append(neighbor)
                parent[neighbor] = current

    return (bfs_backtrack(self, others, parent), others)

proc why_val*[DepsT](self: Node[DepsT], level: int): seq[DepsT] =
    if level >= 0:
        return self.val.why_equal(@[self.val], level)
    else:
        return self.val.why_equal(@[self.val])

proc why_connect*[DepsT](self, node: Node[DepsT]; level: int = -1): seq[DepsT] =
    let rep = self.rep()
    let equivs = rep.edge_graph[node].keys().toSeq()
    if len(equivs) < 1:
        return nil
    let equiv = equivs[0]
    let dep = rep.edge_graph[node][equiv]
    return [dep] + self.why_equal(equiv, level)

proc why_connect*[DepsT](pairs: seq[tuple[a: Node[DepsT], b: Node[DepsT]]]): seq[DepsT] =
    result = @[]
    for node1, node2 in pairs:
        result += node1.why_connect(node2)
    return result

proc is_equiv*[DepsT](x, y: Node[DepsT]; level: int = -1): bool =
    return x.why_equal([y], level) != nil

proc is_equal*[DepsT](x, y: Node[DepsT], level: int = -1): bool =
    if x == y:
        return true
    if (x.val == nil) or (y.val == nil):
        return false
    if x.val != y.val:
        return false
    return is_equiv(x.val, y.val, level)

# Return the path given BFS trace of parent nodes.
proc bfs_backtrack[DepsT](root: Node[DepsT], leafs: seq[Node[DepsT]], parent: Table[Node[DepsT], Node[DepsT]]): seq[DepsT] =
    var backtracked = [root].toHashSet()
    result = @[]

    for node in leafs:
        if node == nil:
            return @[]

        if node in backtracked:
            continue
        if not (node in parent):
            return @[]
        var nn = node
        while not (nn in backtracked):
            backtracked.incl(nn)
            result.add(nn.merge_graph[parent[nn]])
            nn = parent[nn]

    return result

proc new_val*[DepsT](self: Line[DepsT]): Direction[DepsT] =
    return newNodeType[Direction[DepsT]]()

# Why points are connected to self.
proc why_connected_to_self[DepsT](self: Node[DepsT]; points: seq[Node[DepsT]]; level: int = -1): seq[DepsT] =
    var groups: seq[seq[Node[DepsT]]] = @[]
    for p in points:
        var group: seq[Node[DepsT]] = @[]
        for l, d in self.edge_graph[p].pairs:
            if (d == nil) or (d.level < level):
                group.append(l)

        if len(group) < 1:
            return nil

        groups.append(group)

    var min_deps: seq[DepsT] = @[]
    for line in groups[0]:
        let (deps, others) = line.why_equal_groups(groups[1..<len(groups)], level)
        if deps == @[]:
            continue
        for p, o in zip(points, [line] + others):
            deps.append(self.edge_graph[p][o])
        if (min_deps == @[]) or (len(deps) < len(min_deps)):
            min_deps = deps

    if min_deps == @[]:
        return @[]

    return min_deps.filterIt(len(it) > 0)

# Why points are connected to self.
proc why_coll*[DepsT](self: Line[DepsT]; points: seq[Point[DepsT]]; level: int = -1): seq[DepsT] =
    return why_connected_to_self(self, points, level)

proc new_val*[DepsT](self: Segment[DepsT]): Length[DepsT] =
    return newNodeType[Length[DepsT]]()

# Why points are connected to self.
proc why_cyclic*[DepsT](self: Circle[DepsT]; points: seq[Point], level: int = -1): seq[DepsT] =
    return why_connected_to_self(self, points, level)

proc why_equal*[DepsT](x: Node[DepsT]; y: Node[DepsT]; level: int): seq[DepsT] =
    if x == y:
        return @[]
    if (x.val == nil) or (y.val == nil):
        return @[]
    if x.val == y.val:
        return @[]
    return x.val.why_equal(@[y.val], level)

proc why_equal*[DepsT](x: Node[DepsT]; y: Node[DepsT]): seq[DepsT] =
    if x == y:
        return @[]
    if (x.val == nil) or (y.val == nil):
        return @[]
    if x.val == y.val:
        return @[]
    return x.val.why_equal(@[y.val])

proc get_node_type_thru_all[DepsT](points: seq[Point[DepsT]], T: typedesc): seq[T] =
    let points_num = len(points.toHashSet())
    var line2count = initCountTable[Node[DepsT]]()
    for p in points:
        for l in p.neighbors(T):
            line2count[l] += 1
    return toSeq(line2count.pairs()).filterIt(len(it[1]) == points_num).mapIt(it[0])

proc get_lines_thru_all*[DepsT](points: seq[Point[DepsT]]): seq[Line[DepsT]] =
    return get_node_type_thru_all(points, Line[DepsT])

# Why points are collinear.
proc line_of_and_why*[DepsT](points: seq[Point[DepsT]]; level: int = -1): tuple[a: Line, b: seq[DepsT]] =
    for l0 in get_lines_thru_all(points):
        for l in l0.equivs():
            if all(points, (p) => p in l.edge_graph):
                let (x, y) = l.points
                var points_set = points.toHashSet()
                points_set.incl(x)
                points_set.incl(y)
                let colls = points_set.toSeq()
                # if len(colls) < 3:
                #   return l, []
                let why = l.why_coll(colls, level)
                if (why != nil) and (len(why) > 0):
                    return (l, why)

    return (nil, nil)

proc get_circles_thru_all*[DepsT](points: seq[Point[DepsT]]): seq[Circle[DepsT]] =
    return get_node_type_thru_all(points, Circle[DepsT])

# Why points are concyclic.
proc circle_of_and_why*[DepsT](points: seq[Point]; level: int = -1): tuple[a: Circle, b: seq[DepsT]] =
    for c0 in get_circles_thru_all(points):
        for c in c0.equivs():
            if all(points, (p) => p in c.edge_graph):
                let cycls = points.toHashSet().toSeq()
                let why = c.why_cyclic(cycls, level)
                if (why != nil) and (len(why) > 0):
                    return (c, why)

    return (nil, nil)

proc new_val*[DepsT](self: Angle[DepsT]): Measure[DepsT] =
    return newNodeType[Measure[DepsT]]()

proc set_directions*[DepsT](self: Angle[DepsT]; d1, d2: Direction[DepsT]) =
    self.d = (d1, d2)

proc directions*[DepsT](self: Angle[DepsT]): tuple[d1: Direction[DepsT], d2: Direction[DepsT]] =
    let d1, d2 = self.d
    if (d1 == nil) or (d2 == nil):
        return (d1, d2)
    return (d1.rep(), d2.rep())

proc new_val*[DepsT](self: Ratio[DepsT]): Value[DepsT] =
    return newNodeType[Value[DepsT]]()

proc set_lengths*[DepsT](self: Ratio[DepsT]; l1, l2: Length[DepsT]) =
    self.l = (l1, l2)

proc lengths*[DepsT](self: Ratio[DepsT]): tuple[l1: Length[DepsT], l2: Length[DepsT]] =
    let l1, l2 = self.d
    if (l1 == nil) or (l2 == nil):
        return (l1, l2)
    return (l1.rep(), l2.rep())

iterator all_angles*[DepsT](d1, d2: Direction[DepsT]; level: int = -1): tuple[a: Angle[DepsT], b: Table[Node[DepsT], Node[DepsT]], c: Table[Node[DepsT], Node[DepsT]]] =
    let d1s = d1.equivs_upto(level)
    let d2s = d2.equivs_upto(level)

    for ang in d1.rep().neighbors(Angle[DepsT]):
        let d1, d2 = ang.d
        if (d1 in d1s) and (d2 in d2s):
            yield (ang, d1s, d2s)

iterator all_ratios*[DepsT](d1, d2: Direction[DepsT], level: int = -1): tuple[a: Ratio[DepsT], b: Table[Node[DepsT], Node[DepsT]], c: Table[Node[DepsT], Node[DepsT]]] =
    let d1s = d1.equivs_upto(level)
    let d2s = d2.equivs_upto(level)

    for ratio in d1.rep().neighbors(Ratio[DepsT]):
        let d1, d2 = ratio.l
        if (d1 in d1s) and (d2 in d2s):
            yield (ratio, d1s, d2s)

proc val_type*[DepsT](x: Node[DepsT]): typedesc =
    if x of Line[DepsT]:
        return Direction[DepsT]
    if x of Segment[DepsT]:
        return Length[DepsT]
    if x of Angle[DepsT]:
        return Measure[DepsT]
    if x of Ratio[DepsT]:
        return Value[DepsT]
    raiseAssert(fmt"not supported node {type(x)}")