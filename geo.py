import bmesh
import bpy
from mathutils import Vector, Quaternion
from mathutils import geometry

from . import modeler

use_exact = False


class Vert:
    def __init__(self, v1=None):
        if v1 != None:
            self.co = v1.co.copy()


class Edge:
    def __init__(self, e1=None):
        if e1 != None:
            self.verts = [Vert(e1.verts[0]), Vert(e1.verts[1])]
            self.link_faces = [f1.select for f1 in e1.link_faces]

    def calc_length(self):
        return (self.verts[1].co - self.verts[0].co).length


class Face:
    def __init__(self, f1=None):
        if f1 != None:
            self.edges = [Edge(e1) for e1 in f1.edges]
            self.verts = [Vert(v1) for v1 in f1.verts]
            self.normal = f1.normal.copy()
            c1 = f1.calc_center_bounds()
            line_a = c1
            line_b = c1 - f1.normal
            c2 = geometry.intersect_line_plane(line_a, line_b, f1.verts[0].co, f1.normal)
            if c2 == None:
                self.center = c1
            else:
                self.center = c2

    def calc_center_median(self):
        return self.center


class Projection:
    def __init__(self):
        self.current = None
        self.face = None

    def vert(self, v1):
        v2 = Vert()
        v2.co = v1
        return v2

    def edge(self, pv1, pv2):
        e1 = Edge()
        e1.verts = [pv1, pv2]
        e1.link_faces = [True, False]
        return e1

    def create_projection(self, cen, sn, ve1, ve2, ve3, ve4):
        pv1 = self.vert(ve1)
        pv2 = self.vert(ve2)
        pv3 = self.vert(ve3)
        pv4 = self.vert(ve4)
        pe1 = self.edge(pv1, pv2)
        pe2 = self.edge(pv2, pv3)
        pe3 = self.edge(pv3, pv4)
        pe4 = self.edge(pv4, pv1)
        f1 = Face()
        f1.edges = [pe1, pe2, pe3, pe4]
        f1.verts = [pv1, pv2, pv3, pv4]
        f1.center = cen
        f1.normal = sn
        self.face = f1

    def calc_edge(self, p1, p2, offset, normal):
        sn = normal
        sn.normalized()
        mid = (p1 + p2) / 2
        e2 = p2 - p1
        elen = e2.length
        height = e2.cross(sn).normalized()

        h1 = mid + (height * elen / 2)
        h2 = mid - (height * elen / 2)
        ve1 = h2 + (p2 - mid)
        ve2 = h1 + (p2 - mid)
        ve3 = h1 + (p1 - mid)
        ve4 = h2 + (p1 - mid)
        self.create_projection(mid, sn, ve1, ve2, ve3, ve4)

    def calc_space(self, ax, quad, length):
        x, y, z = ax
        ve1, ve2, ve3, ve4 = quad
        cen = (ve1 + ve2 + ve3 + ve4) / 4
        sn = z.normalized()
        self.create_projection(cen, sn, ve1, ve2, ve3, ve4)

def get_bm(context):
    obj = bpy.context.edit_object
    me = obj.data
    bm = bmesh.from_edit_mesh(me)
    return bm


def same_direction(p1, p2):
    if p1.angle(p2) < 0.03:
        return True
    else:
        return False


def vs_center(vs):
    p = Vector()
    for v1 in vs:
        p += v1.co
    p /= len(vs)
    return p


def create_bevel(front, back, p2, dis, cut):
    e1 = (front - p2).normalized()
    e2 = (back - p2).normalized()
    vn = dis * ((e1 + e2) / 2).normalized()
    cir = p2 + vn
    pe1 = vn.project(e1)
    pe2 = vn.project(e2)
    f1 = pe1 - vn
    f2 = pe2 - vn
    p4 = []
    p4.append(cir + f2)
    rot = f2.rotation_difference(f1)
    rot2 = Quaternion(rot.axis, rot.angle / cut)
    for i in range(cut):
        f2.rotate(rot2)
        p4.append(cir + f2)
    return p4


def process_bevel(p, line_only, dis, cut):
    ret = []
    if line_only:
        if len(p) < 2:
            return None
        firsta = p[0]
        lastb = p[-1]
        p2 = p
        pbevel = p2[1:-1]
        front = p2[2:]
        back = p2[:-2]
        p4 = []
        for ip, ifront, iback in zip(pbevel, front, back):
            ires = create_bevel(ifront.co, iback.co, ip.co, dis, cut)
            p4 += ires
        p4 = [firsta.co] + p4 + [lastb.co]
        ret = [modeler.SPoint(a) for a in p4]
    else:
        if len(p) < 3:
            return None
        p2 = p
        front = p2[1:] + [p2[0]]
        back = [p2[-1]] + p2[0:-1]
        p4 = []
        for ip, ifront, iback in zip(p2, front, back):
            ires = create_bevel(ifront.co, iback.co, ip.co, dis, cut)
            p4 += ires
        ret = [modeler.SPoint(a) for a in p4]
    return ret


def create_edges(context, bm, loops):
    es = []
    for item in loops:
        loop = item.loop
        selected = item.select
        pos = [p1.co for p1 in loop]
        vs = []
        for p in pos:
            if p == None:
                continue
            v1 = bm.verts.new(p)
            vs.append(v1)
        vs2 = vs[1:] + [vs[0]]
        pair1 = list(zip(vs, vs2))

        if item.line_only:
            pair1.pop()

        for v1, v2 in pair1:
            e1 = bm.edges.new((v1, v2))
            es.append(e1)

    bm.verts.index_update()
    bm.edges.index_update()
    bm.normal_update()
    bmesh.update_edit_mesh(context.edit_object.data)
    return es


def create_face(context, bm, source, other, loops, viewer, merge):
    source.select = False
    for f1 in other:
        f1.select = False

    fs = []
    for item in loops:
        loop = item.loop
        selected = item.select
        pos = [p1.co for p1 in loop]
        verts = []
        for p in pos:
            if p == None:
                continue
            v1 = bm.verts.new(p)
            verts.append(v1)
        bm.verts.index_update()
        f1 = bm.faces.new(verts)
        # f1.select = True
        f1.select = False
        fs.append(f1)
    bm.faces.index_update()
    bm.normal_update()
    for f1 in fs:
        if viewer == None:
            if same_direction(f1.normal, source.normal) == False:
                f1.select = True
        else:
            if f1.normal.angle(viewer) > (f1.normal * -1).angle(viewer):
                f1.select = True

    bm.select_flush(False)
    bm.select_flush(True)
    bpy.ops.mesh.flip_normals()

    for f1 in fs:
        f1.select = True

    bm.select_flush(True)
    bmesh.update_edit_mesh(context.edit_object.data)

    if merge:
        bpy.ops.mesh.remove_doubles()

    return fs


def cut_face(context, bm, source, other, loops, dissolve, viewer, merge):
    source2 = [source] + other

    es = [e1 for f1 in source2 for e1 in f1.edges]
    es = [e1 for e1 in es if not all(f1.select for f1 in e1.link_faces)]
    fs = create_face(context, bm, source, other, loops, viewer, merge)

    bmesh.ops.delete(bm, geom=source2, context='FACES')

    for f1 in fs:
        for e1 in f1.edges:
            if (e1 in es) == False:
                es.append(e1)

    es = list(set(es))

    bmesh.ops.triangle_fill(bm, use_beauty=False, use_dissolve=dissolve, edges=es)

    bmesh.update_edit_mesh(context.edit_object.data, loop_triangles=True, destructive=True)

    for e1 in bm.edges:
        e1.select = False
    for f1 in fs:
        f1.select = True


def boolean_item(context, bm, f1, source, other, grid_len, propdepth, propoffset, bool_outer_height):
    if bool_outer_height != 0:
        sn = source.normal
        for v1 in f1.verts:
            v1.co += sn * bool_outer_height * grid_len

    ext = bmesh.ops.extrude_face_region(bm, geom=[f1])
    verts = [v for v in ext["geom"] if isinstance(v, bmesh.types.BMVert)]

    depth = source.normal.normalized() * -1 * grid_len * (propdepth + bool_outer_height)
    bmesh.ops.translate(bm, vec=depth, verts=verts)

    f1.select = True
    # new
    bm.select_flush(False)
    bm.select_flush(True)
    bmesh.update_edit_mesh(context.edit_object.data)

    bpy.ops.mesh.select_linked()
    vs = [v1 for v1 in bm.verts if v1.select]
    scale = (propoffset / 5.0) + 1.0
    cen = vs_center(vs)
    for v1 in vs:
        v1.co -= cen
        v1.co *= scale
        v1.co += cen

    return [f1 for f1 in bm.faces if f1.select]


def boolean_cut(context, bm, source, other, loops, grid_len, propdepth, propoffset, bool_intersect, bool_outer_height):
    source.select = False
    for f1 in other:
        f1.select = False
    fsall = []
    for item in loops:
        if item.line_only and len(item.loop) < 3:
            continue
        loop = item.loop
        pos = [p1.co for p1 in loop]
        verts = []
        for p in pos:
            if p == None:
                continue
            v1 = bm.verts.new(p)
            verts.append(v1)
        bm.verts.index_update()
        f1 = bm.faces.new(verts)
        f1.select = False
        bm.faces.index_update()
        bm.normal_update()
        fs = boolean_item(context, bm, f1, source, other, grid_len, propdepth, propoffset, bool_outer_height)
        fsall += fs
        for f in bm.faces:
            f.select = False

    for f in fsall:
        f.select = True

    bm.select_flush(False)
    bm.select_flush(True)
    bmesh.update_edit_mesh(context.edit_object.data, loop_triangles=True, destructive=True)
    bpy.ops.mesh.normals_make_consistent(inside=False)

    if bool_intersect:
        if bpy.app.version >= (2, 91, 0):
            if use_exact:
                bpy.ops.mesh.intersect_boolean(operation='INTERSECT', solver='EXACT')
            else:
                bpy.ops.mesh.intersect_boolean(operation='INTERSECT', solver='FAST')
        else:
            bpy.ops.mesh.intersect_boolean(operation='INTERSECT')
    else:
        if bpy.app.version >= (2, 91, 0):
            if use_exact:
                bpy.ops.mesh.intersect_boolean(operation='DIFFERENCE', solver='EXACT')
            else:
                bpy.ops.mesh.intersect_boolean(operation='DIFFERENCE', solver='FAST')
        else:
            bpy.ops.mesh.intersect_boolean(operation='DIFFERENCE')

