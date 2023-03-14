import math

import blf
import bmesh
import bpy
from bpy.props import (
    FloatProperty,
    IntProperty,
    BoolProperty,
    EnumProperty,
    StringProperty,
    FloatVectorProperty
)
from bpy_extras import view3d_utils
from mathutils import Matrix, Vector, Quaternion
from . import geo
from . import lining
from . import plane


def draw_text_callback(self, context):
    editing = context.mode == 'EDIT_MESH'
    if editing == False:
        return
    if not self.text_handlers:
        return
    else:
        func = self.text_handlers[-1]
        txtall = func(context)
        if txtall == None:
            return
        else:
            left = self.text_pos_x
            sp = self.text_size * 1.7
            top = len(txtall) * sp + 50
            off = 0
            for p in txtall:
                off += sp
                self.draw_text([left, top - off], p)


def draw_3d(self, context):
    editing = context.mode == 'EDIT_MESH'
    if editing == False:
        return

    obj = context.edit_object
    world = obj.matrix_world
    world2 = world.inverted()
    cursor = self.cursor

    pmat = self.get_matrix()
    if pmat == None:
        return

    for item in self.loops:
        item.create_edges()

    if cursor != None:
        coord2 = self.create_cross(context, cursor)
        lining.draw_line(coord2, (0, 0, 0, 1), blend=True, smooth=True, width=4)
        lining.draw_line(coord2, (1, 1, 1, 1), blend=True, smooth=True, width=2)

    if self.coord != None:
        coord2 = self.coord
        lining.draw_line(coord2, self.color_grid_line, blend=True, smooth=True, width=1)

    if self.loops != None and len(self.loops) > 0:
        coord2 = []
        sel_coord = []
        cons = []

        for item in self.loops:
            loop = item.loop
            loop2 = item.loop2
            selected = item.select
            for p1, p2 in zip(loop, loop2):
                if p1 == None or p2 == None:
                    continue

                if selected and self.selection_mode:
                    sel_coord += [world @ p1.co, world @ p2.co]
                else:
                    if item.construction:
                        cons += [world @ p1.co, world @ p2.co]
                    else:
                        coord2 += [world @ p1.co, world @ p2.co]

        lining.draw_line(cons, (1, 0, 0, 1), blend=True, smooth=True, width=1)
        lining.draw_line(coord2, self.color_shape, blend=True, smooth=True, width=1)
        lining.draw_line(sel_coord, (1, 1, 0, 1), blend=True, smooth=True, width=1)

    if self.loops != None and len(self.loops) > 0 and self.selection_mode:
        if self.editshape == None:
            coord2 = []
            pmat2 = pmat.inverted()
            glen = self.grid_len / 4
            for item in self.loops:
                cen = item.center
                if item.select:
                    if cen == None:
                        continue
                    c2 = pmat2 @ cen
                    x1 = pmat @ (c2 + Vector((glen * -1, 0, 0)))
                    x2 = pmat @ (c2 + Vector((glen, 0, 0)))
                    y1 = pmat @ (c2 + Vector((0, glen * -1, 0)))
                    y2 = pmat @ (c2 + Vector((0, glen, 0)))
                    coord2 = [world @ x1, world @ x2, world @ y1, world @ y2]
                    break
            lining.draw_line(coord2, (0, 0, 1, 1), blend=True, smooth=True, width=2)
        else:
            coord2 = []
            coord3 = []
            pmat2 = pmat.inverted()
            glen = self.grid_len / 4
            for item in self.loops:
                if item.select:
                    for i, pt in enumerate(item.loop):
                        selpoint = (item, i)
                        c2 = pmat2 @ pt.co
                        p1 = pmat @ (c2 + Vector((glen * -1, glen * -1, 0)))
                        p2 = pmat @ (c2 + Vector((glen * -1, glen, 0)))
                        p3 = pmat @ (c2 + Vector((glen, glen, 0)))
                        p4 = pmat @ (c2 + Vector((glen, glen * -1, 0)))
                        if selpoint in self.editshape and self.editshape_bevel == None:
                            coord3 += [world @ p1, world @ p2, world @ p2, world @ p3]
                            coord3 += [world @ p3, world @ p4, world @ p4, world @ p1]
                        else:
                            coord2 += [world @ p1, world @ p2, world @ p2, world @ p3]
                            coord2 += [world @ p3, world @ p4, world @ p4, world @ p1]

            lining.draw_line(coord2, (1, 1, 1, 1), blend=True, smooth=True, width=1)
            lining.draw_line(coord3, (0, 0, 1, 1), blend=True, smooth=True, width=2)

    if self.currentloop != None:
        coord2 = []
        if len(self.currentloop) > 0:
            # print(len(self.currentloop))
            loop1 = self.currentloop
            if self.item != None:
                p1, p2 = self.item
                loop2 = loop1[1:] + [SPoint(p1)]
            else:
                loop2 = loop1[1:] + [None]
            for (p1, p2) in zip(loop1, loop2):
                if p1 != None and p2 != None:
                    coord2 += [world @ p1.co, world @ p2.co]

        if self.item != None:
            p1, p2 = self.item
            coord2 += [world @ p1, world @ p2]

        lining.draw_line(coord2, self.color_shape, blend=True, smooth=True, width=1)

    if self.yellow_rect != None and len(self.yellow_rect) > 0:
        coord2 = [world @ a for a in self.yellow_rect]
        lining.draw_line(coord2, (1, 1, 0, 1), blend=True, smooth=True, width=1)

    if self.paste != None:
        coord2 = []
        offset = self.calc_paste_offset()
        for item in ModOperator.clipboard:
            item.create_edges()
            selected = item.select
            for (p1, p2) in zip(item.loop, item.loop2):
                if p1 != None and p2 != None:
                    coord2 += [world @ (p1.co + offset), world @ (p2.co + offset)]

        lining.draw_line(coord2, (1, 1, 0, 1), blend=True, smooth=True, width=1)

    if self.circle != None:
        _, _, loop = self.circle
        if loop != None:
            loop2 = loop[1:] + [loop[0]]
            coord2 = []
            for (p1, p2) in zip(loop, loop2):
                if p1 != None and p2 != None:
                    coord2 += [world @ p1.co, world @ p2.co]

        lining.draw_line(coord2, (1, 1, 0, 1), blend=True, smooth=True, width=1)

    if self.bevel != None:
        loc1, loc2, cut = self.bevel
        coord2 = [world @ loc1, world @ loc2]
        lining.draw_line(coord2, (28 / 255, 240 / 255, 255 / 255, 1), blend=True, smooth=True, width=1)

    if self.main_edge != None:
        p1, p2 = self.main_edge
        coord2 = [world @ p1, world @ p2]
        lining.draw_line(coord2, (0, 1, 1, 1), blend=True, smooth=True, width=1)

    if self.eng != None:
        coord2 = []
        for (p1, p2) in self.eng:
            coord2 += [world @ p1, world @ p2]
        lining.draw_line(coord2, (1, 1, 0, 1), blend=True, smooth=True, width=1)


class SPoint:
    def __init__(self, v):
        self.co = v

    def __repr__(self):
        return '<SP ' + str(self.co) + '>'

    def copy(self):
        p = SPoint(self.co.copy())
        return p


class Shape:
    def __init__(self):
        self.loop = None
        self.loop2 = None
        self.select = False
        self.backup = None
        self.center = None
        self.construction = False
        self.line_only = False

    def __repr__(self):
        s = "<Shape \n"
        for p in self.loop:
            s += str(p) + "\n"
        s += ">\n"
        return s

    def integrity_check_length(self):
        ps = self.loop
        p2 = None
        res = []
        for p in ps:
            if p2 != None:
                if (p.co - p2).length > 0.00001:
                    res.append(p)
                    p2 = p.co
            else:
                res.append(p)
                p2 = p.co

        self.loop = res

    def integrity_check_colinear(self):
        ps = self.loop
        if self.line_only:
            ps = self.loop
            ps2 = ps[1:] + [None]
            ps3 = [None] + ps[:-1]
        else:
            ps = self.loop
            ps2 = ps[1:] + [ps[0]]
            ps3 = [ps[-1]] + ps[:-1]
        rm = []
        for p1, pf, pb in zip(ps, ps2, ps3):
            if pf != None and pb != None:
                m1 = pf.co - p1.co
                m3 = pb.co - p1.co
                if m1.length == 0 or m3.length == 0:
                    print('zero')
                    continue
                if m1.cross(m3).length < 0.000001:
                    if (p1 in rm) == False:
                        rm.append(p1)
        for p in rm:
            self.loop.remove(p)

    def link_point_edges(self, sn):
        if self.line_only:
            ps = self.loop
            ps2 = ps[1:] + [None]
            ps3 = [None] + ps[:-1]
        else:
            ps = self.loop
            ps2 = ps[1:] + [ps[0]]
            ps3 = [ps[-1]] + ps[:-1]

        for p3, p1, p2 in zip(ps3, ps, ps2):
            p1.back = p3
            p1.next = p2
        quo = Quaternion(sn, math.radians(-90))
        for p in self.loop:
            if p.next != None:
                m1 = p.next.co - p.co
                m1.rotate(quo)
                p.inward = m1.normalized()
            else:
                if p.back != None:
                    p.inward = p.back.inward

        for p in self.loop:
            if p.back != None and p.next != None:
                m1 = p.back.co - p.co
                m2 = p.next.co - p.co
                cr = m1.cross(m2)
                if cr.length == 0 or sn.length == 0:
                    p.small = True
                    continue
                if cr.angle(sn) < 0.1:
                    p.small = True
                else:
                    p.small = False
            else:
                p.small = True

    def create_edges(self):
        p = self.loop
        if p == None:
            return
        if len(p) < 2:
            return

        if self.line_only:
            p2 = p[1:] + [None]
        else:
            p2 = p[1:] + [p[0]]

        self.loop2 = p2

    def save_center(self):
        if self.loop == None:
            return
        p1 = self.loop[0]
        self.center = p1.co

    def copy(self):
        s = Shape()
        s.select = self.select
        s.backup = self.backup
        s.center = self.center
        s.construction = self.construction
        s.line_only = self.line_only
        s.loop = [p1.copy() for p1 in self.loop]
        return s

    def calc_center(self):
        if self.loop == None:
            return
        ps = [p1.co for p1 in self.loop]
        total = Vector()
        for p in ps:
            total += p
        cen = total / len(ps)
        self.center = cen


def get_collection_objects():
    if 'GM' in bpy.data.collections:
        return bpy.data.collections['GM'].objects
    else:
        return []


class ModOperator(bpy.types.Operator):
    bl_idname = "mesh.modeler_operator"
    bl_label = "HardSurface modeler"
    bl_options = {"REGISTER", "UNDO"}

    operation_mode: EnumProperty(
        items=[('ngon', 'N-gon', '', '', 0),
               ('newface', 'Create Face', '', '', 1),
               ('boolcut', 'Boolean Cut', '', '', 2),
               ('linepipe', 'Edge Pipe', '', '', 3),
               ],
        name="Operation",
        default='boolcut')

    bool_triangle: BoolProperty(
        name="Fill by triangles",
        description="Fill by triangles",
        default=False
    )

    prop_edge_extrude: FloatProperty(
        name="Depth of extrude",
        description="Depth extrude",
        default=0.01,
        step=0.05
    )

    prop_edge_no_extrude: BoolProperty(
        name="Disable extrude (Create edges only)",
        description="Create edges only",
        default=False,
    )

    propdepth: FloatProperty(
        name="Depth",
        description="Depth for Boolean",
        default=1.0
    )

    propoffset: FloatProperty(
        name="Boundary offset",
        description="Offset at boundary",
        default=0.03
    )

    prop_new_object: BoolProperty(
        name="Separated object",
        description="Create a new face as separated object",
        default=False
    )

    bool_split_all: BoolProperty(
        name="Knife cut through",
        description="Line Cut - cut through whole mesh",
        default=False
    )

    bool_split_break: BoolProperty(
        name="Split geometry",
        description="Split and disconnect the geometry",
        default=False
    )

    bool_flip_normal: BoolProperty(
        name="Flip Normal",
        description="Flip Normal of new face",
        default=False
    )

    bool_merge_edges: BoolProperty(
        name="Merge Edges",
        description="Merge overlapping edges",
        default=False
    )

    bool_isolation: BoolProperty(
        name="Mesh Isolation",
        description="Isolate the mesh (avoid blender's bug)",
        default=True
    )

    bool_outer_height: FloatProperty(
        name="Outer Height",

        default=0
    )

    bool_intersect: BoolProperty(
        name="Intersection",
        default=False
    )

    bool_exact_solver: BoolProperty(
        name="Blender Exact Solver",
        default=False
    )

    create_edge_bevel: IntProperty(
        name="Edges Bevel",
        default=4
    )

    create_edge_bevel_dis: FloatProperty(
        name="Edges Bevel Width",
        default=0.1
    )

    bool_make_pipe: BoolProperty(
        name="Create Pipe",
        default=True
    )

    pipe_bevel_depth: FloatProperty(
        name="Pipe Thickness",
        description="Pipe Thickness",
        default=3.0,
        step=100
    )

    bool_pipe_smooth: BoolProperty(
        name="Smooth Pipe",
        description="Smooth the Pipe by Nurbs",
        default=False
    )

    pipe_bevel_resolution: IntProperty(
        name="Pipe Bevel Resolution",
        description="Pipe Bevel Resolution",
        default=4
    )

    pipe_render_resolution: IntProperty(
        name="Pipe Render Resolution",
        description="Pipe Render Resolution",
        default=16
    )

    bool_pipe_curve: BoolProperty(
        name="Create new curve object",
        description="Create new curve object",
        default=False
    )

    bool_pipe_fill: BoolProperty(
        name="Fill the geometry",
        description="Fill the geometry",
        default=True
    )

    bool_surface_only: BoolProperty(
        name="Slice Surface Only",
        description="Slice Surface Only",
        default=False
    )

    prop_add_text: StringProperty(
        name="Text value",
        description="Text value",
        default='Text'
    )

    prop_add_text_size: FloatProperty(
        name="Text font size",
        description="Text font size",
        default=0.5,
        max=10,
        min=0.001,
        step=0.1
    )

    prop_add_text_extrude: FloatProperty(
        name="Extrude depth",
        description="Extrude depth",
        default=0
    )

    prop_add_text_offset: FloatVectorProperty(
        name="Offset",
        description="Offset",
        default=Vector()
    )

    prop_add_text_font_path: StringProperty(
        name="Font file",
        description="The font file",
        subtype="FILE_PATH",
        default=bpy.context.preferences.filepaths.font_directory
    )

    prop_add_text_convert_mesh: BoolProperty(
        name="Merge text mesh",
        description="Convert text to mesh and join",
        default=False
    )

    prop_add_text_bool_cut: BoolProperty(
        name="Boolean Cut the Text",
        description="Boolean Cut the Text",
        default=False
    )

    prop_add_text_bool_exact: BoolProperty(
        name="Blender's Exact Solver",
        description="Blender's Exact Solver",
        default=False
    )

    prop_object_adder_index: IntProperty(
        name="Object Index (GM collection)",
        description="Object Index (GM collection)",
        default=0,
        min=0,
    )

    prop_object_adder_rotation: FloatVectorProperty(
        name="Rotation",
        description="Rotation (normal axis)",
        default=Vector(),
        subtype='XYZ',
        step=100
    )

    prop_object_adder_offset: FloatVectorProperty(
        name="Offset",
        description="Offset",
        default=Vector(),
        subtype='XYZ',
        step=0.5
    )

    prop_object_adder_scale: FloatVectorProperty(
        name="Scale",
        description="Scale",
        default=Vector((1, 1, 1)),
        subtype='XYZ'
    )

    prop_object_adder_separated: BoolProperty(
        name="Separated Object",
        description="Separated Object",
        default=False
    )

    prop_object_adder_separated_linked: BoolProperty(
        name="Linked Object",
        description="Linked Object",
        default=False
    )

    prop_object_adder_fixed_height_value: FloatProperty(
        name="Height (for rectangle input)",
        description="Height (for rectangle input)",
        default=0.3
    )

    prop_object_adder_picker: StringProperty(
        name="Pick",
        description="Pick",
    )

    prop_object_adder_textmode: EnumProperty(
        # (identifier, name, description, icon, number)
        items=[('Object', 'Object', '', '', 0),
               ('Text', 'Text', '', '', 1)
               ],
        name="Tool",
        default='Object')

    prop_object_adder_use_picker: EnumProperty(
        # (identifier, name, description, icon, number)
        items=[('Object Index', 'Object Index', '', '', 0),
               ('Picker', 'Picker', '', '', 1)
               ],
        name="Mode",
        default='Object Index')

    shader = None
    handle3d = None
    handle = None
    mode_A = True
    mode_new_face = False
    size = 10
    size_rel = 5
    circle_cut = 12
    clipboard = []
    clipboard_sn = None
    clipboard_main = None
    tmp_plane_save = None
    coord_mode = False

    def __init__(self):
        self.space = None
        self.space_point = None
        self.space_center = None
        self.points = []
        self.cursor = None
        self.drawing = False
        self.item = None
        self.firstv = None
        self.lastv = None
        self.source = None
        self.guide = None
        self.guide_points = None
        self.geo_points = []
        self.geo_edges = []

        self.coord = None
        self.aligning = None

        self.grid_len = None
        self.scale_up = 0
        self.scale_up_rel = 0

        self.loops = []
        self.currentloop = []
        self.main_edge = None
        self.shift_vert = None

        self.handlers = []
        self.text_handlers = []
        self.selection = None
        self.yellow_rect = None
        self.select_area = None

        self.paste = None
        self.new_center = None
        self.rotate_deg = 0
        self.rotation = None
        self.array = False
        self.array_count = 0
        self.array_circle = False
        self.scale = None
        self.text_input = ''

        self.selection_mode = False

        self.circle = None
        self.other_source = []
        self.eng = None
        self.text_size = 14

        self.bevel = None
        self.editshape = None
        self.editshape_loc = None
        self.editshape_bevel = None

        self.scale_mode = None

        self.virtualface = False
        self.projection = None
        self.snap_disable = False

        self.draw_down = False

        self.elevation_num = 0

        self.con_lines = []
        self.con_mode = None
        self.intersection = []
        self.helping = False
        self.help_scroll = 0
        self.show_keys = True
        self.keys_list = []
        self.keys_index = 0

        self.snap_to_shape = False
        self.base_location = None
        self.projection_backup = None
        self.pkey_count = 0

        pass

    def execute(self, context):
        self.finish_action(context)
        return {'FINISHED'}

    def draw(self, context):
        layout = self.layout
        layout.label(text="Select the operation: ")
        row = layout.row()
        layout.prop(self, "operation_mode", expand=False)
        layout.label(text="")

        if self.operation_mode == 'boolcut' and self.virtualface == False:
            layout.prop(self, "bool_isolation")

        if self.operation_mode == 'boolcut':
            layout.prop(self, "bool_exact_solver")
            layout.prop(self, "bool_intersect")
            layout.prop(self, "propdepth")
            layout.prop(self, "bool_outer_height")
            layout.prop(self, "propoffset")

        elif self.operation_mode == 'ngon':
            layout.prop(self, "bool_triangle")

        elif self.operation_mode == 'lineart':
            layout.label(text="This function is designed for Edge Shader.")
            layout.prop(self, 'prop_edge_extrude')
            layout.prop(self, 'prop_edge_no_extrude')


        elif self.operation_mode == 'boolslice':
            layout.prop(self, "bool_exact_solver")
            layout.prop(self, "bool_surface_only")
            layout.prop(self, "propdepth")
            layout.prop(self, "propoffset")

        elif self.operation_mode == 'linesplit':
            layout.prop(self, "bool_split_all")
            layout.prop(self, "bool_split_break")

        elif self.operation_mode == 'newface':
            layout.prop(self, "prop_new_object")
            layout.prop(self, "bool_flip_normal")
            layout.prop(self, "bool_merge_edges")

        elif self.operation_mode == 'linepipe':
            layout.prop(self, "create_edge_bevel")
            layout.prop(self, "create_edge_bevel_dis")
            layout.prop(self, "bool_make_pipe")
            layout.prop(self, "bool_pipe_smooth")
            layout.prop(self, "pipe_bevel_depth")

            layout.prop(self, "pipe_bevel_resolution")
            layout.prop(self, "bool_pipe_curve")

        elif self.operation_mode == 'addtext':
            row = layout.row()
            row.prop(self, "prop_object_adder_textmode", expand=False)

            if self.prop_object_adder_textmode == 'Object':
                index = self.prop_object_adder_index
                gms = get_collection_objects()
                if index < len(gms):
                    obj = gms[index]
                    layout.label(text=obj.name)
                else:
                    layout.label(text='')

                row = layout.row()
                row.prop(self, "prop_object_adder_use_picker", expand=False)

                if self.prop_object_adder_use_picker == "Picker":
                    if 'GM' in bpy.data.collections:
                        gmc = bpy.data.collections['GM']
                        row = layout.row()
                        row.prop_search(self, "prop_object_adder_picker", gmc, "objects")
                else:
                    layout.prop(self, "prop_object_adder_index")

                layout.prop(self, "prop_object_adder_rotation")
                layout.prop(self, "prop_object_adder_offset")
                layout.prop(self, "prop_object_adder_scale")
                row = layout.row()
                row.prop(self, "prop_object_adder_separated")
                row.prop(self, "prop_object_adder_separated_linked")
                layout.prop(self, "prop_object_adder_fixed_height_value")
            else:
                layout.prop(self, "prop_add_text")
                layout.prop(self, "prop_add_text_size")
                layout.prop(self, "prop_add_text_extrude")
                layout.prop(self, "prop_add_text_offset")
                layout.prop(self, "prop_add_text_font_path")
                layout.prop(self, "prop_add_text_convert_mesh")
                layout.prop(self, "prop_add_text_bool_cut")
                layout.prop(self, "prop_add_text_bool_exact")

    @classmethod
    def remove_draw(cls):
        if ModOperator.handle3d != None:
            bpy.types.SpaceView3D.draw_handler_remove(
                ModOperator.handle3d, 'WINDOW')
            ModOperator.handle3d = None

        if ModOperator.handle != None:
            bpy.types.SpaceView3D.draw_handler_remove(
                ModOperator.handle, 'WINDOW')
            ModOperator.handle = None

    @classmethod
    def poll(cls, context):
        active_object = context.active_object
        selecting = active_object is not None and active_object.type == 'MESH'
        editing = context.mode == 'EDIT_MESH'
        is_vert_mode, is_edge_mode, is_face_mode = context.tool_settings.mesh_select_mode
        return selecting and editing

    def get_3d_cursor(self, context, event):
        mouse_pos = event.mouse_region_x, event.mouse_region_y
        region = bpy.context.region
        region3D = bpy.context.space_data.region_3d
        view_vector = view3d_utils.region_2d_to_vector_3d(
            region, region3D, mouse_pos)
        view_point = view3d_utils.region_2d_to_origin_3d(
            region, region3D, mouse_pos)
        world_loc = view3d_utils.region_2d_to_location_3d(region,
                                                          region3D, mouse_pos, view_vector)
        return view_vector, view_point, world_loc

    def get_source_real(self, context, bm):
        source = None
        other_source = []
        found = None

        for f1 in bm.faces:
            if f1.select:
                found = True
                if source == None:
                    source = f1
                else:
                    if source.normal.angle(f1.normal) < 0.03:
                        other_source.append(f1)
                    else:
                        print('non-coplanar')
                        return None, None

        if found:
            return source, other_source
        else:
            return None, None

    def get_view_matrix(self):
        return bpy.context.space_data.region_3d.view_matrix

    def is_ortho(self):
        return bpy.context.space_data.region_3d.view_perspective == 'ORTHO'

    def save_source(self, context, bm):
        self.source = None
        self.other_source = []
        found = False

        for f1 in bm.faces:
            if f1.select:
                found = True
                if self.source == None:
                    self.source = geo.Face(f1)
                else:
                    if self.source.normal.angle(f1.normal) < 0.03:
                        self.other_source.append(geo.Face(f1))
                    else:
                        print('non-coplanar')
                        return None
        return found

    def save_source_edge(self, context, bm):
        self.source = None
        self.other_source = []

        sel = None
        for e1 in bm.edges:
            if e1.select:
                if len(e1.verts) != 2:
                    return None
                else:
                    sel = e1
                    break

        if sel != None:
            bpy.ops.mesh.select_linked()
            bpy.ops.mesh.normals_make_consistent(inside=False)
            self.deselect_all_thing(bm, context)
            sel.select_set(True)
            bm.select_flush_mode()
            context.edit_object.data.update()

        if sel != None:
            e1 = sel
            sn = Vector()
            if len(e1.link_faces) == 2:
                f1 = e1.link_faces[0]
                f2 = e1.link_faces[1]
                sn = (f1.normal + f2.normal) / 2
                sn.normalize()
            else:
                sn = (e1.verts[0].normal + e1.verts[1].normal) / 2
                sn.normalize()

            main = e1.verts[1].co - e1.verts[0].co
            if main.length == 0:
                return None

            if sn.length == 0:
                sn = main.cross(Vector((0, 0, 1)))
                if sn.length == 0:
                    sn = main.cross(Vector((1, 0, 0)))
                sn.normalize()

            return self.save_plane(e1.verts[0].co, e1.verts[1].co, sn)
        else:
            return None

    def save_plane(self, p1, p2, sn):
        self.main_edge = (p1, p2)
        self.projection = geo.Projection()
        self.projection.calc_edge(p1, p2, self.propoffset, sn)
        self.source = self.projection.face
        return True

    def save_source_vert(self, context, bm):
        self.source = None
        self.other_source = []

        vs = [v1 for v1 in bm.verts if v1.select]

        if len(vs) != 0:
            bpy.ops.mesh.select_linked()
            bpy.ops.mesh.normals_make_consistent(inside=False)
            self.deselect_all_thing(bm, context)
            for v1 in vs:
                v1.select_set(True)
            bm.select_flush_mode()
            context.edit_object.data.update()

        if len(vs) == 1:
            v1 = vs[0]
            sn = v1.normal
            p1 = Vector((1, 0, 0))
            p2 = Vector((-1, 0, 0))
            pz = Vector((0, 0, 1))
            deg2 = p1.rotation_difference(sn.cross(pz))
            deg = pz.rotation_difference(sn)
            p1.rotate(deg2)
            p1.rotate(deg)
            p2.rotate(deg2)
            p2.rotate(deg)
            return self.save_plane(v1.co + p1, v1.co + p2, sn)

        elif len(vs) == 2:
            v1 = vs[0]
            v2 = vs[1]
            sn = (v1.normal + v2.normal) / 2
            main = v2.co - v1.co
            if main.length == 0:
                return None

            sn.normalize()

            if sn.length == 0:
                sn = main.cross(Vector((0, 0, 1)))
                if sn.length == 0:
                    sn = main.cross(Vector((1, 0, 0)))
                sn.normalize()

            p1 = main.cross(sn)
            sn = main.cross(p1)

            sn.normalize()

            self.save_plane(v1.co, v2.co, sn)
            self.update_guide(context)
            return True

        elif len(vs) >= 3:
            vs = vs[:3]
            vs2 = [vs[-1]] + vs[:-1]
            vslen = [(v1, v2) for (v1, v2) in zip(vs, vs2)]
            vslen = [(v1, v2, (v2.co - v1.co).length) for (v1, v2) in vslen]

            maxp = max(vslen, key=lambda e: e[2])

            v1, v2, _ = maxp
            vset = set(vs)
            vset2 = set([v1, v2])
            vremind = vset - vset2
            v3 = list(vremind)[0]

            sn = (v2.co - v1.co).cross(v3.co - v1.co)
            main = v2.co - v1.co
            if main.length == 0:
                return None

            if sn.length == 0:
                sn = main.cross(Vector((0, 0, 1)))
                if sn.length == 0:
                    sn = main.cross(Vector((1, 0, 0)))
                sn.normalize()
            return self.save_plane(v1.co, v2.co, sn)
        return None

    def save_projection_backup(self):
        if self.projection_backup == None:
            main_hori, main_vert = self.space_point
            n1 = self.source.normal.copy()
            cen = self.source.center
            norm = (cen - n1, cen + n1)

            h1, h2 = main_hori
            v1, v2 = main_vert
            m1 = h2 - h1
            m2 = v2 - v1
            mlen = max(m1.length, m2.length)
            mlen = max(mlen, n1.length)
            mlen = mlen / 2

            h1 = (h1 + h2) / 2 - m1.normalized() * mlen
            h2 = (h1 + h2) / 2 + m1.normalized() * mlen
            v1 = (v1 + v2) / 2 - m2.normalized() * mlen
            v2 = (v1 + v2) / 2 + m2.normalized() * mlen
            norm = (cen - n1 * mlen, cen + n1 * mlen)

            main_hori = (h1, h2)
            main_vert = (v1, v2)

            self.projection_backup = (main_hori, main_vert, norm)

            return self.projection_backup
        else:
            return self.projection_backup

    def deselect_all_thing(self, bm, context):
        for v1 in bm.faces:
            v1.select_set(False)
        for e1 in bm.edges:
            e1.select_set(False)
        for v1 in bm.verts:
            v1.select_set(False)
        bm.select_flush_mode()
        context.edit_object.data.update()

    def viewplane_source(self, context, event):
        obj = bpy.context.object

        vs = [Vector(f[:]) for f in obj.bound_box]
        length = max((vs[0] - vs[1]).length, (vs[0] - vs[2]).length,
                     (vs[0] - vs[3]).length)

        if length == 0:
            return None

        length *= 0.6
        length = length - (length % 0.5)

        vcen = self.get_view_near_bound(context, event)
        ax, quad = self.get_view_space(context, event, vcen, length)

        obj = context.edit_object
        world = obj.matrix_world

        self.projection = geo.Projection()
        self.projection.calc_space(ax, quad, length)
        self.source = self.projection.face
        return True

    def get_view_near_bound(self, context, event):
        vmat = self.get_view_matrix()
        vm2 = vmat.inverted()
        obj = context.edit_object
        world = obj.matrix_world
        world2 = world.inverted()

        obj = bpy.context.object
        vs = [Vector(f[:]) for f in obj.bound_box]
        vs2 = [vmat @ (world @ v) for v in vs]
        mv = max(vs2, key=lambda e: e.z)
        vcen = mv.project(Vector((0, 0, 1)))

        return vcen

    def get_view_space(self, context, event, vcen, length):
        vmat = self.get_view_matrix()
        vm2 = vmat.inverted()
        obj = context.edit_object
        world = obj.matrix_world
        world2 = world.inverted()

        p = length

        ve1 = vcen + Vector((p, p, 0))
        ve2 = vcen + Vector((p, p * -1, 0))
        ve3 = vcen + Vector((p * -1, p * -1, 0))
        ve4 = vcen + Vector((p * -1, p, 0))

        core = world2 @ (vm2 @ Vector((0, 0, 0)))

        x = world2 @ (vm2 @ Vector((1, 0, 0))) - core
        y = world2 @ (vm2 @ Vector((0, 1, 0))) - core
        z = world2 @ (vm2 @ Vector((0, 0, 1))) - core

        quad = (world2 @ (vm2 @ ve1), world2 @ (vm2 @ ve2),
                world2 @ (vm2 @ ve3), world2 @ (vm2 @ ve4))

        ax = (x, y, z)
        return ax, quad

    def update_guide_exe(self, context, ret):
        guide, gridlen, main_hori, main_vert, cen = ret
        m1, m2 = main_hori
        v1, v2 = main_vert
        self.grid_len = gridlen
        self.guide = guide
        self.space = (m2 - m1, v2 - v1)
        self.space_center = cen
        self.space_point = (main_hori, main_vert)
        self.coord = []
        obj = context.edit_object
        world = obj.matrix_world
        for (p1, p2) in self.guide:
            e1 = world @ p1
            e2 = world @ p2
            self.coord += [e1, e2]

    def update_guide(self, context):
        source = [self.source] + self.other_source
        sn = self.source.normal.copy()
        if ModOperator.mode_A:
            size = ModOperator.size
        else:
            size = ModOperator.size_rel

        mode_A = ModOperator.mode_A

        coord_mode = ModOperator.coord_mode

        ret = plane.get_plane(source, sn, self.main_edge,
                              size, mode_A, self.scale_up, self.scale_up_rel, self.shift_vert,
                              coord_mode, self)

        if ret == None:
            return None

        self.update_guide_exe(context, ret)
        return True

    def lay_on(self, cur, p1, p2):
        a1 = cur - p1
        a2 = p2 - p1
        e1 = a1.project(a2)
        n1 = (e1 - a2).length
        n2 = e1.length
        if (n1 + n2) - a2.length < 0.0001:
            return True
        else:
            return False

    def snap(self, context, cursor):
        cursor2 = cursor.copy()

        if self.grid_len == None:
            return
        limit = self.grid_len / 2

        nearedge = []
        for (p1, p2) in self.guide:
            a1 = cursor2 - p1
            a2 = p2 - p1
            e1 = a1.project(a2)
            arm = a1 - e1
            if arm.length < limit:
                nearedge.append((p1, p2, arm.length))

        for (p1, p2) in self.geo_edges:
            if self.lay_on(cursor2, p1, p2) == False:
                continue
            a1 = cursor2 - p1
            a2 = p2 - p1
            e1 = a1.project(a2)
            arm = a1 - e1
            if arm.length < limit:
                nearedge.append((p1, p2, arm.length))

        nearpoints = []
        if self.guide_points != None:
            for p1 in self.guide_points:
                if (p1 - cursor2).length < limit:
                    nearpoints.append(p1)

        edges2 = sorted(nearedge, key=lambda e: e[2])

        for (p1, p2, arm) in reversed(edges2):
            a1 = cursor2 - p1
            a2 = p2 - p1
            e1 = a1.project(a2)
            cursor2 = e1 + p1

        for p1 in nearpoints:
            cursor2 = p1

        c3 = cursor.copy()

        if self.geo_points != None and len(self.geo_points) > 0:
            nearsec = sorted(self.geo_points, key=lambda e: (e - c3).length)
            topn = nearsec[0]
            if (topn - c3).length < limit:
                cursor2 = topn

        if len(self.intersection) > 0:
            nearsec = sorted(self.intersection, key=lambda e: (e - c3).length)
            topn = nearsec[0]
            if (topn - c3).length < limit:
                cursor2 = topn

        self.aligning = edges2
        return cursor2

    def draw_circle_point(self, context, cursor):
        p1, p2, loop = self.circle
        if p1 == None:
            p1 = cursor
            self.circle = (p1, p2, None)
            return
        else:
            if loop != None:
                p2 = cursor
                self.circle = (p1, cursor, None)
                self.draw_circle(context)
                sp = Shape()
                sp.loop = loop
                sp.center = p1

                if self.con_mode != None:
                    sp.construction = True

                self.loops.append(sp)
                self.circle = (None, None, None)
                self.calc_intersection()
                return

    def rotate(self, l, n):
        return l[n:] + l[:n]

    def draw_circle(self, context):
        if self.circle == None:
            return

        p1, p2, _ = self.circle
        if p1 == None or p2 == None:
            return

        if self.space == None:
            return
        main, vertical = self.space
        z_axis = main.cross(vertical)
        size = ModOperator.circle_cut
        deg = 360 / size
        ang = math.radians(deg)
        e1 = p2 - p1
        e1 = e1.normalized() * (e1.length / math.cos(ang / 2))
        pp = []
        for i in range(size):
            a2 = (i + 0.5) * ang
            rot = Quaternion(z_axis, a2)
            e2 = e1.copy()
            e2.rotate(rot)
            p3 = p1 + e2
            pp.append(p3)

        pp3 = [SPoint(d1) for d1 in pp]
        self.circle = (p1, p2, pp3)

    def insert_point(self, context, cursor, pencil):
        if self.drawing:
            if self.item == None:
                return

            if pencil == False:
                if (cursor - self.firstv).length < self.grid_len / 5:
                    self.close_edge(context)
                    return

            v1, v2 = self.item
            self.currentloop.append(SPoint(v1))
            v3 = cursor
            self.lastv = v2
            self.item = [v2, v3]

        else:
            self.drawing = True
            v1 = cursor
            v2 = cursor
            self.firstv = v1
            self.item = (v1, v2)
            self.guide_points = [v1]

    def move_edge(self, context, cursor):
        if self.item == None:
            return
        v1, v2 = self.item
        v2 = cursor
        self.item = (v1, v2)

    def integrity(self, line_only):
        loop = self.currentloop
        if line_only == False and len(loop) < 3:
            return False

        for item in loop:
            if item == None:
                return False
            if item.co == None:
                return False
        return True

    def close_edge(self, context):
        if self.item == None or self.lastv == None or self.firstv == None:
            pass
        else:
            p = SPoint(self.lastv)
            self.currentloop.append(p)

            if self.integrity(False):
                sp = Shape()
                sp.loop = self.currentloop
                # sp.calc_center()
                sp.save_center()
                self.loops.append(sp)

        self.currentloop = []
        self.item = None
        self.firstv = None
        self.lastv = None
        self.drawing = False

    def end_without_close_edge(self, context):
        if self.item == None or self.lastv == None or self.firstv == None:
            pass
        else:
            if self.integrity(True):
                p = SPoint(self.lastv)
                self.currentloop.append(p)

                sp = Shape()
                sp.line_only = True
                sp.loop = self.currentloop
                sp.save_center()
                self.loops.append(sp)

        self.currentloop = []
        self.item = None
        self.firstv = None
        self.lastv = None
        self.drawing = False

    def add_conline(self, context):
        if self.item == None:
            return

        if len(self.currentloop) > 0:
            p = SPoint(self.lastv)
            self.currentloop.append(p)

            sp = Shape()
            sp.loop = self.currentloop
            sp.save_center()
            sp.construction = True
            self.loops.append(sp)
            self.calc_intersection()

        self.currentloop = []
        self.item = None
        self.firstv = None
        self.lastv = None
        self.drawing = False

    def inter_test(self, p, q, r, s):
        rs = r.cross(s)
        qp = q - p
        if rs == 0:
            return None
        t = qp.cross(s) / rs
        u = qp.cross(r) / rs
        if t >= 0 and t <= 1 and u >= 0 and u <= 1:
            return p + t * r

    def intersect(self, mat, mat2, sp1, sp2):
        secs = []
        sp1.create_edges()
        sp2.create_edges()

        for (p1, p2) in zip(sp1.loop, sp1.loop2):
            if p1 == None or p2 == None:
                continue
            for (p3, p4) in zip(sp2.loop, sp2.loop2):
                if p3 == None or p4 == None:
                    continue
                ap1 = mat2 @ p1.co
                ap2 = mat2 @ p2.co
                ap3 = mat2 @ p3.co
                ap4 = mat2 @ p4.co
                bp1 = Vector((ap1.x, ap1.y))
                bp2 = Vector((ap2.x, ap2.y))
                bp3 = Vector((ap3.x, ap3.y))
                bp4 = Vector((ap4.x, ap4.y))
                res = self.inter_test(bp1, bp3, bp2 - bp1, bp4 - bp3)
                if res != None:
                    sec = mat @ Vector((res.x, res.y, 0))
                    secs.append(sec)
        return secs

    def calc_intersection(self):
        mat = self.get_matrix()
        if mat == None:
            return
        mat2 = mat.inverted()
        cons = [item for item in self.loops if item.construction]
        secs = []
        for i, item in enumerate(cons):
            for i2, item2 in enumerate(cons):
                if i < i2:
                    secs += self.intersect(mat, mat2, item, item2)
        self.intersection = secs

        for item in cons:
            for p1 in item.loop:
                self.intersection.append(p1.co)

    def finish_action(self, context):
        ModOperator.remove_draw()
        bm = geo.get_bm(context)

        is_vert_mode, is_edge_mode, is_face_mode = context.tool_settings.mesh_select_mode

        if is_vert_mode:
            source = None
        elif is_edge_mode:
            source = None
        else:
            if self.virtualface:
                source = None
            else:
                source, other = self.get_source_real(context, bm)
                if source == None and other == None:
                    return

        if source == None:
            if self.source == None:
                return
            else:
                source = self.source
                other = []

        loops2 = [item for item in self.loops if item.construction == False and item.line_only == False]

        if self.operation_mode == 'newface':
            self.deselect_all_thing(bm, context)

            vm2 = self.get_view_matrix().inverted()
            vn = vm2 @ Vector((0, 0, -1))
            vn.normalize()

            fs = geo.create_face(context, bm, source, other, loops2, vn, self.bool_merge_edges)
            self.source.select = False
            if self.bool_flip_normal:
                for face in fs:
                    face.normal_flip()
                bm.normal_update()
                bmesh.update_edit_mesh(context.edit_object.data)

            if self.prop_new_object:
                self.make_new_object(bm, fs)

        elif self.operation_mode == 'ngon':
            if self.virtualface:
                return

            if self.bool_triangle:
                geo.cut_face(context, bm, source, other, loops2, False, None, False)
            else:
                geo.cut_face(context, bm, source, other, loops2, True, None, False)

        elif self.operation_mode == 'boolcut':
            bpy.ops.mesh.select_linked()
            bpy.ops.mesh.normals_make_consistent(inside=False)

            if self.bool_isolation and self.virtualface == False:

                prev_show = []
                for face in bm.faces:
                    if face.select == False:
                        if face.hide == False:
                            face.hide_set(True)
                            prev_show.append(face)

                self.deselect_all_thing(bm, context)
                geo.boolean_cut(context, bm, source, other, loops2, self.grid_len, self.propdepth, self.propoffset,
                                self.bool_intersect, self.bool_outer_height)

                for face in prev_show:
                    face.hide_set(False)

                bmesh.update_edit_mesh(context.edit_object.data)
            else:
                self.deselect_all_thing(bm, context)
                geo.boolean_cut(context, bm, source, other, loops2, self.grid_len, self.propdepth, self.propoffset,
                                self.bool_intersect, self.bool_outer_height)

        elif self.operation_mode == 'linepipe':
            loops3 = [item for item in self.loops if item.construction == False]
            if len(loops3) == 0:
                return

            self.deselect_all_thing(bm, context)
            if self.bool_make_pipe:
                bm2 = bmesh.new()
                fill = any([item.line_only for item in loops3])
                if self.create_edge_bevel > 0:
                    if self.bool_pipe_smooth == False:
                        loops3 = self.finish_action_bevel(loops3)
                geo.create_edges(context, bm2, loops3)
                self.make_pipe(context, bm2, fill)
            else:
                if self.create_edge_bevel > 0:
                    loops3 = self.finish_action_bevel(loops3)
                geo.create_edges(context, bm, loops3)

    def make_new_object(self, bm, fs):
        sel = bpy.context.selected_editable_objects.copy()
        for f1 in fs:
            f1.select_set(True)
        bm.select_flush_mode()
        bpy.ops.mesh.separate(type='SELECTED')
        obj = None
        for p in bpy.context.selected_editable_objects:
            if p in sel:
                pass
            else:
                obj = p
        if obj == None:
            return
        obj.modifiers.clear()
        bpy.ops.object.editmode_toggle()
        bpy.context.view_layer.objects.active = obj
        obj.select_set(True)
        bpy.ops.object.editmode_toggle()
        me = obj.data
        bm2 = bmesh.from_edit_mesh(me)
        for f1 in bm2.faces:
            f1.select = True
        bmesh.update_edit_mesh(me)

    def finish_action_bevel(self, loops2):
        cut = self.create_edge_bevel
        dis = self.create_edge_bevel_dis
        loops3 = []
        for item in loops2:
            item2 = item.copy()
            ret = geo.process_bevel(item2.loop, item2.line_only, dis, cut)
            if ret == None:
                continue
            item2.loop = ret
            loops3.append(item2)
        return loops3

    def make_pipe_filling(self, obj):
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.edge_face_add()
        bpy.ops.object.mode_set(mode='OBJECT')

    def make_pipe(self, context, bm, fill):
        me = bpy.data.meshes.new('modeler_pipe')
        new_obj = bpy.data.objects.new('modeler_pipe', me)
        new_obj.location = (0, 0, 0)
        bpy.context.collection.objects.link(new_obj)
        new_obj.matrix_world = context.edit_object.matrix_world
        bm.to_mesh(me)
        current = context.edit_object
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.select_all(action='DESELECT')
        new_obj.select_set(True)
        context.view_layer.objects.active = new_obj
        bpy.ops.object.convert(target='CURVE')

        obj_data = new_obj.data
        obj_data.fill_mode = 'FULL'
        obj_data.bevel_depth = self.pipe_bevel_depth * 0.01
        obj_data.resolution_u = self.pipe_render_resolution
        obj_data.render_resolution_u = self.pipe_render_resolution
        obj_data.bevel_resolution = self.pipe_bevel_resolution

        if self.bool_pipe_smooth:
            obj_data.twist_smooth = 10
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.curve.spline_type_set(type='NURBS')

            for sp in obj_data.splines:
                sp.use_endpoint_u = True

            bpy.ops.object.mode_set(mode='OBJECT')

        if self.bool_pipe_curve == False:
            bpy.ops.object.convert(target='MESH')
            if fill:
                self.make_pipe_filling(new_obj)
            ctx = bpy.context.copy()
            ctx['active_object'] = current
            ctx['selected_editable_objects'] = [current, new_obj]
            bpy.ops.object.join(ctx)
            current.select_set(True)
            context.view_layer.objects.active = current
            bpy.ops.object.mode_set(mode='EDIT')
        else:
            obj_data.use_fill_caps = True
            new_obj.select_set(False)
            current.select_set(True)
            context.view_layer.objects.active = current
            bpy.ops.object.mode_set(mode='EDIT')


    def paste_commit(self):
        offset = self.calc_paste_offset()

        for item in ModOperator.clipboard:
            loop = item.loop
            loop2 = []
            for p1 in loop:
                p3 = p1.copy()
                p3.co += offset
                loop2.append(p3)
            item2 = Shape()
            item2.loop = loop2
            item2.select = True
            item2.center = item.center + offset
            item2.construction = item.construction
            item2.line_only = item.line_only
            self.loops.append(item2)

    def pasting_control(self, context, event):
        if event.type == 'LEFTMOUSE':
            if event.value == 'PRESS':
                self.paste_commit()
                self.cursor = None
                self.paste = None
                self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'MOUSEMOVE':
            if self.paste != None:
                self.snap_to_grid(context, event)
                # loc = self.calc_pos(context, event)
                loc = self.cursor
                base, _ = self.paste
                self.paste = (base, loc)
                return {'RUNNING_MODAL'}

        elif event.type == 'ESC':
            if event.value == 'PRESS':
                self.paste = None
                self.handlers.pop()
                return {'RUNNING_MODAL'}

        elif event.type == 'U':
            if self.snap_disable == False:
                self.snap_disable = True
            else:
                self.snap_disable = False
            return {'RUNNING_MODAL'}

    def is_connected(self, a, b):
        item1, id1 = a
        item2, id2 = b
        if item1 == item2:
            first = min(id1, id2)
            second = max(id1, id2)
            total = len(item1.loop)
            if first == 0 and second == total - 1 and total > 2:
                if item1.line_only:
                    return False
                else:
                    return True
            else:
                if first + 1 == second:
                    return True
                else:
                    return False
        else:
            return False

    def bevel_control(self, context, event):
        obj = context.edit_object
        if event.type == 'LEFTMOUSE':
            if event.value == 'PRESS':
                self.bevel = None
                self.loop_backup_clear()
                self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'MOUSEMOVE':
            loc = self.calc_pos(context, event)
            self.set_bevel_move(loc)
            return {'RUNNING_MODAL'}

        elif event.type == 'ESC':
            if event.value == 'PRESS':
                self.bevel = None
                self.loop_backup_load()
                self.loop_backup_clear()
                self.handlers.pop()
                return {'RUNNING_MODAL'}

        elif event.type == 'WHEELUPMOUSE' or (event.type == 'RIGHT_BRACKET' and event.value == 'PRESS') \
                or (event.type == 'EIGHT' and event.value == 'PRESS'):
            self.set_bevel_cut(1)
            return {'RUNNING_MODAL'}

        elif event.type == 'WHEELDOWNMOUSE' or (event.type == 'LEFT_BRACKET' and event.value == 'PRESS') \
                or (event.type == 'SEVEN' and event.value == 'PRESS'):
            self.set_bevel_cut(-1)
            return {'RUNNING_MODAL'}

    def set_bevel_cut(self, cut):
        if self.bevel == None:
            return
        else:
            loc1, loc2, cut2 = self.bevel
            cut2 = cut2 + cut
            if cut2 < 1:
                cut2 = 1
            self.bevel = (loc1, loc2, cut2)
        self.update_bevel()

    def set_bevel_move(self, loc):
        if self.bevel == None:
            self.bevel = (loc, loc, 1)
        else:
            loc1, loc2, cut = self.bevel
            self.bevel = (loc1, loc, cut)
        self.update_bevel()

    def update_bevel(self):
        loc1, loc2, cut = self.bevel
        glen = self.grid_len
        dis = (loc2 - loc1).length * glen * 5
        if dis == 0:
            return
        sel = [e for e in self.loops if e.select]
        for item in sel:
            buf = []
            loop = item.loop
            select = item.select
            backup = item.backup
            if backup == None:
                item.backup = loop
                backup = loop
            p = backup
            ret = geo.process_bevel(p, item.line_only, dis, cut)
            if ret == None:
                continue
            item.loop = ret

    def get_matrix(self):
        if self.space == None or self.space_center == None:
            return None
        m1, m2 = self.space
        m1 = m1.copy().normalized()
        m2 = m2.copy().normalized()
        z = m1.cross(m2)
        m = Matrix.Identity(4)
        cen = self.space_center
        m[0][0:3] = m1
        m[1][0:3] = m2
        m[2][0:3] = z
        m[3][0:3] = cen
        m = m.transposed()
        return m

    def get_bounding(self, loop, mat):
        ps = [mat @ p1.co for p1 in loop]
        xs = [p.x for p in ps]
        ys = [p.y for p in ps]
        min_x = min(xs)
        max_x = max(xs)
        min_y = min(ys)
        max_y = max(ys)
        return (min_x, max_x, min_y, max_y)

    def contain_point(self, bound, loc2):
        x = loc2.x
        y = loc2.y
        x1, x2, y1, y2 = bound
        return x >= x1 and x <= x2 and y >= y1 and y <= y2

    def select_contain_shape(self, loc):
        mat = self.get_matrix().inverted()
        if mat == None:
            return
        loc2 = mat @ loc
        for item in self.loops:
            loop = item.loop
            bound = self.get_bounding(loop, mat)
            if self.contain_point(bound, loc2):
                item.select = True
                return

    def edit_select_shapes(self, loc, loc2, context):
        sel = []
        mat = self.get_matrix().inverted()
        vloc = mat @ loc
        vloc2 = mat @ loc2
        x1 = min(vloc.x, vloc2.x)
        x2 = max(vloc.x, vloc2.x)
        y1 = min(vloc.y, vloc2.y)
        y2 = max(vloc.y, vloc2.y)
        for item in self.loops:
            loop = item.loop
            for i, pt in enumerate(loop):
                p1 = pt.co
                v1 = mat @ p1
                if v1.x > x1 and v1.x < x2:
                    if v1.y > y1 and v1.y < y2:
                        self.editshape.append((item, i))

    def edit_select_drag_control(self, context, event):
        if event.type == 'LEFTMOUSE':
            if event.value == 'RELEASE':
                loc, loc2 = self.selection
                if event.shift == False:
                    self.editshape = []
                if loc == None or loc2 == None:
                    pass
                else:
                    self.edit_select_shapes(loc, loc2, context)

                self.selection = None
                self.yellow_rect = []
                self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'MOUSEMOVE':
            if self.selection != None:
                loc2 = self.calc_pos(context, event)
                loc, _ = self.selection
                self.selection = (loc, loc2)
                self.create_selection_rect(context)
                return {'RUNNING_MODAL'}

    def select_drag_control(self, context, event):
        if event.type == 'LEFTMOUSE':
            if event.value == 'RELEASE':
                loc, loc2 = self.selection
                if event.shift == False:
                    self.deselect_all()
                if loc == None or loc2 == None:
                    loc2 = self.calc_pos(context, event)
                    self.select_contain_shape(loc2)
                else:
                    self.select_shapes(context)

                self.selection = None
                self.yellow_rect = []
                self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'MOUSEMOVE':
            if self.selection != None:
                loc2 = self.calc_pos(context, event)
                loc, _ = self.selection
                self.selection = (loc, loc2)
                self.create_selection_rect(context)
                return {'RUNNING_MODAL'}

    def new_center_commit(self, loc):
        for item in self.loops:
            cen = item.center
            if item.select:
                if cen == None:
                    continue
                item.center = loc

    def new_center_control(self, context, event):
        if event.type == 'LEFTMOUSE':
            if event.value == 'PRESS':
                self.snap_to_grid(context, event)
                self.new_center_commit(self.cursor)
                self.cursor = None
                self.guide_points = None
                self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'MOUSEMOVE':
            self.snap_to_grid(context, event)
            return {'RUNNING_MODAL'}

        if event.value == 'PRESS':
            if event.type == 'ESC':
                self.cursor = None
                self.handlers.pop()
                return {'RUNNING_MODAL'}

            elif event.type == 'U':
                if self.snap_disable == False:
                    self.snap_disable = True
                else:
                    self.snap_disable = False
                return {'RUNNING_MODAL'}

    def get_numbers(self, event):
        if event.type == 'ONE':
            self.text_input += '1'
            return True
        elif event.type == 'TWO':
            self.text_input += '2'
            return True
        elif event.type == 'THREE':
            self.text_input += '3'
            return True
        elif event.type == 'FOUR':
            self.text_input += '4'
            return True
        elif event.type == 'FIVE':
            self.text_input += '5'
            return True
        elif event.type == 'SIX':
            self.text_input += '6'
            return True
        elif event.type == 'SEVEN':
            self.text_input += '7'
            return True
        elif event.type == 'EIGHT':
            self.text_input += '8'
            return True
        elif event.type == 'NINE':
            self.text_input += '9'
            return True
        elif event.type == 'ZERO':
            self.text_input += '0'
            return True
        elif event.type == 'PERIOD':
            self.text_input += '.'
            return True
        elif event.type == 'MINUS':
            self.text_input += '-'
            return True
        elif event.type == 'BACK_SPACE':
            self.text_input = self.text_input[:-1]
            return True
        else:
            return False

    def get_text_float(self, old):
        try:
            return float(self.text_input)
        except ValueError:
            return old

    def rotation_control(self, context, event):

        if event.type == 'MOUSEMOVE':
            loc = self.calc_pos(context, event)
            self.rotation_by_mouse(loc)
            return {'RUNNING_MODAL'}

        if event.value != 'PRESS':
            return

        if self.get_numbers(event) == True:
            self.rotate_deg = self.get_text_float(self.rotate_deg)
            self.loop_backup_load()
            self.rotate_selected_round_free()
            return {'RUNNING_MODAL'}

        if event.type == 'R' or event.type == 'LEFTMOUSE':
            self.rotation = None
            self.loop_backup_clear()
            self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'ESC':
            self.loop_backup_load()
            self.rotation = None
            self.loop_backup_clear()
            self.handlers.pop()
            return {'RUNNING_MODAL'}

    def scale_control(self, context, event):

        if event.type == 'MOUSEMOVE':
            loc = self.calc_pos(context, event)
            self.scale_by_mouse(loc)
            return {'RUNNING_MODAL'}

        if event.value != 'PRESS':
            return

        if event.type == 'LEFTMOUSE':
            self.scale = None
            self.loop_backup_clear()
            self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'ESC':
            self.loop_backup_load()
            self.scale = None
            self.loop_backup_clear()
            self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'X':
            self.scale_mode = 'X'
            loc = self.calc_pos(context, event)
            self.scale_by_mouse(loc)
            return {'RUNNING_MODAL'}

        elif event.type == 'Y':
            self.scale_mode = 'Y'
            loc = self.calc_pos(context, event)
            self.scale_by_mouse(loc)
            return {'RUNNING_MODAL'}

    def arraycopy_control(self, context, event):

        if event.type == 'MOUSEMOVE':
            if self.paste != None:
                self.snap_to_grid(context, event)
                loc = self.cursor
                base, _ = self.paste
                self.paste = (base, loc)
                self.array_copy_update()
                return {'RUNNING_MODAL'}

        if event.value != 'PRESS':
            return

        if event.type == 'LEFTMOUSE':
            self.paste_commit()
            self.cursor = None
            self.paste = None
            self.array_count = 0
            self.array = False
            self.array_circle = False
            self.handlers.pop()
            ModOperator.clipboard = None
            return {'RUNNING_MODAL'}

        elif event.type == 'ESC':
            self.cursor = None
            self.paste = None
            self.array_count = 0
            self.array = False
            self.handlers.pop()
            ModOperator.clipboard = None
            return {'RUNNING_MODAL'}

        elif event.type == 'WHEELUPMOUSE' or (event.type == 'RIGHT_BRACKET' and event.value == 'PRESS') \
                or (event.type == 'EIGHT' and event.value == 'PRESS'):
            self.set_array_count(1)
            self.array_copy_update()
            return {'RUNNING_MODAL'}

        elif event.type == 'WHEELDOWNMOUSE' or (event.type == 'LEFT_BRACKET' and event.value == 'PRESS') \
                or (event.type == 'SEVEN' and event.value == 'PRESS'):
            self.set_array_count(-1)
            self.array_copy_update()
            return {'RUNNING_MODAL'}

        elif event.type == 'C':
            if self.array_circle:
                self.array_circle = False
            else:
                self.array_circle = True
                if self.array_count < 2:
                    self.array_count = 2
            self.array_copy_update()
            return {'RUNNING_MODAL'}

    def array_copy_update(self):
        if self.paste == None:
            return
        cen, loc = self.paste
        main = loc - cen
        ModOperator.clipboard = []
        sel = []
        for item in self.loops:
            if item.select:
                sel.append(item)

        cut = self.array_count + 1

        for i in range(cut):
            k = i + 1
            for item in sel:
                item2 = Shape()
                item2.select = item.select
                item2.construction = item.construction
                item2.line_only = item.line_only
                loop2 = []

                if self.array_circle:
                    if i == cut - 1:
                        continue
                    deg = math.radians(360) / cut
                    d2 = deg * k
                    rot = Quaternion(self.source.normal, d2)
                    c2 = item.center - cen
                    c2.rotate(rot)
                    c2 += cen
                    item2.center = c2
                    for p1 in item.loop:
                        a1 = p1.copy()
                        a1.co -= cen
                        a1.co.rotate(rot)
                        a1.co += cen
                        loop2.append(a1)
                else:
                    m2 = main / cut
                    item2.center = item.center + (m2 * k)
                    for p1 in item.loop:
                        a1 = p1.copy()
                        a1.co += (m2 * k)
                        loop2.append(a1)

                item2.loop = loop2
                ModOperator.clipboard.append(item2)

    def get_selected(self):
        if not self.loops:
            return None

        for item in self.loops:
            if item.select:
                return item
        return None

    def get_all_selected(self):
        if self.loops == None:
            return None
        sel = []
        for item in self.loops:
            if item.select:
                sel.append(item)
        if len(sel) == 0:
            return None
        else:
            return sel

    def array_pending(self, loc):
        sel = self.get_all_selected()
        ModOperator.clipboard = sel
        cen = self.calc_paste_cen()
        self.paste = (cen, loc)

    def set_array_count(self, num):
        self.array_count = self.array_count + num
        if self.array_count < 0:
            self.array_count = 0

    def scale_by_mouse(self, loc):
        if self.scale == None:
            return
        loc1, _ = self.scale
        if loc1 == None:
            self.scale = (loc, None)
            return
        else:
            self.scale = (loc1, loc)
            item = self.get_selected()
            if item != None:
                cen = item.center
                self.scale_free(cen, loc1, loc)

    def get_aligned_axis(self):
        main_hori, main_vert = self.space_point
        h1, h2 = main_hori
        v1, v2 = main_vert
        haxis = h2 - h1
        vaxis = v2 - v1
        obj = bpy.context.object
        px = obj.matrix_world.inverted() @ Vector((1, 0, 0))
        if haxis.angle(px) < vaxis.angle(px):
            return haxis, vaxis
        else:
            return vaxis, haxis

    def scale_free(self, cen, loc1, loc2):
        vcen = self.vm_convert(cen)
        vloc1 = self.vm_convert(loc1)
        vloc2 = self.vm_convert(loc2)
        p1 = vloc1 - vcen
        p2 = vloc2 - vcen
        scale = p2.length / p1.length

        haxis, vaxis = self.get_aligned_axis()

        self.loop_backup_load()
        for item in self.loops:
            if item.select:
                loop2 = []
                for p1 in item.loop:
                    s1 = p1.copy()
                    if self.scale_mode == None:
                        s1.co = (s1.co - cen) * scale + cen
                    else:
                        s2 = s1.co - cen
                        s2_h = s2.project(haxis)
                        s2_v = s2.project(vaxis)
                        if self.scale_mode == 'X':
                            s2_h *= scale
                        elif self.scale_mode == 'Y':
                            s2_v *= scale
                        s1.co = cen + s2_h + s2_v

                    loop2.append(s1)
                item.loop = loop2

    def rotation_by_mouse(self, loc):
        if self.rotation == None:
            return
        loc1, _ = self.rotation
        if loc1 == None:
            self.rotation = (loc, None)
            return
        else:
            self.rotation = (loc1, loc)
            item = self.get_selected()
            if item != None:
                cen = item.center
                self.rotate_free(cen, loc1, loc)

    def rotate_free(self, cen, loc1, loc2):
        vcen = self.vm_convert(cen)
        vloc1 = self.vm_convert(loc1)
        vloc2 = self.vm_convert(loc2)
        p1 = vloc1 - vcen
        p2 = vloc2 - vcen
        deg = math.atan2(p1.y, p1.x) - math.atan2(p2.y, p2.x)
        deg = math.degrees(deg)
        deg = deg - (deg % 5)
        self.text_input = str(deg)
        self.rotate_deg = deg
        self.loop_backup_load()
        self.rotate_selected_round_free()

    def vm_convert(self, loc):
        obj = bpy.context.edit_object
        world = obj.matrix_world
        vm = self.get_view_matrix()
        loc2 = world @ loc
        loc3 = vm @ loc2
        return loc3

    def loop_backup_load(self):
        sel = [e for e in self.loops if e.select]
        for item in sel:
            item.loop = item.backup

    def loop_backup_clear(self):
        sel = [e for e in self.loops if e.select]
        for item in sel:
            item.backup = None

    def loop_backup(self):
        sel = [e for e in self.loops if e.select]
        for item in sel:
            item.backup = item.loop

    def points_backup_load(self):
        sel = [e for e in self.loops if e.select]
        for item in sel:
            item.loop = []
            item.line_only = item.line_only2
            for p1 in item.backup:
                item.loop.append(p1.copy())

    def points_backup_clear(self):
        sel = [e for e in self.loops if e.select]
        for item in sel:
            item.backup = None

    def points_backup(self):
        sel = [e for e in self.loops if e.select]
        for item in sel:
            item.backup = []
            item.line_only2 = item.line_only
            for p1 in item.loop:
                item.backup.append(p1.copy())

    def snap_only_guide_point(self, loc):
        limit = self.grid_len / 2
        loc2 = loc.copy()
        if self.guide_points != None:
            for p1 in self.guide_points:
                if (p1 - loc).length < limit:
                    loc2 = p1
        return loc2

    def snap_to_grid(self, context, event):
        if self.snap_disable:
            loc = self.calc_pos(context, event)
            loc2 = self.snap_only_guide_point(loc)
            self.cursor = loc2

        else:
            loc = self.calc_pos(context, event)
            loc2 = self.snap(context, loc)
            self.cursor = loc2

    def selection_control(self, context, event):
        obj = context.edit_object

        if event.value != 'PRESS':
            return

        if event.type == 'RIGHTMOUSE' or event.type == 'W':
            self.text_handlers.pop()
            self.handlers.pop()
            self.selection_mode = False
            self.snap_to_grid(context, event)
            self.calc_intersection()
            return {'RUNNING_MODAL'}

        elif event.type == 'LEFTMOUSE' and event.alt == False:
            if event.ctrl:
                loc = self.calc_pos(context, event)
                p, item = plane.get_nearby(loc, [self.source] + self.other_source,
                                           self.loops)
                self.main_edge = p
                self.main_edge_item = item
                self.update_guide(context)
                return {'RUNNING_MODAL'}
            else:
                loc = self.calc_pos(context, event)
                self.selection = (loc, None)
                self.handlers.append(self.select_drag_control)
                return {'RUNNING_MODAL'}

        elif event.type == 'RET' or event.type == 'Q':
            if self.drawing:
                return {'RUNNING_MODAL'}
            else:
                self.finish_action(context)
                self.clean_up(context)
                return {'FINISHED'}

        elif event.type == 'ESC':
            self.text_handlers.pop()
            self.handlers.pop()
            self.selection_mode = False
            return {'RUNNING_MODAL'}

        elif event.type == 'A':
            for item in self.loops:
                item.select = True
            return {'RUNNING_MODAL'}

        elif event.type == 'U':
            if self.snap_disable == False:
                self.snap_disable = True
            else:
                self.snap_disable = False
            return {'RUNNING_MODAL'}

        elif event.type == 'ONE':
            self.operation_mode = 'ngon'
            self.finish_action(context)
            self.clean_up(context)
            return {'FINISHED'}

        elif event.type == 'TWO':
            self.operation_mode = 'newface'
            self.finish_action(context)
            self.clean_up(context)
            return {'FINISHED'}

        elif event.type == 'THREE':
            self.operation_mode = 'boolcut'
            self.finish_action(context)
            self.clean_up(context)
            return {'FINISHED'}

        elif event.type == 'FOUR':
            self.operation_mode = 'linepipe'
            self.finish_action(context)
            self.clean_up(context)
            return {'FINISHED'}


        elif event.type == 'V' and event.ctrl:
            loc = self.calc_pos(context, event)
            self.deselect_all()
            if self.pending_paste(loc):
                self.handlers.append(self.pasting_control)
            return {'RUNNING_MODAL'}

        if self.get_all_selected() == None:
            return
        else:
            if event.type == 'B':
                self.handlers.append(self.bevel_control)
                return {'RUNNING_MODAL'}

            elif event.type == 'C' and event.ctrl:
                self.copy_loop()
                return {'RUNNING_MODAL'}

            elif event.type == 'X' and event.ctrl == False:
                self.guide_points = [p1.co for item in self.loops for p1 in item.loop]
                self.snap_to_grid(context, event)
                self.handlers.append(self.new_center_control)
                return {'RUNNING_MODAL'}

            elif event.type == 'X' and event.ctrl:
                self.copy_loop()
                if ModOperator.clipboard:
                    self.delete_shape()
                return {'RUNNING_MODAL'}

            elif event.type == 'G':
                self.copy_loop()
                if ModOperator.clipboard:
                    self.delete_shape()
                    loc = self.calc_pos(context, event)
                    self.pending_paste(loc)
                    self.handlers.append(self.pasting_control)
                return {'RUNNING_MODAL'}

            elif event.type == 'M':
                if event.ctrl:
                    self.symmetric()
                else:
                    self.rotate_selected_hori()
                return {'RUNNING_MODAL'}


            elif event.type == 'N':
                self.rotate_selected_vert()
                return {'RUNNING_MODAL'}

            elif event.type == 'R':
                self.text_input = ''
                self.rotation = (None, None)
                self.loop_backup()
                self.handlers.append(self.rotation_control)
                return {'RUNNING_MODAL'}

            elif event.type == 'S':
                self.scale = (None, None)
                self.scale_mode = None
                self.loop_backup()
                self.handlers.append(self.scale_control)
                return {'RUNNING_MODAL'}

            elif event.type == 'D':
                self.array = True
                loc = self.calc_pos(context, event)
                self.array_pending(loc)
                self.handlers.append(self.arraycopy_control)
                return {'RUNNING_MODAL'}

            elif event.type == 'T':
                self.rotate_selected_round()
                return {'RUNNING_MODAL'}

            elif event.type == 'I':
                loc = self.calc_pos(context, event)
                self.base_location = loc
                self.points_backup()
                self.handlers.append(self.inset_control)
                return {'RUNNING_MODAL'}

            elif event.type == 'DEL' or event.type == 'BACK_SPACE':
                self.delete_shape()
                return {'RUNNING_MODAL'}

    def inset_shape(self, loc1, loc2, context):
        loc2b = self.vm_convert(loc1)
        loc1b = self.vm_convert(loc2)
        dis = loc2b.x - loc1b.x
        dis *= -1

        sel = [item for item in self.loops if item.select]
        self.clocktise(sel)
        sn = self.source.normal.copy()

        for item in sel:
            item.integrity_check_length()
            item.integrity_check_colinear()
            item.link_point_edges(sn)
            for p in item.loop:
                p.move = Vector()
            for p in item.loop:
                if p.back != None and p.next != None:
                    s1 = self.calc_sliding(p.back.inward, dis, p, p.back, p.next)
                    s2 = self.calc_sliding(p.inward, dis, p, p.next, p.back)
                    p.move = s1 + s2

            if item.line_only == False:
                for p in item.loop:
                    p.co += p.move
            else:
                p = item.loop[0]
                p.move = p.inward * dis
                p = item.loop[-1]
                p.move = p.inward * dis
                loop1 = []
                loop2 = []
                for p in item.loop:
                    p2 = p.copy()
                    p2.co += p.move
                    loop1.append(p2)
                for p in item.loop:
                    p2 = p.copy()
                    p2.co -= p.move
                    loop2.append(p2)
                loop2 = list(reversed(loop2))
                item.loop = loop1 + loop2
                item.line_only = False
                item.integrity_check_length()
                item.integrity_check_colinear()

    def calc_sliding(self, arm, dis, p, pback, pnext):
        m1 = pback.co - p.co
        m2 = pnext.co - p.co
        m2p = m2.project(arm)
        if m2p.length == 0:
            return None
        slide = (dis / m2p.length) * m2.length
        s1 = m2.normalized() * slide
        if p.small:
            return s1
        else:
            return s1 * -1

    def inset_control(self, context, event):
        obj = context.edit_object
        if event.type == 'LEFTMOUSE':
            if event.value == 'PRESS':
                self.base_location = None
                self.points_backup_clear()
                self.handlers.pop()
            return {'RUNNING_MODAL'}

        elif event.type == 'MOUSEMOVE':
            self.points_backup_load()
            loc = self.calc_pos(context, event)
            self.inset_shape(self.base_location, loc, context)
            return {'RUNNING_MODAL'}

        elif event.type == 'ESC':
            if event.value == 'PRESS':
                self.base_location = None
                self.points_backup_load()
                self.points_backup_clear()
                self.handlers.pop()
                return {'RUNNING_MODAL'}

    def symmetric(self):
        for item in self.loops:
            if item.line_only and len(item.loop) >= 3:
                head = item.loop[0]
                last = item.loop[-1]
                main = last.co - head.co
                qua = Quaternion(main, math.radians(180))
                added = []
                for p in item.loop:
                    if p != head and p != last:
                        p2 = p.copy()
                        p2.co -= head.co
                        p2.co.rotate(qua)
                        p2.co += head.co
                        added.append(p2)
                added = reversed(added)
                item.loop += added
                item.line_only = False

    def rotate_selected_hori(self):
        if self.space == None:
            return
        main, vertical = self.space
        rot = Quaternion(vertical, math.radians(180.0))
        self.rotate_loops(rot)

    def rotate_selected_vert(self):
        if self.space == None:
            return
        main, vertical = self.space
        rot = Quaternion(main, math.radians(180.0))
        self.rotate_loops(rot)

    def rotate_selected_round(self):
        rot = Quaternion(self.source.normal, math.radians(-90.0))
        self.rotate_loops(rot)

    def rotate_selected_round_free(self):
        if self.rotate_deg == None:
            return
        deg = self.rotate_deg
        rot = Quaternion(self.source.normal, math.radians(deg * -1))
        self.rotate_loops(rot)

    def rotate_loops(self, rot):
        selected = None
        for item in self.loops:
            loop = item.loop
            if item.select:
                selected = item
                break

        if selected == None:
            return
        cen = selected.center

        for item in self.loops:
            loop = item.loop
            if item.select:
                loop2 = self.rotate_loop_single(rot, loop, cen)
                item.loop = loop2

    def rotate_loop_single(self, rot, loop, cen):
        loop2 = []
        for p1 in loop:
            a1 = p1.copy()
            a1.co = a1.co - cen
            a1.co.rotate(rot)
            a1.co = a1.co + cen
            loop2.append(a1)
        return loop2

    def copy_loop(self):
        ModOperator.clipboard = [e.copy() for e in self.loops if e.select]
        ModOperator.clipboard_sn = self.source.normal.copy()
        main_hori, main_vert = self.space_point
        m1, m2 = main_hori
        ModOperator.clipboard_main = (m2 - m1).normalized()

    def calc_paste_rotation(self):
        cn = None
        if ModOperator.clipboard_sn == None:
            cn = self.source.normal.copy()
        else:
            cn = ModOperator.clipboard_sn
        rot = cn.rotation_difference(self.source.normal)
        return rot

    def calc_paste_cen(self):
        item = ModOperator.clipboard[0]
        return item.center

    def pending_paste(self, loc):
        if ModOperator.clipboard == None or (not ModOperator.clipboard):
            return False

        cen = self.calc_paste_cen()
        cn = ModOperator.clipboard_sn

        if cn != None and geo.same_direction(cn, self.source.normal) == False:
            rot = self.calc_paste_rotation()
            for item in ModOperator.clipboard:
                loop = item.loop
                item.loop = self.rotate_loop_single(rot, loop, cen)
            ModOperator.clipboard_sn = self.source.normal.copy()

            if ModOperator.clipboard_main != None:
                mat = rot.to_matrix().to_4x4()
                ModOperator.clipboard_main = mat @ ModOperator.clipboard_main

        main_hori, _ = self.space_point
        m1, m2 = main_hori
        main2 = (m2 - m1).normalized()
        gmain = ModOperator.clipboard_main
        if gmain != None and geo.same_direction(gmain, main2) == False:
            rot = gmain.rotation_difference(main2)
            for item in ModOperator.clipboard:
                loop = item.loop
                item.loop = self.rotate_loop_single(rot, loop, cen)
            ModOperator.clipboard_main = main2

        self.paste = (cen, loc)
        return True

    def delete_shape(self):
        while True:
            finish = True
            for i, item in enumerate(self.loops):
                if item.select:
                    del self.loops[i]
                    finish = False
                    break
            if finish:
                break

    def is_between(self, p1, p2, k):
        e = (p2 - p1).length - ((k - p1).length + (p2 - k).length)
        return abs(e) < 0.05

    def deselect_all(self):
        for item in self.loops:
            item.select = False

    def select_shapes(self, context):
        rect = self.select_area
        if rect == None or len(rect) != 4:
            return None
        a1, a2, a3, a4 = rect
        e1 = a2 - a1
        e2 = a4 - a1
        for item in self.loops:
            loop = item.loop
            pos = [p1.co for p1 in loop]
            total = Vector()
            for p in pos:
                total += p
            cen = total / len(pos)
            cen2 = cen - a1
            p1 = cen2.project(e1)
            p2 = cen2.project(e2)
            if self.is_between(Vector(), e1, p1) and self.is_between(Vector(), e2, p2):
                item.select = True

    def create_selection_rect(self, context):
        if self.selection == None:
            return
        loc1, loc2 = self.selection
        self.yellow_rect = []
        obj = context.edit_object
        world = obj.matrix_world

        if self.space_point == None:
            return
        main_hori, main_vert = self.space_point
        m1, m2 = main_hori

        main = m2 - m1
        p1 = (loc1 - m1).project(main) + m1
        p2 = (loc2 - m1).project(main) + m1
        off1 = loc1 - p1
        off2 = loc2 - p2
        fp1 = p2 + off1
        fp2 = p1 + off2
        w1 = loc1
        w2 = fp1
        w3 = loc2
        w4 = fp2
        self.select_area = [w1, w2, w3, w4]
        self.yellow_rect = [w1, w2, w2, w3, w3, w4, w4, w1]

    def pencil_draw(self, context):
        cursor = self.cursor

        if self.item == None:
            return
        if (cursor - self.firstv).length < self.grid_len / 5:
            if len(self.currentloop) > 3:
                self.close_edge(context)
                return

        draw = False
        if self.lastv == None:
            draw = True
        else:
            p = self.lastv - cursor
            if p.length > self.grid_len / 2:
                draw = True
        if draw:
            self.insert_point(context, cursor, True)

    def add_current_points(self):
        self.guide_points = []
        for item in self.loops:
            for p in item.loop:
                self.guide_points.append(p.co)

    def overall(self, context, event):

        if event.type == 'LEFTMOUSE' and event.alt == False:
            if event.value == 'PRESS':
                if event.ctrl:
                    loc = self.calc_pos(context, event)
                    p, item = plane.get_nearby(loc, [self.source] + self.other_source,
                                               self.loops)
                    self.main_edge = p
                    self.main_edge_item = item
                    self.update_guide(context)
                    return {'RUNNING_MODAL'}
                elif event.shift:
                    loc = self.calc_pos(context, event)
                    v = plane.get_nearby_vert(loc, [self.source] + self.other_source,
                                              self.loops)
                    self.shift_vert = v
                    self.update_guide(context)
                    return {'RUNNING_MODAL'}

                else:

                    if self.circle != None:
                        self.draw_circle_point(context, self.cursor)
                    elif self.con_mode != None:
                        self.add_con_point(context, self.cursor)
                        return {'RUNNING_MODAL'}
                    else:
                        self.draw_down = True
                        self.insert_point(context, self.cursor, False)
                    return {'RUNNING_MODAL'}
            elif event.value == 'RELEASE':
                self.draw_down = False
                return {'RUNNING_MODAL'}

        elif event.type == 'MOUSEMOVE':

            if self.draw_down and self.snap_disable:
                self.snap_to_grid(context, event)
                self.pencil_draw(context)
                return {'RUNNING_MODAL'}
            else:
                if self.snap_to_shape:
                    self.add_current_points()
                else:
                    self.guide_points = []
                self.snap_to_grid(context, event)
                self.move_edge(context, self.cursor)
                if self.circle != None:
                    p1, p2, _ = self.circle
                    if p1 != None:
                        self.circle = (p1, self.cursor, None)
                        self.draw_circle(context)
                return {'RUNNING_MODAL'}

        elif event.type == 'WHEELUPMOUSE':
            return self.basic_scroll(context, event, 1)

        elif event.type == 'WHEELDOWNMOUSE':
            return self.basic_scroll(context, event, -1)

        if event.value == 'PRESS':
            if event.type == 'ESC':
                if self.circle != None:
                    self.circle = None
                    return {'RUNNING_MODAL'}
                elif self.item != None:
                    self.reset_draw()
                    return {'RUNNING_MODAL'}
                else:
                    self.clean_up(context)
                    return {'CANCELLED'}

            elif event.type == 'RET' or event.type == 'Q':
                if self.drawing:
                    return {'RUNNING_MODAL'}
                else:
                    self.finish_action(context)
                    self.clean_up(context)
                    return {'FINISHED'}

            elif event.type == 'C':
                if self.circle == None:
                    self.circle = (None, None, None)
                    self.reset_draw()
                else:
                    self.circle = None
                return {'RUNNING_MODAL'}

            elif event.type == 'G':
                cen = self.source.calc_center_median()
                m1, v1 = self.space
                m2 = m1 / 2
                self.main_edge = (cen + m2, cen - m2)
                self.update_guide(context)
                return {'RUNNING_MODAL'}

            elif event.type == 'U':
                if event.ctrl:
                    if self.snap_to_shape == False:
                        self.snap_to_shape = True
                    else:
                        self.snap_to_shape = False
                else:
                    if self.snap_disable == False:
                        self.snap_disable = True
                    else:
                        self.snap_disable = False
                return {'RUNNING_MODAL'}

            elif event.type == 'ONE':
                self.operation_mode = 'ngon'
                self.finish_action(context)
                self.clean_up(context)
                return {'FINISHED'}

            elif event.type == 'TWO':
                self.operation_mode = 'lineart'
                self.finish_action(context)
                self.clean_up(context)
                return {'FINISHED'}

            elif event.type == 'THREE':
                self.operation_mode = 'newface'
                self.finish_action(context)
                self.clean_up(context)
                return {'FINISHED'}

            elif event.type == 'FOUR':
                self.operation_mode = 'boolcut'
                self.finish_action(context)
                self.clean_up(context)
                return {'FINISHED'}

            elif event.type == 'NINE':
                self.operation_mode = 'linepipe'
                self.finish_action(context)
                self.clean_up(context)
                return {'FINISHED'}

            elif event.type == 'SPACE':
                self.end_without_close_edge(context)
                return {'RUNNING_MODAL'}

    def basic_scroll(self, context, event, num):
        if event.ctrl:
            if ModOperator.mode_A:
                ModOperator.size += num
                if ModOperator.size < 2:
                    ModOperator.size = 2
            else:
                ModOperator.size_rel += num
                if ModOperator.size_rel < 2:
                    ModOperator.size_rel = 2

            self.update_guide(context)
            return {'RUNNING_MODAL'}
        elif event.shift:
            ModOperator.circle_cut += num
            if ModOperator.circle_cut < 3:
                ModOperator.circle_cut = 3
            self.draw_circle(context)
            return {'RUNNING_MODAL'}
        elif event.alt:
            if self.mode_A:
                self.scale_up += num
                if self.scale_up < 0:
                    self.scale_up = 0
            else:
                self.scale_up_rel += num
                if self.scale_up_rel < 0:
                    self.scale_up_rel = 0
            self.update_guide(context)
            return {'RUNNING_MODAL'}

    def add_con_point(self, context, cursor):
        if self.con_mode == None:
            return

        s1, s2 = self.con_mode
        self.insert_point(context, self.cursor, False)
        if s1 != None:
            self.add_conline(context)
            self.con_mode = (None, None)
        else:
            self.con_mode = (self.cursor, None)

    def reset_draw(self):
        self.lastv = None
        self.item = None
        self.drawing = False
        self.currentloop = []

    def modal(self, context, event):
        editing = context.mode == 'EDIT_MESH'
        if editing == False:
            self.clean_up(context)
            return {'CANCELLED'}

        context.area.tag_redraw()
        obj = context.edit_object

        if self.show_keys:
            self.update_keys(context, event)

        wheel = ('WHEEL' in event.type)
        mac_key = event.type == 'SEVEN' or event.type == 'EIGHT'
        mac_key2 = event.type == 'RIGHT_BRACKET' or event.type == 'LEFT_BRACKET'
        grid_enlarge = wheel or mac_key or mac_key2

        if not self.handlers:
            return {'PASS_THROUGH'}
        elif 'NUM' in event.type:
            return {'PASS_THROUGH'}
        elif event.alt and grid_enlarge == False and (event.type != 'M' and event.type != 'N'):
            return {'PASS_THROUGH'}
        else:
            func = self.handlers[-1]
            p = func(context, event)
            if p == None:
                if 'MOUSE' in event.type:
                    if event.type == 'RIGHTMOUSE':
                        return {'RUNNING_MODAL'}
                    else:
                        return {'PASS_THROUGH'}
                else:
                    return {'RUNNING_MODAL'}
            else:
                return p

    def save_geo_point(self):
        ps = []
        es = []
        sc = [self.source] + self.other_source
        for f1 in sc:
            for v1 in f1.verts:
                ps.append(v1.co)
        self.geo_points = ps

        for f1 in sc:
            for e1 in f1.edges:
                p1 = e1.verts[0].co
                p2 = e1.verts[1].co
                es.append((p1, p2))
        self.geo_edges = es

    def invoke(self, context, event):
        if context.edit_object:

            is_vert_mode, is_edge_mode, is_face_mode = context.tool_settings.mesh_select_mode
            bm = geo.get_bm(context)

            self.virtualface = True

            if is_face_mode:
                found = self.save_source(context, bm)
                if found:
                    self.virtualface = False
                elif found == None:
                    return {'CANCELLED'}
                elif found == False:
                    found2 = self.viewplane_source(context, event)
                    if found2 == None:
                        return {'CANCELLED'}
            elif is_edge_mode:
                found = self.save_source_edge(context, bm)
                if found == None:
                    return {'CANCELLED'}

            elif is_vert_mode:
                found = self.save_source_vert(context, bm)
                if found == None:
                    return {'CANCELLED'}

            self.save_geo_point()

            self.text_size = 10
            self.operation_mode = 'boolcut'
            ModOperator.mode_A = True
            self.show_keys = False
            self.text_pos_x = 120
            self.color_grid_line = (1, 1, 1, 0.2)
            self.color_shape = (1, 1, 1, 1)

            self.handlers.append(self.overall)
            self.text_handlers.append(self.base_text)

            args = (self, context)

            if self.update_guide(context) == None:
                print('guide error')
                return {'CANCELLED'}

            ModOperator.remove_draw()
            ModOperator.handle3d = bpy.types.SpaceView3D.draw_handler_add(
                draw_3d, args, 'WINDOW', 'POST_VIEW')
            ModOperator.handle = bpy.types.SpaceView3D.draw_handler_add(
                draw_text_callback, args, 'WINDOW', 'POST_PIXEL')

            context.area.tag_redraw()

            context.window_manager.modal_handler_add(self)

            return {'RUNNING_MODAL'}
        else:
            self.report({'WARNING'}, "No active object, could not finish")
            return {'CANCELLED'}

        wm = context.window_manager
        return wm.invoke_props_dialog(self)

    def draw_text(self, pos, text):
        if pos == None:
            return
        font_id = 0
        blf.color(font_id, 1, 1, 1, 1)
        blf.position(font_id, pos[0], pos[1], 0)
        blf.size(font_id, self.text_size, 72)
        blf.draw(font_id, text)

    def base_text(self, context):
        if ModOperator.mode_A:
            mode = 'Absolute size'
        else:
            mode = 'Relative size'

        if self.operation_mode == 'lineart':
            cut = 'Line art'
        elif self.operation_mode == 'ngon':
            cut = 'Cut and fill (n-gon)'
        elif self.operation_mode == 'newface':
            cut = 'New face'
        elif self.operation_mode == 'boolcut':
            cut = 'Boolean cut'

        txtall = []

        if ModOperator.mode_A:
            size = ModOperator.size * 2
        else:
            size = ModOperator.size_rel

        txtall.append('Grids : ' + str(size))
        txtall.append('Ctrl + Mouse Wheel : Set number of grids')
        txtall.append('Alt + Mouse Wheel : Enlarge grid plane')
        txtall.append('Right-Click : Switch to Selection Mode')
        return txtall

    def calc_pos(self, context, event):
        obj = context.edit_object
        view_vector, view_point, world_loc = self.get_3d_cursor(context, event)

        scenter = self.source.calc_center_median()
        sn = self.source.normal.copy()

        world = obj.matrix_world
        world2 = world.inverted()
        viewer = world2 @ view_point
        tar = world2 @ world_loc

        w = viewer - scenter
        u = tar - viewer
        snd = sn.dot(u)

        if snd != 0:
            s1 = sn.dot(w) * -1 / snd
            s2 = viewer + u * s1
            return s2
        else:
            return scenter

    def create_cross(self, context, cursor):
        obj = context.edit_object
        world = obj.matrix_world
        coord = []
        m = Vector((0, 0, 1)).normalized().rotation_difference(self.source.normal).to_matrix().to_4x4()

        s_len = self.grid_len / 5

        c1 = (world @ cursor) + (m @ Vector((s_len, -1 * s_len, 0)))
        c2 = (world @ cursor) + (m @ Vector((-1 * s_len, s_len, 0)))
        c3 = (world @ cursor) + (m @ Vector((s_len, s_len, 0)))
        c4 = (world @ cursor) + (m @ Vector((-1 * s_len, -1 * s_len, 0)))
        coord += [c1, c2]
        coord += [c3, c4]

        return coord

    def calc_paste_offset(self):
        if self.array:
            return Vector()
        else:
            base, cursor = self.paste
            offset = cursor - base
            return offset

    def get_mode(self, context):
        is_vert_mode, is_edge_mode, is_face_mode = context.tool_settings.mesh_select_mode
        return is_vert_mode, is_edge_mode, is_face_mode

    def clean_up(self, context):
        ModOperator.remove_draw()

    def get_center(self, fs):
        total = Vector()
        for s in fs:
            total += s.calc_center_median()
        return total / len(fs)

    def update_keys(self, context, event):
        if self.keys_list == None:
            self.keys_list = []

        p = event.type
        if 'MOUSEMOVE' in p or p == 'MIDDLEMOUSE':
            return

        if 'SHIFT' in p or 'CTRL' in p or 'ALT' in p:
            return

        if event.value != 'PRESS' and event.value != 'CLICK':
            return

        if p == 'LEFTMOUSE':
            p = 'Left Mouse Click'

        if p == 'RIGHTMOUSE':
            p = 'Right Mouse Click'

        if p == 'WHEELUPMOUSE':
            p = 'Scroll Mouse Wheel (up)'

        if p == 'WHEELDOWNMOUSE':
            p = 'Scroll Mouse Wheel (down)'

        if event.ctrl:
            p = 'Ctrl ' + p

        if event.shift:
            p = 'Shift ' + p

        if event.alt:
            p = 'Alt ' + p

        if len(self.keys_list) == 0:
            self.keys_index += 1
            self.keys_list.append((self.keys_index, p, 1))
            self.keys_list = self.keys_list[-15:]
        else:
            last = self.keys_list[-1]
            index, p2, count = last
            if p == p2:
                count += 1
                self.keys_list[-1] = (index, p, count)
            else:
                self.keys_index += 1
                self.keys_list.append((self.keys_index, p, 1))
                self.keys_list = self.keys_list[-15:]


    def check_clockwise(self, item):
        sn = self.source.normal
        ps1 = item.loop
        ps2 = ps1[1:] + [ps1[0]]
        ps3 = ps2[1:] + [ps2[0]]
        total = Vector()
        for a, b, c in zip(ps1, ps2, ps3):
            m1 = a.co - b.co
            m2 = c.co - b.co
            cr1 = m1.cross(m2)
            total += cr1
        if total.length == 0:
            return True
        if total.normalized() == sn.normalized():
            return True

    def clocktise(self, loops):
        for item in loops:
            if self.check_clockwise(item):
                pass
            else:
                item.loop = list(reversed(item.loop))

    def colinear(self, p1, p2, p3, p4):
        mainlen = (p2.co - p1.co).length
        if (p3.co - p1.co).length + (p2.co - p3.co).length - mainlen < 0.001:
            return p3
        if (p4.co - p1.co).length + (p2.co - p4.co).length - mainlen < 0.001:
            return p4
        return False