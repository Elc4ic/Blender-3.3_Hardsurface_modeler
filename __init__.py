import importlib

import bpy
from . import geo
from . import lining
from . import modeler
from . import plane
from . import pref

bl_info = {
    "name": "Modeler",
    "description": "",
    "author": "Andrey",
    "version": (0, 0, 1),
    "blender": (3, 4, 0),
    "location": "",
    "category": "Mesh",
    "wiki_url": ""
}


def menu_func(self, context):
    self.layout.operator_context = "INVOKE_DEFAULT";
    self.layout.operator(modeler.ModOperator.bl_idname)


def register():
    importlib.reload(modeler)
    importlib.reload(geo)
    importlib.reload(lining)
    importlib.reload(plane)
    importlib.reload(pref)

    bpy.utils.register_class(modeler.ModOperator)
    bpy.types.VIEW3D_MT_edit_mesh_context_menu.append(menu_func)


def unregister():
    bpy.utils.unregister_class(pref.GridModelerPreferences)
    bpy.types.VIEW3D_MT_edit_mesh_context_menu.remove(menu_func)
    bpy.utils.unregister_class(modeler.ModOperator)


if __name__ == "__main__":
    register()


def test():
    print(__package__)

    try:
        unregister()
    except:
        pass

    try:
        register()
    except:
        pass

    print('test loaded')
