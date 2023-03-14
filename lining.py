import bgl
import gpu
from gpu_extras.batch import batch_for_shader

shader = gpu.shader.from_builtin('3D_UNIFORM_COLOR')


def draw_line(points, color, blend=False, smooth=False, width=1):
    draw_line_gl(points, color, blend=blend, smooth=smooth, width=width)
    return

    global shader
    gpu.state.blend_set('ALPHA')
    gpu.state.line_width_set(width)

    shader.bind()
    shader.uniform_float("color", color)
    batch = batch_for_shader(shader, 'LINES', {"pos": points})
    batch.draw(shader)


def draw_line_gl(points, color, blend=False, smooth=False, width=1):
    global shader

    if blend:
        bgl.glEnable(bgl.GL_BLEND)
    else:
        bgl.glDisable(bgl.GL_BLEND)

    if smooth:
        bgl.glEnable(bgl.GL_LINE_SMOOTH)
    else:
        bgl.glDisable(bgl.GL_LINE_SMOOTH)

    bgl.glLineWidth(width)

    shader.bind()
    shader.uniform_float("color", color)
    batch = batch_for_shader(shader, 'LINES', {"pos": points})
    batch.draw(shader)

    bgl.glDisable(bgl.GL_BLEND)
    bgl.glDisable(bgl.GL_LINE_SMOOTH)
    bgl.glLineWidth(1)


