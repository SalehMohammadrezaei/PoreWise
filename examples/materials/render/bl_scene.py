"""Blender scene builder for the CV renders (Cycles, GPU). Run:

    LD_LIBRARY_PATH=/usr/lib/wsl/lib blender -b --factory-startup -P bl_scene.py -- spec.json

spec.json:
  out            PNG path (RGBA, transparent background; compose_*.py puts it on white or dark)
  res, samples   [w, h], Cycles samples
  meshes         [{ply | stl, color [r,g,b] sRGB 0..1, rough, metal, smooth, scale, offset}]
  cylinders      {centres [[y, z]...], radius, x0, x1, color, rough}  fibres along x
  boxes          [{size [x,y,z], loc [x,y,z], rot_z deg, color, rough, metal}]   (illustrative detail)
  curves         [{npz (pts, sizes, col), radius, emit, rough}]  streamlines as hair curves
  camera         {loc, target, lens, fstop, focus}
  lights         {key, fill, rim}: strengths in W (area lights placed around target)
  ground         {z, size}  shadow-catcher plane
"""
import json, math, sys
import bpy, numpy as np
from mathutils import Vector

spec = json.load(open(sys.argv[sys.argv.index("--") + 1]))
sc = bpy.context.scene
for o in list(bpy.data.objects):
    bpy.data.objects.remove(o, do_unlink=True)

# ---- render settings
sc.render.engine = "CYCLES"
prefs = bpy.context.preferences.addons["cycles"].preferences
try:
    prefs.compute_device_type = "OPTIX"; prefs.get_devices()
    for d in prefs.devices:
        d.use = d.type in ("OPTIX", "CUDA")
    sc.cycles.device = "GPU"
except Exception as e:                                                    # noqa: BLE001
    print("GPU unavailable, CPU render:", e)
sc.cycles.samples = spec.get("samples", 256); sc.cycles.use_denoising = True
sc.cycles.max_bounces = 8; sc.cycles.transparent_max_bounces = 16
sc.render.resolution_x, sc.render.resolution_y = spec["res"]; sc.render.resolution_percentage = 100
sc.render.film_transparent = True
sc.render.image_settings.file_format = "PNG"; sc.render.image_settings.color_mode = "RGBA"; sc.render.image_settings.color_depth = "16"
sc.view_settings.view_transform = spec.get("view_transform", "AgX"); sc.view_settings.look = spec.get("look", "AgX - Base Contrast")
sc.render.filepath = spec["out"]
sc.view_settings.exposure = spec.get("exposure", 0.0)
world = bpy.data.worlds.new("w"); sc.world = world; world.use_nodes = True
bg = world.node_tree.nodes["Background"]; bg.inputs[0].default_value = (*spec.get("world", [1, 1, 1]), 1); bg.inputs[1].default_value = spec.get("world_strength", 0.5)


def srgb2lin(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def material(name, color=(0.8, 0.8, 0.8), rough=0.5, metal=0.0, attr=None, emit=0.0, coat=0.0, sss=0.0, fibre=None, glass=None, spec=None, fade=None, flat_attr=False):
    m = bpy.data.materials.new(name); m.use_nodes = True
    nt = m.node_tree; b = nt.nodes["Principled BSDF"]
    b.inputs["Roughness"].default_value = rough; b.inputs["Metallic"].default_value = metal
    if spec is not None:                                    # matte look: weaker specular layer
        b.inputs["Specular IOR Level"].default_value = spec
    if coat: b.inputs["Coat Weight"].default_value = coat
    if fibre:                                               # fine striations along the fibre direction (bump only)
        tc = nt.nodes.new("ShaderNodeTexCoord"); wv = nt.nodes.new("ShaderNodeTexWave")
        wv.wave_type = "BANDS"; wv.bands_direction = fibre.get("across", "Y"); wv.inputs["Scale"].default_value = fibre.get("scale", 400)
        wv.inputs["Distortion"].default_value = fibre.get("distortion", 2.0); wv.inputs["Detail"].default_value = 2
        bp = nt.nodes.new("ShaderNodeBump"); bp.inputs["Strength"].default_value = fibre.get("strength", 0.25)
        nt.links.new(tc.outputs["Object"], wv.inputs["Vector"]); nt.links.new(wv.outputs["Fac"], bp.inputs["Height"])
        nt.links.new(bp.outputs["Normal"], b.inputs["Normal"])
        b.inputs["Anisotropic"].default_value = fibre.get("aniso", 0.6)
    if glass:
        b.inputs["Transmission Weight"].default_value = 1.0; b.inputs["IOR"].default_value = glass.get("ior", 1.45)
        m.blend_method = "BLEND" if hasattr(m, "blend_method") else None
    if sss:
        b.inputs["Subsurface Weight"].default_value = sss; b.inputs["Subsurface Radius"].default_value = (0.4, 0.4, 0.4)
    if attr and flat_attr:                                  # colour-mapped data: emission only, so the colours are not shaded
        a = nt.nodes.new("ShaderNodeAttribute"); a.attribute_name = attr
        b.inputs["Base Color"].default_value = (0, 0, 0, 1); b.inputs["Specular IOR Level"].default_value = 0.0
        nt.links.new(a.outputs["Color"], b.inputs["Emission Color"]); b.inputs["Emission Strength"].default_value = emit or 1.0
    elif attr:
        a = nt.nodes.new("ShaderNodeAttribute"); a.attribute_name = attr
        nt.links.new(a.outputs["Color"], b.inputs["Base Color"])
        if emit:
            nt.links.new(a.outputs["Color"], b.inputs["Emission Color"]); b.inputs["Emission Strength"].default_value = emit
    else:
        b.inputs["Base Color"].default_value = (*srgb2lin(color), 1)
    if fade:                                                # display only: opaque below fade[0], fully transparent above fade[1] (world z)
        geo = nt.nodes.new("ShaderNodeNewGeometry"); sep = nt.nodes.new("ShaderNodeSeparateXYZ"); mr = nt.nodes.new("ShaderNodeMapRange")
        mr.inputs["From Min"].default_value = fade[0]; mr.inputs["From Max"].default_value = fade[1]; mr.clamp = True
        mr.interpolation_type = "SMOOTHSTEP"
        tr = nt.nodes.new("ShaderNodeBsdfTransparent"); mix = nt.nodes.new("ShaderNodeMixShader")
        nt.links.new(geo.outputs["Position"], sep.inputs[0]); nt.links.new(sep.outputs["Z"], mr.inputs["Value"])
        nt.links.new(mr.outputs["Result"], mix.inputs["Fac"]); nt.links.new(b.outputs["BSDF"], mix.inputs[1]); nt.links.new(tr.outputs["BSDF"], mix.inputs[2])
        nt.links.new(mix.outputs["Shader"], nt.nodes["Material Output"].inputs["Surface"])
    return m


def link(o):
    sc.collection.objects.link(o); return o


for i, m in enumerate(spec.get("meshes", [])):
    before = set(bpy.data.objects)
    if "ply" in m:
        bpy.ops.wm.ply_import(filepath=m["ply"])
    else:
        bpy.ops.wm.stl_import(filepath=m["stl"])
    for o in set(bpy.data.objects) - before:
        o.scale = [m.get("scale", 1.0)] * 3
        o.location = m.get("offset", [0, 0, 0])
        if m.get("smooth", True):
            for p in o.data.polygons:
                p.use_smooth = True
        o.data.materials.clear()
        attr = None
        if m.get("vertex_colors"):
            attr = o.data.color_attributes[0].name if len(o.data.color_attributes) else None
        o.data.materials.append(material(f"m{i}", m.get("color", (0.8, 0.8, 0.8)), m.get("rough", 0.45), m.get("metal", 0.0),
                                         coat=m.get("coat", 0.0), sss=m.get("sss", 0.0), fibre=m.get("fibre"), glass=m.get("glass"), spec=m.get("specular"),
                                         fade=m.get("fade"), attr=attr, emit=m.get("emit", 0.0), flat_attr=bool(attr)))
        if m.get("flat"):
            for p in o.data.polygons:
                p.use_smooth = False
        if m.get("no_shadow"):                              # lifted, display-only parts of exploded views
            o.visible_shadow = False

if "cylinders" in spec:
    cy = spec["cylinders"]; mat = material("fibre", cy.get("color", (0.82, 0.82, 0.82)), cy.get("rough", 0.35), coat=cy.get("coat", 0.0))
    L = cy["x1"] - cy["x0"]
    for k, (y, z) in enumerate(cy["centres"]):
        bpy.ops.mesh.primitive_cylinder_add(vertices=96, radius=cy["radius"], depth=L, location=((cy["x0"] + cy["x1"]) / 2, y, z),
                                            rotation=(0, math.pi / 2, 0))
        o = bpy.context.active_object; o.data.materials.append(mat)
        for p in o.data.polygons:
            p.use_smooth = len(p.vertices) == 4                           # flat end caps, smooth sides
        bev = o.modifiers.new("b", "BEVEL"); bev.width = cy["radius"] * 0.06; bev.segments = 3; bev.limit_method = "ANGLE"

_box_mats = {}
for i, bx in enumerate(spec.get("boxes", [])):
    bpy.ops.mesh.primitive_cube_add(size=1, location=bx["loc"], rotation=(0, 0, math.radians(bx.get("rot_z", 0))))
    o = bpy.context.active_object; o.scale = bx["size"]
    key = (tuple(bx.get("color", (0.6, 0.6, 0.6))), bx.get("rough", 0.5), bx.get("metal", 0.0))
    if key not in _box_mats:
        _box_mats[key] = material(f"box{len(_box_mats)}", key[0], key[1], key[2])
    o.data.materials.append(_box_mats[key])
    if bx.get("bevel"):
        b = o.modifiers.new("b", "BEVEL"); b.width = bx["bevel"]; b.segments = 2

_clip = []
def _clip_box():
    if not _clip:
        x0, x1, y0, y1, z0, z1 = spec["clip"]
        bpy.ops.mesh.primitive_cube_add(size=1, location=((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2))
        cb = bpy.context.active_object; cb.scale = (x1 - x0, y1 - y0, z1 - z0); cb.hide_render = True; cb.hide_viewport = True
        _clip.append(cb)
    return _clip[0]


_cyl_mats = {}
for i, c in enumerate(spec.get("cyls", [])):
    key = (tuple(c.get("color", (0.6, 0.6, 0.6))), c.get("rough", 0.5), c.get("metal", 0.0), c.get("coat", 0.0))
    if key not in _cyl_mats:
        _cyl_mats[key] = material(f"cyl{len(_cyl_mats)}", key[0], key[1], key[2], coat=key[3])
    rot = {"z": (0, 0, 0), "x": (0, math.pi / 2, 0), "y": (math.pi / 2, 0, 0)}[c.get("axis", "z")]
    bpy.ops.mesh.primitive_cylinder_add(vertices=c.get("vertices", 48), radius=c["radius"], depth=c["depth"], location=c["loc"], rotation=rot)
    o = bpy.context.active_object; o.data.materials.append(_cyl_mats[key])
    if "basis" in c:                                        # columns u, v, n: local x -> u, local z -> n
        from mathutils import Matrix
        u_, v_, n_ = (Vector(q) for q in c["basis"])
        o.rotation_mode = "QUATERNION"; o.rotation_quaternion = Matrix((u_, v_, n_)).transposed().to_quaternion()
    if "clip" in spec:
        bm = o.modifiers.new("clip", "BOOLEAN"); bm.operation = "INTERSECT"; bm.object = _clip_box(); bm.solver = "EXACT"
    for p in o.data.polygons:
        p.use_smooth = len(p.vertices) == 4

for i, cv in enumerate(spec.get("curves", [])):
    d = np.load(cv["npz"]); pts, sizes, col = d["pts"], d["sizes"], d["col"]
    hc = bpy.data.hair_curves.new(f"lines{i}"); hc.add_curves([int(s) for s in sizes])
    hc.attributes["position"].data.foreach_set("vector", pts.astype(np.float32).ravel())
    rad = hc.attributes.get("radius") or hc.attributes.new("radius", "FLOAT", "POINT")
    r = np.full(len(pts), cv["radius"], np.float32)
    if "taper" in d:
        r *= d["taper"].astype(np.float32)
    rad.data.foreach_set("value", r)
    ca = hc.attributes.new("col", "FLOAT_COLOR", "POINT")
    rgba = np.c_[srgb2lin(col), np.ones(len(col))].astype(np.float32); ca.data.foreach_set("color", rgba.ravel())
    o = link(bpy.data.objects.new(f"lines{i}", hc))
    hc.materials.append(material(f"lines{i}", attr="col", rough=cv.get("rough", 0.35), emit=cv.get("emit", 0.0), coat=cv.get("coat", 0.0), spec=cv.get("specular")))
    print(f"curves {len(sizes)} lines, {len(pts)} points")

if "ground" in spec:
    g = spec["ground"]
    bpy.ops.mesh.primitive_plane_add(size=g.get("size", 100), location=(g.get("x", 0), g.get("y", 0), g["z"]))
    o = bpy.context.active_object; o.is_shadow_catcher = True
    o.data.materials.append(material("ground", g.get("color", (1, 1, 1)), 0.9))

cam = spec["camera"]
cd = bpy.data.cameras.new("cam"); cd.lens = cam.get("lens", 50); cd.clip_end = 1e5; cd.clip_start = cam.get("clip_start", 0.01)
if cam.get("fstop"):
    cd.dof.use_dof = True; cd.dof.aperture_fstop = cam["fstop"]; cd.dof.focus_distance = cam["focus"]
co = link(bpy.data.objects.new("cam", cd)); co.location = cam["loc"]
co.rotation_euler = (Vector(cam["target"]) - Vector(cam["loc"])).to_track_quat("-Z", "Y").to_euler()
if cam.get("shift"):
    cd.shift_x, cd.shift_y = cam["shift"]
sc.camera = co

T = Vector(cam["target"]); C = Vector(cam["loc"]); dist = (C - T).length
fwd = (T - C).normalized(); right = fwd.cross(Vector((0, 0, 1))).normalized(); up = right.cross(fwd)
lt = spec.get("lights", {})
for name, (d_right, d_up, d_back, size, power) in {
        "key": (-0.9, 1.1, 0.6, 0.6, lt.get("key", 1.0)), "fill": (1.1, 0.2, 0.7, 1.0, lt.get("fill", 0.35)),
        "rim": (0.3, 0.9, -1.2, 0.5, lt.get("rim", 0.6))}.items():
    ld = bpy.data.lights.new(name, "AREA"); ld.size = size * dist; ld.energy = power * dist ** 2 * 120
    lo = link(bpy.data.objects.new(name, ld))
    lo.location = T + (right * d_right + up * d_up - fwd * d_back) * dist
    lo.rotation_euler = (T - lo.location).to_track_quat("-Z", "Y").to_euler()

# ---- annotation objects (annotated figures only): scale bars and flow arrows as real 3D objects,
# so their length is correct in perspective; their screen positions go to <out>.anchors.json
anchors = {}
_ann = material("annotation", (0.12, 0.12, 0.13), 0.6)
for i, item in enumerate(spec.get("bars", []) + spec.get("arrows", [])):
    a_, b_ = Vector(item["start"]), Vector(item["end"]); v = b_ - a_; L_ = v.length; r_ = item.get("radius", L_ * 0.01)
    head = item.get("head", 0.0)                          # arrows: cone head length (0 = plain bar)
    body = L_ - head
    bpy.ops.mesh.primitive_cylinder_add(vertices=24, radius=r_, depth=body, location=a_ + v.normalized() * body / 2)
    o = bpy.context.active_object; o.rotation_mode = "QUATERNION"; o.rotation_quaternion = v.to_track_quat("Z", "Y"); o.data.materials.append(_ann)
    if head:
        bpy.ops.mesh.primitive_cone_add(vertices=24, radius1=r_ * 3, depth=head, location=a_ + v.normalized() * (body + head / 2))
        c_ = bpy.context.active_object; c_.rotation_mode = "QUATERNION"; c_.rotation_quaternion = v.to_track_quat("Z", "Y"); c_.data.materials.append(_ann)
    from bpy_extras.object_utils import world_to_camera_view
    pa, pb = world_to_camera_view(sc, co, a_), world_to_camera_view(sc, co, b_)
    anchors[item.get("name", f"item{i}")] = dict(start=[pa.x, 1 - pa.y], end=[pb.x, 1 - pb.y], label=item.get("label", ""))

bpy.ops.render.render(write_still=True)
if anchors:
    json.dump(anchors, open(spec["out"] + ".anchors.json", "w"), indent=1)
print("saved", spec["out"])
