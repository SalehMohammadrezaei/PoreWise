"""Put a transparent Blender render on its background.

    python compose.py <rgba.png> <out.png> cv|dark

cv:   pure white background (the CV page style).
dark: near-black background, plus a soft glow of the bright, saturated lines (screen blend).
"""
import sys
import numpy as np
from PIL import Image, ImageFilter

src, out, style = sys.argv[1:4]
im = Image.open(src).convert("RGBA")
if style == "cv":
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255)); bg.alpha_composite(im); bg.convert("RGB").save(out, optimize=True)
else:
    bg = Image.new("RGBA", im.size, (8, 10, 14, 255)); bg.alpha_composite(im)
    base = np.asarray(bg.convert("RGB")).astype(np.float32) / 255
    a = np.asarray(im)[..., 3:4].astype(np.float32) / 255
    rgb = np.asarray(im.convert("RGB")).astype(np.float32) / 255
    sat = rgb.max(-1, keepdims=True) - rgb.min(-1, keepdims=True)            # glow from coloured lines, not grey solids
    src_glow = Image.fromarray((np.clip(rgb * a * np.clip(sat * 2.5, 0, 1), 0, 1) * 255).astype(np.uint8))
    glow = sum(np.asarray(src_glow.filter(ImageFilter.GaussianBlur(r))).astype(np.float32) / 255 * w for r, w in ((6, 0.22), (24, 0.16), (70, 0.10)))
    res = 1 - (1 - base) * (1 - np.clip(glow, 0, 1))
    Image.fromarray((np.clip(res, 0, 1) * 255).astype(np.uint8)).save(out, optimize=True)
print("saved", out, im.size)
