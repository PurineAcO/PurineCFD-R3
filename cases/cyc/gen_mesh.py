import argparse
import array
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description='Generate the Ma=0.2 cylinder mesh, 3D with symmetry span')
parser.add_argument('output', nargs='?', type=Path, default=Path(__file__).resolve().parent)
out = parser.parse_args().output.resolve()
out.mkdir(parents=True, exist_ok=True)
nt, nr, nz = 100, 80, 30
radius, outer_radius, first_height = 0.5, 30.0, 0.02
span_length = 3.0
span_bc = 'SYM'
lo, hi = 1.0, 1.2
for _ in range(80):
  q = (lo + hi) / 2
  if first_height * (q**nr - 1) / (q - 1) > outer_radius - radius:
    hi = q
  else:
    lo = q
q = (lo + hi) / 2
radii = [radius + first_height * (q**j - 1) / (q - 1) for j in range(nr + 1)]
cos_t = [math.cos(2 * math.pi * i / nt) for i in range(nt)]
sin_t = [math.sin(2 * math.pi * i / nt) for i in range(nt)]
zs = [span_length * k / nz for k in range(nz + 1)]
nnode = nt * (nr + 1) * (nz + 1)
ncells = nt * nr * nz


def node(j, i, k):
  return (k * (nr + 1) + j) * nt + (i % nt) + 1


def cell(i, j, k):
  return (i % nt) + nt * (j + nr * k) + 1


def table(typecode):
  return array.array(typecode, [0]) * (nnode + 1)


NI, NJ, NK = table('h'), table('h'), table('h')
XS, YS, ZS = table('d'), table('d'), table('d')
row_i = array.array('h', range(nt))
for k in range(nz + 1):
  for j in range(nr + 1):
    base = (k * (nr + 1) + j) * nt + 1
    NI[base:base + nt] = row_i
    NJ[base:base + nt] = array.array('h', [j]) * nt
    NK[base:base + nt] = array.array('h', [k]) * nt
    XS[base:base + nt] = array.array('d', [radii[j] * c for c in cos_t])
    YS[base:base + nt] = array.array('d', [radii[j] * s for s in sin_t])
    ZS[base:base + nt] = array.array('d', [zs[k]]) * nt

n_i = nt * nr * nz
n_j = nt * (nr - 1) * nz
n_k = nt * nr * (nz - 1)
n_cyl = nt * nz
n_far = nt * nz
n_sym = 2 * nt * nr
nfaces = n_i + n_j + n_k + n_cyl + n_far + n_sym


def face_i(i, j, k):
  return (i % nt) + nt * (j + nr * k) + 1


def face_j(i, j, k):
  return n_i + (i % nt) + nt * (j - 1 + (nr - 1) * k) + 1


def face_k(i, j, k):
  return n_i + n_j + (i % nt) + nt * (j + nr * (k - 1)) + 1


def face_cyl(i, k):
  return n_i + n_j + n_k + (i % nt) + nt * k + 1


def face_far(i, k):
  return n_i + n_j + n_k + n_cyl + (i % nt) + nt * k + 1


def face_sym(i, j, side):
  return n_i + n_j + n_k + n_cyl + n_far + side * nt * nr + (i % nt) + nt * j + 1


def loop_i(i, j, k):
  return (node(j, i, k), node(j, i, k + 1), node(j + 1, i, k + 1), node(j + 1, i, k))


def loop_j(i, j, k):
  return (node(j, i, k), node(j, i + 1, k), node(j, i + 1, k + 1), node(j, i, k + 1))


def loop_k(i, j, k):
  return (node(j, i, k), node(j + 1, i, k), node(j + 1, i + 1, k), node(j, i + 1, k))


def reverse(loop):
  return (loop[0],) + loop[:0:-1]


def cell_sides(i, j, k):
  return (
    (face_i(i, j, k), loop_i(i, j, k)),
    (face_i(i + 1, j, k), loop_i(i + 1, j, k)),
    (face_cyl(i, k), reverse(loop_j(i, 0, k)))
    if j == 0 else (face_j(i, j, k), loop_j(i, j, k)),
    (face_far(i, k), loop_j(i, nr, k))
    if j == nr - 1 else (face_j(i, j + 1, k), loop_j(i, j + 1, k)),
    (face_sym(i, j, 0), reverse(loop_k(i, j, 0)))
    if k == 0 else (face_k(i, j, k), loop_k(i, j, k)),
    (face_sym(i, j, 1), loop_k(i, j, nz))
    if k == nz - 1 else (face_k(i, j, k + 1), loop_k(i, j, k + 1)),
  )


def normal(loop):
  a, b, c, d = loop
  ax, ay, az = XS[c] - XS[a], YS[c] - YS[a], ZS[c] - ZS[a]
  bx, by, bz = XS[d] - XS[b], YS[d] - YS[b], ZS[d] - ZS[b]
  return (
    0.5 * (ay * bz - az * by),
    0.5 * (az * bx - ax * bz),
    0.5 * (ax * by - ay * bx),
  )


def cell_center(c):
  n = c - 1
  i = n % nt
  rest = n // nt
  j = rest % nr
  k = rest // nr
  rm = 0.5 * (radii[j] + radii[j + 1])
  return (
    rm * 0.5 * (cos_t[i] + cos_t[(i + 1) % nt]),
    rm * 0.5 * (sin_t[i] + sin_t[(i + 1) % nt]),
    0.5 * (zs[k] + zs[k + 1]),
  )


def midpoint(loop):
  a, b, c, d = loop
  return (
    (XS[a] + XS[b] + XS[c] + XS[d]) * 0.25,
    (YS[a] + YS[b] + YS[c] + YS[d]) * 0.25,
    (ZS[a] + ZS[b] + ZS[c] + ZS[d]) * 0.25,
  )


def face_error(loop):
  a, b, c, d = loop
  if a == b or a == c or a == d or b == c or b == d or c == d:
    return '顶点重复'
  for p, r in ((a, b), (b, c), (c, d), (d, a)):
    if (NI[p] != NI[r]) + (NJ[p] != NJ[r]) + (NK[p] != NK[r]) != 1:
      return '相邻顶点不是同一条网格边, 环序被破坏'
  for p, r in ((a, c), (b, d)):
    if (NI[p] != NI[r]) + (NJ[p] != NJ[r]) + (NK[p] != NK[r]) != 2:
      return '对角线上的顶点被写成相邻'
  return None


def cell_error(loops):
  sides = [set(loop) for loop in loops]
  shared = set().union(*sides)
  if len(shared) != 8:
    return f'不是八节点六面体(共 {len(shared)} 个顶点)'
  for p in (0, 2, 4):
    if not sides[p].isdisjoint(sides[p + 1]):
      return f'f{p + 1} 与 f{p + 2} 不是相对面'
  for p in range(6):
    for r in range(p + 1, 6):
      if r == p + 1 and p % 2 == 0:
        continue
      if len(sides[p] & sides[r]) != 2:
        return f'f{p + 1} 与 f{r + 1} 未共用一条棱'
  return None


def face_row(face, loop, left, right):
  bad = face_error(loop)
  if bad is not None:
    raise SystemExit(f'face {face} (cell {left} / {right}): {bad}')
  nor, mid = normal(loop), midpoint(loop)
  s = cell_center(left)
  if nor[0] * (mid[0] - s[0]) + nor[1] * (mid[1] - s[1]) + nor[2] * (mid[2] - s[2]) <= 0.0:
    raise SystemExit(f'face {face}: 法向未由 cell {left} 指向外侧')
  if right:
    s = cell_center(right)
    if nor[0] * (mid[0] - s[0]) + nor[1] * (mid[1] - s[1]) + nor[2] * (mid[2] - s[2]) >= 0.0:
      raise SystemExit(f'face {face}: 法向未由 cell {right} 指向内侧')
  return f'{face} {loop[0]} {loop[1]} {loop[2]} {loop[3]} {left} {right}\n'


def cell_row(i, j, k):
  sides = cell_sides(i, j, k)
  bad = cell_error([loop for _, loop in sides])
  if bad is not None:
    raise SystemExit(f'cell {cell(i, j, k)} (i={i}, j={j}, k={k}): {bad}')
  return ' '.join(str(v) for v in (cell(i, j, k), *(face for face, _ in sides))) + '\n'


with (out / 'mesh.txt').open('w') as f:
  f.write(
    f'{nnode} {nfaces} {ncells} 4\n'
    'interior = INTER\ncylinder = WALL\nfarfield = FAR\n'
    f'span = {span_bc}\n(node)\n'
  )
  f.writelines(
    f'{node(j, i, k)} {radii[j] * cos_t[i]:.16e} {radii[j] * sin_t[i]:.16e} {zs[k]:.16e}\n'
    for k in range(nz + 1)
    for j in range(nr + 1)
    for i in range(nt)
  )
  f.write('(end)\n(edge)\ninterior\n')
  f.writelines(
    face_row(face_i(i, j, k), loop_i(i, j, k), cell(i - 1, j, k), cell(i, j, k))
    for k in range(nz)
    for j in range(nr)
    for i in range(nt)
  )
  f.writelines(
    face_row(face_j(i, j, k), loop_j(i, j, k), cell(i, j - 1, k), cell(i, j, k))
    for k in range(nz)
    for j in range(1, nr)
    for i in range(nt)
  )
  f.writelines(
    face_row(face_k(i, j, k), loop_k(i, j, k), cell(i, j, k - 1), cell(i, j, k))
    for k in range(1, nz)
    for j in range(nr)
    for i in range(nt)
  )
  f.write('(end)\ncylinder\n')
  f.writelines(
    face_row(face_cyl(i, k), reverse(loop_j(i, 0, k)), cell(i, 0, k), 0)
    for k in range(nz)
    for i in range(nt)
  )
  f.write('(end)\nfarfield\n')
  f.writelines(
    face_row(face_far(i, k), loop_j(i, nr, k), cell(i, nr - 1, k), 0)
    for k in range(nz)
    for i in range(nt)
  )
  f.write('(end)\nspan\n')
  f.writelines(
    face_row(face_sym(i, j, 0), reverse(loop_k(i, j, 0)), cell(i, j, 0), 0)
    for j in range(nr)
    for i in range(nt)
  )
  f.writelines(
    face_row(face_sym(i, j, 1), loop_k(i, j, nz), cell(i, j, nz - 1), 0)
    for j in range(nr)
    for i in range(nt)
  )
  f.write('(end)\n(end)\n(cell)\n')
  f.writelines(
    cell_row(i, j, k) for k in range(nz) for j in range(nr) for i in range(nt)
  )
  f.write('(end)\n(all end)\n')
with (out / 'mesh_struct.txt').open('w') as f:
  f.write(f'{nt} {nr} {nz}\n')
  f.writelines(
    f'{i + 1},{j + 1},{k + 1},{cell(i, j, k)}\n'
    for k in range(nz)
    for j in range(nr)
    for i in range(nt)
  )
print(f'verified {nfaces} faces (环序/无对角相邻) and {ncells} cells (f1-f2/f3-f4/f5-f6 相对, 法向朝外)')
T, Ma, Re, diameter = 300.0, 0.2, 1000.0, 1.0
R = 287.05
mu = 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)
U = Ma * math.sqrt(1.4 * R * T)
rho = Re * mu / (U * diameter)
p = rho * R * T
(out / 'parameters.json').write_text(
  json.dumps(
    dict(
      Ma=Ma,
      Re=Re,
      T_K=T,
      D_m=diameter,
      U_m_s=U,
      rho_kg_m3=rho,
      p_Pa=p,
      mu_Pa_s=mu,
      nt=nt,
      nr=nr,
      nz=nz,
      outer_radius_D=outer_radius,
      first_height_D=first_height,
      growth_ratio=q,
      span_length_D=span_length,
      span_spacing_D=span_length / nz,
      span_boundary=span_bc.lower(),
    ),
    indent=2,
  )
  + '\n'
)
print((out / 'parameters.json').read_text())
