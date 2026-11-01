import argparse
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description='Generate the Ma=0.2 cylinder mesh, 3D with symmetry span')
parser.add_argument('output', nargs='?', type=Path, default=Path(__file__).resolve().parent)
out = parser.parse_args().output.resolve()
out.mkdir(parents=True, exist_ok=True)
nt, nr, nz = 192, 128, 10
radius, outer_radius, first_height = 0.5, 30.0, 0.0015
span_length = 0.5
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


# 节点编号: i(环向) 周期共享, j(径向) 与 k(展向) 均为独立平面
def node(j, i, k):
  return (k * (nr + 1) + j) * nt + (i % nt) + 1


def coord(index):
  n = index - 1
  k, rest = divmod(n, (nr + 1) * nt)
  j, i = divmod(rest, nt)
  return radii[j] * cos_t[i], radii[j] * sin_t[i], zs[k]


def cell(i, j, k):
  return (i % nt) + nt * (j + nr * k) + 1


# 面按块连续编号, 顺序与写出顺序一致: 内部 i 向 -> 内部 j 向 -> 内部 k 向 -> 壁面 -> 远场 -> 展向对称
n_i = nt * nr * nz  # i 向面 (θ 向法向), 全部内部
n_j = nt * (nr - 1) * nz  # j 向面 (r 向法向), 仅内部
n_k = nt * nr * (nz - 1)  # k 向面 (z 向法向), 仅内部
n_cyl = nt * nz  # 圆柱壁面
n_far = nt * nz  # 远场
n_sym = 2 * nt * nr  # 展向两侧对称面
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


# 节点环按右手定则排列, 面法向由 0.5*(p2-p0)x(p3-p1) 给出, 由 cell_1 指向 cell_2
def loop_i(i, j, k):
  return (node(j, i, k), node(j, i, k + 1), node(j + 1, i, k + 1), node(j + 1, i, k))


def loop_j(i, j, k):
  return (node(j, i, k), node(j, i + 1, k), node(j, i + 1, k + 1), node(j, i, k + 1))


def loop_k(i, j, k):
  return (node(j, i, k), node(j + 1, i, k), node(j + 1, i + 1, k), node(j, i + 1, k))


def reverse(loop):
  return (loop[0],) + loop[:0:-1]


def normal(loop):
  p = [coord(n) for n in loop]
  d1 = [p[2][t] - p[0][t] for t in range(3)]
  d2 = [p[3][t] - p[1][t] for t in range(3)]
  return [
    0.5 * (d1[1] * d2[2] - d1[2] * d2[1]),
    0.5 * (d1[2] * d2[0] - d1[0] * d2[2]),
    0.5 * (d1[0] * d2[1] - d1[1] * d2[0]),
  ]


def dot(a, b):
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# 自检: 内部面法向由 cell_1 指向 cell_2, 边界面由流体指向外部
for i, j, k in ((0, 0, 0), (nt // 3, nr // 2, nz // 2), (nt - 1, nr - 1, nz - 1)):
  radial = (cos_t[i], sin_t[i], 0.0)
  tangent = (-sin_t[i], cos_t[i], 0.0)
  assert dot(normal(loop_i(i, j, k)), tangent) > 0, 'i 向面法向应朝 +θ'
  assert dot(normal(loop_j(i, j, k)), radial) > 0, 'j 向面法向应朝 +r'
  assert dot(normal(loop_k(i, j, k)), (0.0, 0.0, 1.0)) > 0, 'k 向面法向应朝 +z'
  assert dot(normal(reverse(loop_j(i, 0, k))), radial) < 0, '圆柱面法向应朝 -r'
  assert dot(normal(reverse(loop_k(i, j, 0))), (0.0, 0.0, 1.0)) < 0, '低侧对称面法向应朝 -z'
  assert dot(normal(loop_k(i, j, nz)), (0.0, 0.0, 1.0)) > 0, '高侧对称面法向应朝 +z'


def row(face, loop, left, right):
  return f'{face} {loop[0]} {loop[1]} {loop[2]} {loop[3]} {left} {right}\n'


with (out / 'mesh.txt').open('w') as f:
  f.write(
    f'{nt * (nr + 1) * (nz + 1)} {nfaces} {nt * nr * nz} 4\n'
    'interior = INTER\ncylinder = VIL\nfarfield = FAR\nspan = SYM\n(node)\n'
  )
  f.writelines(
    f'{node(j, i, k)} {radii[j] * cos_t[i]:.16e} {radii[j] * sin_t[i]:.16e} {zs[k]:.16e}\n'
    for k in range(nz + 1)
    for j in range(nr + 1)
    for i in range(nt)
  )
  f.write('(end)\n(edge)\ninterior\n')
  f.writelines(
    row(face_i(i, j, k), loop_i(i, j, k), cell(i - 1, j, k), cell(i, j, k))
    for k in range(nz)
    for j in range(nr)
    for i in range(nt)
  )
  f.writelines(
    row(face_j(i, j, k), loop_j(i, j, k), cell(i, j - 1, k), cell(i, j, k))
    for k in range(nz)
    for j in range(1, nr)
    for i in range(nt)
  )
  f.writelines(
    row(face_k(i, j, k), loop_k(i, j, k), cell(i, j, k - 1), cell(i, j, k))
    for k in range(1, nz)
    for j in range(nr)
    for i in range(nt)
  )
  f.write('(end)\ncylinder\n')
  f.writelines(
    row(face_cyl(i, k), reverse(loop_j(i, 0, k)), cell(i, 0, k), 0)
    for k in range(nz)
    for i in range(nt)
  )
  f.write('(end)\nfarfield\n')
  f.writelines(
    row(face_far(i, k), loop_j(i, nr, k), cell(i, nr - 1, k), 0)
    for k in range(nz)
    for i in range(nt)
  )
  f.write('(end)\nspan\n')
  f.writelines(
    row(face_sym(i, j, 0), reverse(loop_k(i, j, 0)), cell(i, j, 0), 0)
    for j in range(nr)
    for i in range(nt)
  )
  f.writelines(
    row(face_sym(i, j, 1), loop_k(i, j, nz), cell(i, j, nz - 1), 0)
    for j in range(nr)
    for i in range(nt)
  )
  f.write('(end)\n(end)\n(cell)\n')
  f.writelines(
    f'{cell(i, j, k)} '
    f'{face_i(i, j, k)} {face_i(i + 1, j, k)} '
    f'{face_cyl(i, k) if j == 0 else face_j(i, j, k)} '
    f'{face_far(i, k) if j == nr - 1 else face_j(i, j + 1, k)} '
    f'{face_sym(i, j, 0) if k == 0 else face_k(i, j, k)} '
    f'{face_sym(i, j, 1) if k == nz - 1 else face_k(i, j, k + 1)}\n'
    for k in range(nz)
    for j in range(nr)
    for i in range(nt)
  )
  f.write('(end)\n(all end)\n')
# 结构化邻接表: 首行 nt nr nz, 随后每行 i,j,k,cell
with (out / 'mesh_struct.txt').open('w') as f:
  f.write(f'{nt} {nr} {nz}\n')
  f.writelines(
    f'{i + 1},{j + 1},{k + 1},{cell(i, j, k)}\n'
    for k in range(nz)
    for j in range(nr)
    for i in range(nt)
  )
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
      span_boundary='symmetry',
    ),
    indent=2,
  )
  + '\n'
)
print((out / 'parameters.json').read_text())
