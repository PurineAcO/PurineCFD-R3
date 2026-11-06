#include "geometry.hpp"
#include "implict.hpp"
#include "readmesh.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

int failures = 0;

void check(bool cond, const char* what) {
    if (cond) return;
    std::printf("FAIL: %s\n", what);
    ++failures;
}

double outward_dot(const cc::face_class& f, const cc::cell_class& c) {
    const cc::vec3 p0 = f.node[0]->p, p1 = f.node[1]->p;
    const cc::vec3 p2 = f.node[2]->p, p3 = f.node[3]->p;
    const cc::vec3 n = 0.5 * cc::cross(p2 - p0, p3 - p1);
    const cc::vec3 mid = (p0 + p1 + p2 + p3) * 0.25;
    return (f.nei[0] == &c ? 1.0 : -1.0) * cc::dot(mid, n);
}

double divergence_volume(const cc::cell_class& c) {
    double v = 0.0;
    for (int k = 0; k < c.fcnt; ++k) v += outward_dot(*c.faces[k], c);
    return v / 3.0;
}

}


static cc::vec5 phys_flux_n(const cc::vec3& nor, const cc::vec5& q) {
    const double rho = q.c, u = q.x / rho, v = q.y / rho, w = q.z / rho;
    const double ke = 0.5 * (u * u + v * v + w * w);
    const double p = (cc::gamma - 1.0) * (q.e - q.c * ke);
    const double un = u * nor.x + v * nor.y + w * nor.z;
    return cc::vec5(rho * un, rho * u * un + p * nor.x, rho * v * un + p * nor.y,
                    rho * w * un + p * nor.z, un * (q.e + p));
}

static void check_euler_jacobi() {
    cc::face_class f;
    f.nor = cc::vec3(1.3, -0.7, 0.45);
    f.phy = cc::vecp(1.7, 210.0, -55.0, 33.0, 412.0);
    f.otphy = form_otphy(f.phy);
    const double rho = f.phy.rho;
    const cc::vec5 q(rho, rho * f.phy.u, rho * f.phy.v, rho * f.phy.w,
                     rho * get_energy(f.phy));
    const cc::mat5 A = impl::euler_jacobi_n(f);
    double worst = 0.0;
    for (int k = 0; k < 5; k++) {
        cc::vec5 dir;
        (&dir.c)[k] = 1.0;
        const double h = 1e-6 * (k == 0 ? q.c : (k == 4 ? q.e : 1.0));
        const cc::vec5 an = A * dir;
        const cc::vec5 fd = (1.0 / (2.0 * h)) * (phys_flux_n(f.nor, q + h * dir) -
                                                 phys_flux_n(f.nor, q - h * dir));
        const double scale = std::max(1.0, std::max(std::abs(an.e), std::abs(fd.e)));
        worst = std::max(worst, std::abs(an.e - fd.e) / scale);
    }
    check(worst < 1e-6, "euler jacobi matches finite difference of physical flux");
    std::printf("jacobi: worst relative diff=%.3e\n", worst);
}

int main() {
    if (!readmesh("cases/cyc/mesh.txt")) {
        std::printf("check: FAIL\n");
        return 1;
    }

    const auto nodes = static_cast<std::size_t>(cc::node_num);
    const auto faces = static_cast<std::size_t>(cc::face_num);
    const auto cells = static_cast<std::size_t>(cc::cell_num);
    check(cc::NodeList.size() == nodes, "node count matches header");
    check(cc::FaceList.size() == faces, "face count matches header");
    check(cc::CellList.size() == cells, "cell count matches header");
    check(cc::bounds.size() == static_cast<std::size_t>(cc::type_num), "bound count matches header");
    int bad_order = 0;
    for (std::size_t i = 0; i < faces; ++i) {
        bad_order += cc::FaceList[i].index != static_cast<int>(i) + 1;
    }
    for (std::size_t i = 0; i < cells; ++i) {
        bad_order += cc::CellList[i].index != static_cast<int>(i) + 1;
    }
    check(bad_order == 0, "list position matches file index");

    std::vector<int> refs(faces + 1, 0);
    int bad_link = 0;
    int bad_nei = 0;
    allcell {
        cc::cell_class& c = icell(i);
        for (int k = 0; k < c.fcnt; ++k) {
            const cc::face_class* f = c.faces[k];
            if (f->nei[0] != &c && f->nei[1] != &c) {
                ++bad_link;
                continue;
            }
            ++refs[f->index];
            const cc::cell_class* other = c.nei[k];
            const bool interior = f->type == cc::INTER;
            if (interior ? (other == nullptr || other == &c) : other != nullptr) ++bad_nei;
        }
    }
    check(bad_link == 0, "every cell face lists that cell");
    check(bad_nei == 0, "cell neighbour agrees with face type");

    std::vector<char> used(nodes, 0);
    int bad_ref = 0;
    int bad_node = 0;
    int boundary = 0;
    for (const cc::face_class& f : cc::FaceList) {
        const int want = f.type == cc::INTER ? 2 : 1;
        if (refs[f.index] != want) ++bad_ref;
        boundary += want == 1;
        for (const cc::node_class* n : f.node) {
            if (n->index < 1 || n->index > cc::node_num || &cc::NodeList[n->index - 1] != n) {
                ++bad_node;
                continue;
            }
            used[n->index - 1] = 1;
        }
    }
    check(bad_ref == 0, "interior faces shared by 2 cells, boundary faces by 1");
    check(bad_node == 0, "face nodes point into node list");
    int unused = 0;
    for (char u : used) unused += u == 0;
    check(unused == 0, "every node is referenced by a face");

    std::size_t owned = 0;
    for (const cc::bound& b : cc::bounds) owned += b.fs.size();
    check(owned == static_cast<std::size_t>(boundary), "bounds own every boundary face");

    double worst = 0.0;
    double sum_tet = 0.0;
    double sum_div = 0.0;
    int bad_vol = 0;
    allcell {
        cc::cell_class& c = icell(i);
        const double v_tet = cube_volumn_cen(c).first;
        const double v_div = divergence_volume(c);
        sum_tet += v_tet;
        sum_div += v_div;
        const double rel = std::abs(v_tet - v_div) / std::abs(v_div);
        if (!(rel <= 1e-9)) ++bad_vol;
        if (rel > worst) worst = rel;
    }
    check(bad_vol == 0, "hex volume: tetrahedron split matches divergence theorem");

    double sum_edge = 0.0;
    for (const cc::bound& b : cc::bounds) {
        if (b.type == cc::INTER) continue;
        for (const cc::face_class* f : b.fs) {
            sum_edge += outward_dot(*f, *f->nei[0]);
        }
    }
    sum_edge /= 3.0;
    check(std::abs(sum_tet - sum_edge) / std::abs(sum_edge) < 1e-9,
          "total volume matches boundary-only integral");

    double wall_r2 = 0.0, far_r2 = 0.0, z_lo = 0.0, z_hi = 0.0;
    int wall_faces = 0;
    std::vector<double> zs;
    for (const cc::bound& b : cc::bounds) {
        for (const cc::face_class* f : b.fs) {
            for (const cc::node_class* n : f->node) {
                const double r2 = n->p.x * n->p.x + n->p.y * n->p.y;
                if (b.type == cc::WALL) wall_r2 = r2;
                else if (b.type == cc::FAR) far_r2 = r2;
                if (zs.empty()) { z_lo = z_hi = n->p.z; }
                if (n->p.z < z_lo) z_lo = n->p.z;
                if (n->p.z > z_hi) z_hi = n->p.z;
                bool seen = false;
                for (double z : zs) seen = seen || std::abs(z - n->p.z) < 1e-12;
                if (!seen) zs.push_back(n->p.z);
            }
            wall_faces += b.type == cc::WALL;
        }
    }
    const int nz = static_cast<int>(zs.size()) - 1;
    const int nt = wall_faces > 0 && nz > 0 ? wall_faces / nz : 0;
    const double pi = 3.14159265358979323846;
    const double closed = nt > 0 ? 0.5 * nt * std::sin(2.0 * pi / nt) * (far_r2 - wall_r2) * (z_hi - z_lo) : 0.0;
    check(nt > 0 && std::abs(sum_tet - closed) / std::abs(closed) < 1e-9,
          "total volume matches polygon prism formula");

    std::printf("volume: tet=%.9f div=%.9f edge=%.9f closed=%.9f (nt=%d nz=%d)\n", sum_tet, sum_div,
                sum_edge, closed, nt, nz);
    std::printf("volume: worst per-cell relative diff=%.3e\n", worst);

    std::printf("nodes=%zu faces=%zu cells=%zu bounds=%zu boundary_faces=%d\n", nodes, faces, cells,
                cc::bounds.size(), boundary);
    check_euler_jacobi();

    std::printf("check: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
