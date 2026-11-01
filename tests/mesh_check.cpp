#include "readmesh.hpp"
#include <cstdio>

namespace {

int failures = 0;

void check(bool cond, const char* what) {
    if (!cond) {
        std::printf("FAIL: %s\n", what);
        ++failures;
    }
}

}

int main() {
    const bool ok = readmesh("cases/cyc/mesh.txt");
    check(ok, "readmesh returns true");

    check(cc::NodeList.size() == static_cast<std::size_t>(cc::node_num), "node count matches header");
    check(cc::FaceList.size() == static_cast<std::size_t>(cc::face_num), "face count matches header");
    check(cc::CellList.size() == static_cast<std::size_t>(cc::cell_num), "cell count matches header");
    check(cc::bounds.size() == static_cast<std::size_t>(cc::type_num), "bound count matches header");

    int bad_node = 0;
    int bad_nei = 0;
    int boundary_faces = 0;
    for (const auto& f : cc::FaceList) {
        if (f.node[0] == nullptr || f.node[1] == nullptr ||
            f.node[2] == nullptr || f.node[3] == nullptr) {
            ++bad_node;
        }
        if (f.type == cc::INTER) {
            if (f.nei[0] == nullptr || f.nei[1] == nullptr) ++bad_nei;
        } else {
            ++boundary_faces;
            if (f.nei[0] == nullptr || f.nei[1] != nullptr) ++bad_nei;
        }
    }
    check(bad_node == 0, "every face has 4 valid nodes");
    check(bad_nei == 0, "face to cell links are consistent");

    std::size_t linked = 0;
    for (const auto& b : cc::bounds) linked += b.fs.size();
    check(linked == static_cast<std::size_t>(boundary_faces), "bounds own every boundary face");

    int bad_cell = 0;
    for (const auto& c : cc::CellList) {
        for (int i = 0; i < 6; i++) {
            if (c.faces[i] == nullptr) ++bad_cell;
        }
    }
    check(bad_cell == 0, "every cell has 6 faces");

    std::printf("nodes=%zu faces=%zu cells=%zu bounds=%zu boundary_faces=%d\n",
                cc::NodeList.size(), cc::FaceList.size(), cc::CellList.size(),
                cc::bounds.size(), boundary_faces);
    std::printf("mesh_check: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
