#include "readmesh.hpp"
#include <cstdio>
#include <vector>

namespace {

int failures = 0;

void check(bool cond, const char* what) {
    if (cond) return;
    std::printf("FAIL: %s\n", what);
    ++failures;
}

} // namespace

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

    // 单元侧: 每个面确实属于该单元, 且邻接指针与面类型一致
    std::vector<int> refs(faces + 1, 0);
    int bad_link = 0;
    int bad_nei = 0;
    allcell {
        cc::cell_class& c = icell(i);
        for (int k = 0; k < c.ecnt; ++k) {
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

    // 面侧: 内部面被两个单元共享, 边界面被一个单元共享; 节点指针有效并统计使用
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

    std::printf("nodes=%zu faces=%zu cells=%zu bounds=%zu boundary_faces=%d\n", nodes, faces, cells,
                cc::bounds.size(), boundary);
    std::printf("check: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
