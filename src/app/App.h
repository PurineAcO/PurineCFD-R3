#pragma once

#include "data_structure/DataStructure.h"
#include <string>
#include <vector>

// 扫描 root 下已有的算例目录编号, 返回下一个可用目录名
std::string next_case_dir(const std::string& root, const std::string& tag);

// 在圆柱附近按幅值 amp 叠加初始扰动, 并刷新引申物理量与守恒量
void add_seed(double amp);

// 保存当前流场快照, 供定常迭代发散时回退
void save_state();

// 回退到最近一次快照
void restore_state();

extern std::vector<cc::vec5> snap_q;
extern std::vector<double> snap_nu;
extern double cfl_scale;
