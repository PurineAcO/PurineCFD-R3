#pragma once

// 读取非结构网格
// 注：设计初期会在所有文件中硬编码六面体棱柱层网格，以后会逐步修改
bool readmesh(const char* meshpath);
// 读取二维非结构网格
bool readmesh_2D(const char* meshpath);
// 建立面与单元的邻接关系
bool linkmesh();
// 返回一个标识
void typetoken(char* a,char* b);
