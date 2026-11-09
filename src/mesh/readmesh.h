#pragma once



/*
在设计初期,我们会在所有文件中硬编码六面体棱柱层网格,在以后会逐步修改
*/

// 读取非结构网格
bool readmesh(const char* meshpath);
// 建立面与单元的邻接关系
bool linkmesh();
// 返回一个标识
void typetoken(char* a,char* b);



