#include "physic.h"

/*
TODO:本文件中的内容有重复,将在新的一次重构中将其修正
*/

void cc::cell_class::copyconver() {
    conserformer = conser;
    if (!cc::urans) { tur->miubl_former = tur->miubl; }
}

void cc::cell_class::copyconver_time() {
    lastconser = conser;
    if (cc::urans) { tur->miubl_former = tur->miubl; }
}