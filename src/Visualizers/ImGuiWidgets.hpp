#pragma once
#include <imgui.h>
#include <vector>
#include <imgui_internal.h>


void DrawLegendScale(ImVec2 size, float min, float max);
int DrawPieChart(const std::vector<float>& values, const std::vector<ImU32>& colors, float radius);
int StackedProgressBar(const std::vector<float>& values, const ImVec2& size = ImVec2(-FLT_MIN, 0.0f));
bool EditProbe(ImVec2 pos, bool* open, float p1[2], float p2[2], char* buf, size_t size);
