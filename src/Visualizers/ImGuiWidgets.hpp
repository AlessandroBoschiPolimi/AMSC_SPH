#pragma once
#include <imgui.h>
#include <vector>
#include <imgui_internal.h>

#include "Utility.hpp"



void DrawLegendScale(ImVec2 size, float min, float max);
int DrawPieChart(const std::vector<float>& values, const std::vector<ImU32>& colors, float radius);
int DrawStackedProgressBar(const std::vector<float>& values, const ImVec2& size = ImVec2(-FLT_MIN, 0.0f));

ImVec2 ToScreenSpace(const coord<float, 2>& p, ImVec2 min, ImVec2 range, ImVec2 canvasPos, ImVec2 canvasSize);
GraphRange DrawGraph(const std::vector<coord<float, 2>>& points, ImVec2 size, const std::function<void(ImDrawList* drawlist, ImVec2 min, ImVec2 max, ImVec2 range, ImVec2 canvasPos, ImVec2 canvasSize)>& custom_draw);

bool BeginEditProbeWindow(ImVec2 pos, bool* open, float p1[2], float p2[2], char* buf, size_t size, bool copy, bool paste);
