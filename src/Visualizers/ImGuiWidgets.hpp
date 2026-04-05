#pragma once
#include <imgui.h>
#include <vector>
#include <imgui_internal.h>

#include "Utility.hpp"


void DrawLegendScale(ImVec2 size, float min, float max);
int DrawPieChart(const std::vector<float>& values, const std::vector<ImU32>& colors, float radius);
int DrawStackedProgressBar(const std::vector<float>& values, const ImVec2& size = ImVec2(-FLT_MIN, 0.0f));

enum class DrawGraphMode
{
	NONE = 0,
	LINE = 1 << 0,
	POINT = 1 << 1
};
std::underlying_type_t<DrawGraphMode> operator&(DrawGraphMode a, DrawGraphMode b);
DrawGraphMode operator|(DrawGraphMode a, DrawGraphMode b);

struct AxisRange
{
	float min, max;
};
struct GraphRange
{
	AxisRange x, y;
};

GraphRange DrawGraph(const std::vector<coord<float, 2>>& points, ImVec2 size, DrawGraphMode mode, const consumer<ImDrawList*>& custom_draw);

bool BeginEditProbeWindow(ImVec2 pos, bool* open, float p1[2], float p2[2], char* buf, size_t size, bool copy, bool paste);
