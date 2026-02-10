#include "ImGuiWidgets.hpp"

float Lerp(float a, float b, float t)
{
	return a + (b - a) * t;
}
void DrawLegendScale(ImVec2 size, float min, float max)
{
	ImGui::Text("%.2f", min);
	ImGui::SameLine();

	ImDrawList* draw_list = ImGui::GetWindowDrawList();
	ImVec2 p0 = ImGui::GetCursorScreenPos();
	ImVec2 p1 = ImVec2(p0.x + size.x, p0.y + size.y);

	const int steps = 100; // more = smoother

	for (int i = 0; i < steps; ++i)
	{
		float t0 = (float)i / steps;
		float t1 = (float)(i + 1) / steps;

		// Long-way hue interpolation: 2/3 -> 0
		float h0 = (2.0f / 3.0f) * (1.0f - t0);
		float h1 = (2.0f / 3.0f) * (1.0f - t1);

		ImU32 c0 = ImColor::HSV(h0, 1.0f, 1.0f);
		ImU32 c1 = ImColor::HSV(h1, 1.0f, 1.0f);

		float x0 = Lerp(p0.x, p1.x, t0);
		float x1 = Lerp(p0.x, p1.x, t1);

		draw_list->AddRectFilledMultiColor(
			ImVec2(x0, p0.y),
			ImVec2(x1, p1.y),
			c0, c1, c1, c0
		);
	}

	// Advance layout cursor
	ImGui::Dummy(size);
	ImGui::SameLine();
	ImGui::Text("%.2f", max);
}


float AngleFromCenter(ImVec2 center, ImVec2 point)
{
	return atan2f(point.y - center.y, point.x - center.x);
}
int DrawPieChart(const std::vector<float>& values, const std::vector<ImU32>& colors, float radius)
{
	ImDrawList* draw_list = ImGui::GetWindowDrawList();

	ImVec2 pos = ImGui::GetCursorScreenPos();
	ImVec2 center = { pos.x + radius, pos.y + radius };

	ImVec2 mouse = ImGui::GetIO().MousePos;
	float mouse_dist = (mouse.x - center.x) * (mouse.x - center.x) + (mouse.y - center.y) * (mouse.y - center.y);
	float mouse_angle = AngleFromCenter(center, mouse);

	float total = 0.0f;
	for (float v : values)
		total += v;

	float start_angle = -IM_PI / 2.0f;
	int hovered_index = -1;

	for (int i = 0; i < values.size(); i++)
	{
		float slice_angle = (values[i] / total) * IM_PI * 2.0f;
		float end_angle = start_angle + slice_angle;

		bool hovered = false;

		if (mouse_dist <= radius * radius)
		{
			// Normalize angle range
			float a = mouse_angle;
			if (a < start_angle) a += IM_PI * 2.0f;

			hovered = (a >= start_angle && a <= end_angle);
		}

		ImU32 color = colors[i];

		// Brighten hovered slice
		if (hovered)
		{
			hovered_index = i;
			ImVec4 c = ImGui::ColorConvertU32ToFloat4(color);
			c.x = ImMin(c.x * 1.2f, 1.0f);
			c.y = ImMin(c.y * 1.2f, 1.0f);
			c.z = ImMin(c.z * 1.2f, 1.0f);
			color = ImGui::ColorConvertFloat4ToU32(c);
		}

		draw_list->PathClear();
		draw_list->PathArcTo(center, radius, start_angle, end_angle, 32);
		draw_list->PathLineTo(center);
		draw_list->PathFillConvex(color);

		start_angle = end_angle;
	}

	ImGui::Dummy(ImVec2(radius * 2, radius * 2));
	return hovered_index;
}



int StackedProgressBar(const std::vector<float>& values, const ImVec2& size)
{
	ImGuiWindow* window = ImGui::GetCurrentWindow();
	if (window->SkipItems)
		return -1;

	float sum = 0.0f;
	for (float v : values)
		sum += ImMax(v, 0.0f);

	if (sum <= 0.0f)
		return -1;

	ImVec2 pos = ImGui::GetCursorScreenPos();
	ImVec2 bar_size = ImGui::CalcItemSize(size, ImGui::GetContentRegionAvail().x, ImGui::GetFrameHeight());

	ImRect bb(pos, ImVec2(pos.x + bar_size.x, pos.y + bar_size.y));
	ImGui::ItemSize(bb);
	if (!ImGui::ItemAdd(bb, 0))
		return -1;

	// Background
	ImGui::RenderFrame(bb.Min, bb.Max, ImGui::GetColorU32(ImGuiCol_FrameBg), true, ImGui::GetStyle().FrameRounding);

	ImDrawList* draw = ImGui::GetWindowDrawList();

	std::vector<ImU32> colors =
	{
		IM_COL32(255, 80, 80, 255),
		IM_COL32(80, 200, 80, 255),
		IM_COL32(80, 120, 255, 255)
	};

	int hovering = -1;

	float x = bb.Min.x;
	for (int i = 0; i < values.size(); i++)
	{
		float v = ImMax(values[i], 0.0f);
		if (v == 0.0f)
			continue;

		float w = (v / sum) * bar_size.x;
		if (w <= 0.0f)
			continue;

		ImU32 col = colors[i % colors.size()];

		draw->AddRectFilled(
			ImVec2(x, bb.Min.y), ImVec2(x + w, bb.Max.y),
			col, ImGui::GetStyle().FrameRounding,
			(x == bb.Min.x) ? ImDrawFlags_RoundCornersLeft : 0);

		if (ImGui::IsMouseHoveringRect(ImVec2(x, bb.Min.y), ImVec2(x + w, bb.Max.y)))
		{
			hovering = i;
		}

		x += w;
	}

	return hovering;
}