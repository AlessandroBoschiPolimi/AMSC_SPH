#include "ImGuiWidgets.hpp"
#include <string>
#include "Probe.hpp"


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



int DrawStackedProgressBar(const std::vector<float>& values, const ImVec2& size)
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


bool BeginEditProbeWindow(ImVec2 pos, bool* open, float p1[2], float p2[2], char* buf, size_t size, bool copy, bool paste)
{
	ImGui::SetNextWindowPos(pos, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
	ImGui::SetNextWindowSize(ImVec2(400, 150), ImGuiCond_FirstUseEver);

	static stdclock::time_point error_time;
	static std::string error_msg = "";
	static bool error_show = false;

	ImGui::Begin("Add Probe", open);

	ImGui::TextUnformatted("Name: "); ImGui::SameLine();
	ImGui::InputText("##name", buf, size);

	float lw = std::max(ImGui::CalcTextSize("Top Left: ").x, ImGui::CalcTextSize("Bottom Right: ").x) * 1.1;

	ImGui::TextUnformatted("Top Left: "); ImGui::SameLine();
	ImGui::SetCursorPosX(lw);
	ImGui::InputFloat2("##apc1", p1); // add probe coord

	ImGui::TextUnformatted("Bottom Right: "); ImGui::SameLine();
	ImGui::SetCursorPosX(lw);
	ImGui::InputFloat2("##apc2", p2);

	if (ImGui::Button("Confirm"))
	{
		ImGui::End();
		return true;
	}
	ImGui::SameLine();

	if (ImGui::Button("Copy") || copy)
	{
		Probe p;
		p.Name = buf;
		p.TL = { p1[0], p1[1] };
		p.BR = { p2[0], p2[1] };
		
		ImGui::SetClipboardText(p.ToString().c_str());
	}
	ImGui::SameLine();

	if (ImGui::Button("Paste") || paste)
	{
		const char* clip = ImGui::GetClipboardText();
		if (clip)
		{
			std::string str = clip;
			std::optional<Probe> op = Probe::ParseOne(str);
			if (op.has_value())
			{
				Probe p = op.value();
				p1[0] = p.TL.x; p1[1] = p.TL.y; p2[0] = p.BR.x; p2[1] = p.BR.y;
				strcpy_s(buf, size * sizeof(char), p.Name.c_str());
			}
			else
			{
				error_msg = "Invalid text pasted, use format {name}:{tl.x}:{tl.y}:{br.x}:{br.y}";
				error_time = stdclock::now();
				error_show = true;
			}
		}
	}

	if (error_show)
	{
		ImGui::TextUnformatted(error_msg.c_str());

		if (stdclock::now() - error_time > 2s)
			error_show = false;
	}

	ImGui::End();
	return false;
}



std::underlying_type_t<DrawGraphMode> operator&(DrawGraphMode a, DrawGraphMode b)
{
	using type = std::underlying_type_t<DrawGraphMode>;
	return (to<type>(a) & to<type>(b));
}
DrawGraphMode operator|(DrawGraphMode a, DrawGraphMode b)
{
	using type = std::underlying_type_t<DrawGraphMode>;
	return to<DrawGraphMode>(to<type>(a) | to<type>(b));
}

GraphRange DrawGraph(const std::vector<coord<float, 2>>& points, ImVec2 size, DrawGraphMode mode, const consumer<ImDrawList*>& custom_draw)
{
	// Create a canvas
	ImGui::BeginChild("GraphCanvas", size, true);

	ImDrawList* drawList = ImGui::GetWindowDrawList();

	ImVec2 canvasPos = ImGui::GetCursorScreenPos();
	ImVec2 canvasSize = ImGui::GetContentRegionAvail();

	// Background
	drawList->AddRectFilled(canvasPos,
		ImVec2(canvasPos.x + canvasSize.x, canvasPos.y + canvasSize.y),
		IM_COL32(30, 30, 30, 255));

	if (points.size() < 2)
	{
		ImGui::EndChild();
		return GraphRange{ {0, 0}, {0, 0} };
	}

	// Find bounds
	float minX = points[0].x, maxX = points[0].x;
	float minY = points[0].y, maxY = points[0].y;

	for (const auto& p : points)
	{
		minX = std::min(minX, p.x);
		maxX = std::max(maxX, p.x);
		minY = std::min(minY, p.y);
		maxY = std::max(maxY, p.y);
	}

	float rangeX = (maxX - minX) > 0 ? (maxX - minX) : 1.0f;
	float rangeY = (maxY - minY) > 0 ? (maxY - minY) : 1.0f;

	// Helper: transform point -> screen space
	auto ToScreen = [&](const coord<float, 2>& p) -> ImVec2
		{
			float nx = (p.x - minX) / rangeX;
			float ny = (p.y - minY) / rangeY;

			return ImVec2(
				canvasPos.x + nx * canvasSize.x,
				canvasPos.y + (1.0f - ny) * canvasSize.y // invert Y
			);
		};

	if ((mode & DrawGraphMode::LINE) != 0)
	{
		for (size_t i = 0; i < points.size() - 1; ++i)
		{
			ImVec2 p1 = ToScreen(points[i]);
			ImVec2 p2 = ToScreen(points[i + 1]);

			drawList->AddLine(p1, p2, IM_COL32(0, 255, 0, 255), 2.0f);
		}
	}

	if ((mode & DrawGraphMode::POINT) != 0)
	{
		for (const auto& p : points)
			drawList->AddCircleFilled(ToScreen(p), 3.0f, IM_COL32(255, 255, 255, 255));
	}

	if (custom_draw)
		custom_draw(drawList);

	ImGui::EndChild();

	return GraphRange{ { minX, maxX }, { minY, maxY } };
}