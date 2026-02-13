cmake_minimum_required(VERSION 3.16)
project(AMSC_SPH LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# --------------------------------------------------
# Source files
# --------------------------------------------------
file(GLOB_RECURSE SRC_FILES
	src/*.cpp)

# --------------------------------------------------
# ImGui sources
# --------------------------------------------------
set(IMGUI_DIR vendor/imgui)

set(IMGUI_SOURCES
	${IMGUI_DIR}/imgui.cpp
	${IMGUI_DIR}/imgui_draw.cpp
	${IMGUI_DIR}/imgui_tables.cpp
	${IMGUI_DIR}/imgui_widgets.cpp
	${IMGUI_DIR}/imgui_demo.cpp
	${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
	${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)

# --------------------------------------------------
# GLAD
# --------------------------------------------------
add_library(glad STATIC
	vendor/glad/src/glad.c
)

target_include_directories(glad PUBLIC
	vendor/glad/include
)

# --------------------------------------------------
# GLFW
# --------------------------------------------------
add_subdirectory(vendor/glfw)


# --------------------------------------------------
# Virtual Terminal
# --------------------------------------------------
# file(GLOB_RECURSE VIRTUALTERMINAL_SOURCES
#     vendor/VirtualTerminal/Src/*.cpp
# )
# 
# add_library(virtualterminal STATIC
#     ${VIRTUALTERMINAL_SOURCES}
# )
# 
# target_include_directories(virtualterminal PUBLIC
#     vendor/VirtualTerminal/Src
# )


# --------------------------------------------------
# Executable
# --------------------------------------------------
add_executable(${PROJECT_NAME}
	${SRC_FILES}
	${IMGUI_SOURCES}
)

# --------------------------------------------------
# Includes
# --------------------------------------------------
target_include_directories(${PROJECT_NAME} PRIVATE
	src
	vendor/imgui
	vendor/imgui/backends
	vendor/glad/include
	vendor/glfw/include
	# vendor/VirtualTerminal/Src
)

# --------------------------------------------------
# Link libraries
# --------------------------------------------------
target_link_libraries(${PROJECT_NAME}
	glad
	glfw
	OpenMP::OpenMP_CXX
	# virtualterminal
)

# --------------------------------------------------
# OpenMP
# --------------------------------------------------
find_package(OpenMP REQUIRED)

# --------------------------------------------------
# Platform-specific OpenGL
# --------------------------------------------------
if (APPLE)
	target_link_libraries(${PROJECT_NAME}
		"-framework OpenGL"
	)
elseif (WIN32)
	target_link_libraries(${PROJECT_NAME}
		opengl32
	)
else ()
	find_package(OpenGL REQUIRED)
	find_package(X11 REQUIRED)

	target_link_libraries(${PROJECT_NAME}
		OpenGL::GL
		X11::X11
	)
endif()


target_compile_definitions(AMSC_SPH PRIVATE
	$<$<CONFIG:Debug>:MY_DEBUG>
	$<$<CONFIG:Release>:MY_RELEASE>
	$<$<CONFIG:Release>:MY_NDEBUG>
	$<$<CONFIG:RelWithDebInfo>:MY_DIST>
	$<$<CONFIG:RelWithDebInfo>:MY_DEBUG>
)