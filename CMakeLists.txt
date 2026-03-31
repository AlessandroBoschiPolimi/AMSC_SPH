cmake_minimum_required(VERSION 3.16)

project(AMSC_SPH LANGUAGES C CXX)

# ==================================================
# Options
# ==================================================
option(DISABLE_UI "Disable UI" OFF)
# ==================================================
# C++ Standard
# ==================================================
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==================================================
# Source Files
# ==================================================
file(GLOB_RECURSE CPP_SOURCES
    src/*.cpp
)
if (DISABLE_UI)
	list(REMOVE_ITEM CPP_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/src/Visualizers/ImGuiWidgets.cpp)
endif()


list(REMOVE_ITEM CPP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main_mpi.cpp
)

# ==================================================
# ImGui
# ==================================================
if (NOT DISABLE_UI)
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
endif()

# ==================================================
# GLAD
# ==================================================
if (NOT DISABLE_UI)
	add_library(glad STATIC
	    vendor/glad/src/glad.c
	)
	
	target_include_directories(glad PUBLIC
	    vendor/glad/include
	)
endif()

# ==================================================
# GLFW
# ==================================================
if (NOT DISABLE_UI)
	add_subdirectory(vendor/glfw)
endif()

# ==================================================
# OpenMP
# ==================================================
find_package(OpenMP REQUIRED)

# ==================================================
# MPI (optional)
# ==================================================
find_package(MPI REQUIRED)

# ==================================================
# Executable
# ==================================================
add_executable(${PROJECT_NAME}
    src/main.cpp
    ${CPP_SOURCES}
    ${IMGUI_SOURCES}
)
add_executable(${PROJECT_NAME}_MPI
     src/main_mpi.cpp
     ${CPP_SOURCES}
     ${IMGUI_SOURCES}
)

# ==================================================
# Includes
# ==================================================
if (DISABLE_UI)
	target_include_directories(${PROJECT_NAME} PRIVATE
	    src
	)
	target_include_directories(${PROJECT_NAME}_MPI PRIVATE
	    src
	)
else()
	target_include_directories(${PROJECT_NAME} PRIVATE
	    src
	    vendor/imgui
	    vendor/imgui/backends
	    vendor/glad/include
	    vendor/glfw/include
	)
	target_include_directories(${PROJECT_NAME}_MPI PRIVATE
	    src
	    vendor/imgui
	    vendor/imgui/backends
	    vendor/glad/include
	    vendor/glfw/include
	)
endif()
# ==================================================
# Link Libraries
# ==================================================

if (DISABLE_UI)
	target_link_libraries(${PROJECT_NAME}
	    OpenMP::OpenMP_CXX
	)
	
	target_link_libraries(${PROJECT_NAME}_MPI
	    MPI::MPI_CXX
	    OpenMP::OpenMP_CXX
	)
else()
	target_link_libraries(${PROJECT_NAME}
	    glad
	    glfw
	    OpenMP::OpenMP_CXX
	)
	
	target_link_libraries(${PROJECT_NAME}_MPI
	    glad
	    glfw
	    MPI::MPI_CXX
	    OpenMP::OpenMP_CXX
	)
endif()
# ==================================================
# Platform-specific OpenGL
# ==================================================
if (NOT DISABLE_UI)
	if (APPLE)
	    target_link_libraries(${PROJECT_NAME}
	        "-framework OpenGL"
	    )
	    target_link_libraries(${PROJECT_NAME}_MPI
	        "-framework OpenGL"
	    )
	elseif (WIN32)
	    target_link_libraries(${PROJECT_NAME}
	        opengl32
	    )
	    target_link_libraries(${PROJECT_NAME}_MPI
	        opengl32
	    )
	else ()
	    find_package(OpenGL REQUIRED)
	    find_package(X11 REQUIRED)
	
	    target_link_libraries(${PROJECT_NAME}
	        OpenGL::GL
	        X11::X11
	    )
	    target_link_libraries(${PROJECT_NAME}_MPI
	        OpenGL::GL
	        X11::X11
	    )
        endif()
endif()

# ==================================================
# Build Configuration Macros
# ==================================================
target_compile_definitions(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Debug>:MY_DEBUG>
    $<$<CONFIG:Release>:MY_RELEASE>
    $<$<CONFIG:Release>:MY_NDEBUG>
    $<$<CONFIG:RelWithDebInfo>:MY_DIST>
    $<$<CONFIG:RelWithDebInfo>:MY_DEBUG>
)
target_compile_definitions(${PROJECT_NAME}_MPI PRIVATE
    $<$<CONFIG:Debug>:MY_DEBUG>
    $<$<CONFIG:Release>:MY_RELEASE>
    $<$<CONFIG:Release>:MY_NDEBUG>
    $<$<CONFIG:RelWithDebInfo>:MY_DIST>
    $<$<CONFIG:RelWithDebInfo>:MY_DEBUG>
)
if(DISABLE_UI)
	target_compile_definitions(${PROJECT_NAME} PRIVATE DISABLE_UI)
	target_compile_definitions(${PROJECT_NAME}_MPI PRIVATE DISABLE_UI)
endif()
