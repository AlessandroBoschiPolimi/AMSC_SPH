cmake_minimum_required(VERSION 3.16)

project(AMSC_SPH LANGUAGES C CXX)

# ==================================================
# Options
# ==================================================
option(ENABLE_CUDA "Enable CUDA acceleration" ON)

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

file(GLOB_RECURSE CUDA_SOURCES
    src/*.cu
)

list(REMOVE_ITEM CPP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main_mpi.cpp
)

# ==================================================
# CUDA (Optional)
# ==================================================
set(HAS_CUDA FALSE)

if (ENABLE_CUDA)
    find_package(CUDAToolkit QUIET)

    if (CUDAToolkit_FOUND)

        message(STATUS "CUDA toolkit found. Enabling CUDA.")

        enable_language(CUDA)

        set(HAS_CUDA TRUE)

        set(CMAKE_CUDA_STANDARD 20)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)

        # Set architecture manually (REQUIRED to avoid auto-detect)
        set(CMAKE_CUDA_ARCHITECTURES 75)

    else()
        message(STATUS "CUDA toolkit not found. Building CPU-only version.")
    endif()
endif()

# ==================================================
# ImGui
# ==================================================
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

# ==================================================
# GLAD
# ==================================================
add_library(glad STATIC
    vendor/glad/src/glad.c
)

target_include_directories(glad PUBLIC
    vendor/glad/include
)

# ==================================================
# GLFW
# ==================================================
add_subdirectory(vendor/glfw)

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
if (HAS_CUDA)
    add_executable(${PROJECT_NAME}
	src/main.cpp
        ${CPP_SOURCES}
        ${CUDA_SOURCES}
        ${IMGUI_SOURCES}
    )
else()
    add_executable(${PROJECT_NAME}
	src/main.cpp
        ${CPP_SOURCES}
        ${IMGUI_SOURCES}
    )
endif()
add_executable(${PROJECT_NAME}_MPI
     src/main_mpi.cpp
     ${CPP_SOURCES}
     ${IMGUI_SOURCES}
)

# ==================================================
# Includes
# ==================================================
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
# ==================================================
# Link Libraries
# ==================================================
target_link_libraries(${PROJECT_NAME}
    glad
    glfw
    OpenMP::OpenMP_CXX
)

target_link_libraries(${PROJECT_NAME}_MPI
    glad
    glfw
    MPI::MPI_CXX
)

if (HAS_CUDA)
    target_link_libraries(${PROJECT_NAME}
        CUDA::cudart
    )

    set_target_properties(${PROJECT_NAME} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )

    target_compile_definitions(${PROJECT_NAME} PRIVATE HAS_CUDA)
endif()

# ==================================================
# Platform-specific OpenGL
# ==================================================
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

    add_compile_definitions(NOMINMAX)

    if (MSVC AND HAS_CUDA)
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
            $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
        )
    endif()
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
