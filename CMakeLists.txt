cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0015 NEW)

set (CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# search for pkg-config
include (FindPkgConfig)
if (NOT PKG_CONFIG_FOUND)
    message (FATAL_ERROR "pkg-config not found")
endif ()

pkg_check_modules (LIBPNG libpng16 REQUIRED)
if (NOT LIBPNG_FOUND)
    message(FATAL_ERROR "You don't seem to have libpng16 development libraries installed")
else ()
    include_directories (${LIBPNG_INCLUDE_DIRS})
    link_directories (${LIBPNG_LIBRARY_DIRS})
    link_libraries (${LIBPNG_LIBRARIES})
endif ()

set(appName reconstruction)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLM REQUIRED)

include_directories(${OPENGL_INCLUDE_DIRS})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories(${GLM_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/gl3w")
include_directories("${CMAKE_SOURCE_DIR}/eigen")
include_directories("${CMAKE_SOURCE_DIR}/pugixml")

link_directories(${OPENGL_LIBRARY_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})

link_libraries(dl)

set(imgui_files imgui/imconfig.h imgui/imgui.cpp imgui/imgui_demo.cpp imgui/imgui_draw.cpp imgui/imgui.h imgui/imgui_internal.h imgui/imgui_widgets.cpp imgui/imstb_rectpack.h imgui/imstb_textedit.h imgui/imstb_truetype.h imgui/imgui_impl_opengl3.h imgui/imgui_impl_opengl3.cpp imgui/imgui_impl_glut.h imgui/imgui_impl_glut.cpp BiharmonicSolver.h BiharmonicSolver.cpp)

set(gl3w_files gl3w/gl3w.c gl3w/GL/gl3w.h)

set(pugixml_files pugixml/pugixml.hpp pugixml/pugixml.cpp pugixml/pugixml.hpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_executable(${appName} ${imgui_files} ${gl3w_files} ${pugixml_files} main.cpp mesh.h mesh.cpp mesh_io.h mesh_io.cpp points_generator.h points_generator.cpp stb_image_write.h stb_image.h Image.cpp Image.h ScalarField.cpp ScalarField.h timing.h timing.cpp Quadtree.cpp Quadtree.h)

target_link_libraries(${appName} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
