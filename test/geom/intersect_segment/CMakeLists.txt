################################################################################
# Copyright (c) 2013 David D. Marshall <ddmarsha@calpoly.edu>
#
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the Eclipse Public License v1.0
# which accompanies this distribution, and is available at
# http://www.eclipse.org/legal/epl-v10.html
#
# Contributors:
#    David D. Marshall - initial code and implementation
################################################################################

project(INTERSECT_SEGMENT_TEST)

# add the minimum_distance test
set (INTERSECT_SEGMENT_TEST_EXE IntersectSegmentTest)
set (INTERSECT_SEGMENT_TEST_SOURCE ${PROJECT_SOURCE_DIR}/intersect_segment_test.cpp)
set (INTERSECT_SEGMENT_TEST_HEADER ${PROJECT_SOURCE_DIR}/intersect_segment_triangle_test_suite.hpp
                                  ${PROJECT_SOURCE_DIR}/intersect_segment_surface_test_suite.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/constants/math.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/util/tolerance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/mutil/nls/newton_raphson_method.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/general/bounding_box.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/point/distance.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/curve/bezier.hpp
                                  ${ELI_SOURCE_DIR}/include/eli/geom/surface/bezier.hpp)


add_executable(${INTERSECT_SEGMENT_TEST_EXE} ${INTERSECT_SEGMENT_TEST_SOURCE} ${INTERSECT_SEGMENT_TEST_HEADER})
target_link_libraries(${INTERSECT_SEGMENT_TEST_EXE} ${_ELI_TEST_LIBRARIES_LIST})
add_test(NAME ${INTERSECT_SEGMENT_TEST_EXE} COMMAND $<TARGET_FILE:${INTERSECT_SEGMENT_TEST_EXE}> --text)
