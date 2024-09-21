add_library(filtering_kalman "kalman.c")
target_link_libraries(
    filtering_kalman
    PUBLIC
    matrix
    linear_algebra_decompositions)