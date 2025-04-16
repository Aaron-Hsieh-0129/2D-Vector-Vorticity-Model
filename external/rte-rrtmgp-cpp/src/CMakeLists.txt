FILE(GLOB sourcefiles "../src/*.cpp")
# Include Radiation_solver.cpp from src_test
list(APPEND sourcefiles "../src_test/Radiation_solver.cpp")
add_library(rte_rrtmgp STATIC ${sourcefiles} Aerosol_optics.cpp ../include/Aerosol_optics.h)
target_link_libraries(rte_rrtmgp PRIVATE rte_rrtmgp_kernels)
target_include_directories(rte_rrtmgp PUBLIC
    "../include"
    "../include_test"
    "../include_rt"
    "../include_rt_kernels"
    SYSTEM ${INCLUDE_DIRS}
)
