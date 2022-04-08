
################################################################################
# ODE Solver
################################################################################

ADD_LIBRARY(OdeSolver INTERFACE)

TARGET_SOURCES(OdeSolver INTERFACE
    Vector.h
    OdeFunction.h
    OdeSolver.h
    Euler.h
    MidPoint.h
    RungeKutta.h
    VelocityVerlet.h
)

TARGET_INCLUDE_DIRECTORIES(OdeSolver INTERFACE ${CMAKE_CURRENT_LIST_DIR})
