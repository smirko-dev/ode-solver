
################################################################################
# ODE Solver
################################################################################

ADD_LIBRARY(ode INTERFACE)

TARGET_SOURCES(ode INTERFACE
    Vector.h
    OdeFunction.h
    OdeSolver.h
    Euler.h
    MidPoint.h
    RungeKutta.h
    VelocityVerlet.h
)

TARGET_INCLUDE_DIRECTORIES(ode INTERFACE ${CMAKE_CURRENT_LIST_DIR})
