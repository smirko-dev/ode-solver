
################################################################################
# ODE Solver
################################################################################

ADD_LIBRARY(ode INTERFACE)

TARGET_SOURCES(ode INTERFACE
    ode/Vector.h
    ode/Function.h
    ode/Solver.h
    ode/Euler.h
    ode/MidPoint.h
    ode/RungeKutta.h
    ode/VelocityVerlet.h
)

TARGET_INCLUDE_DIRECTORIES(ode INTERFACE ${CMAKE_CURRENT_LIST_DIR})
