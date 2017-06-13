set(HEADER_FILES Constants.h  AdlerFunction.h NumericalMethods.h RunAlpha.h IntegralMomentum.h experimentalData/ExperimentalData.h)
set(SOURCE_FILES NumericalMethods.cpp RunAlpha.cpp Constants.cpp AdlerFunction.cpp IntegralMomentum.cpp
        experimentalData/ExperimentalData.cpp experimentalData/aleph14_vpa.f90 experimentalData/num_const.f90 experimentalData/teubner.f)

add_library(program ${HEADER_FILES} ${SOURCE_FILES})

