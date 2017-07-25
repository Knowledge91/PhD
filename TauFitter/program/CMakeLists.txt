set(HEADER_FILES Constants.h  AdlerFunction.h NumericalMethods.h RunAlpha.h Theory.h experimentalData/ExperimentalData.h Minuit2/RosenBrockFCN.h Chisquared.h)
set(SOURCE_FILES NumericalMethods.cpp RunAlpha.cpp Constants.cpp AdlerFunction.cpp Theory.cpp
        experimentalData/ExperimentalData.cpp experimentalData/aleph14_vpa.f90 experimentalData/num_const.f90 experimentalData/teubner.f Chisquared.cpp)

add_library(program ${HEADER_FILES} ${SOURCE_FILES})

