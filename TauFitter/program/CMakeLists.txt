set(HEADER_FILES Constants.h  AdlerFunction.h NumericalMethods.h RunAlpha.h IntegralMomentum.h)
set(SOURCE_FILES NumericalMethods.cpp RunAlpha.cpp Constants.cpp AdlerFunction.cpp IntegralMomentum.cpp)


add_library(program ${HEADER_FILES} ${SOURCE_FILES})
