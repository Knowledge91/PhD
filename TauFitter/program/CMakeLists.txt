set(HEADER_FILES Constants.h  AdlerFunction.h)
set(SOURCE_FILES NumericalMethods.cpp NumericalMethods.h RunAlpha.cpp RunAlpha.h Constants.cpp AdlerFunction.cpp)


add_library(program ${HEADER_FILES} ${SOURCE_FILES})
