set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -pedantic -std=c99 -ldl")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -D_DEBUG -g -ggdb -ldl")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -D_NDEBUG -O3 -g0 -ldl")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DNDEBUG -D_NDEBUG -O3 -g -ldl")
