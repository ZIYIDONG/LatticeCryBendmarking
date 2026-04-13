# 清理旧缓存（因为改了 CMakeLists.txt）
rm -rf build
mkdir build && cd build

cmake ..
cmake --build . -j

# 然后挨个跑（或者只跑你关心的那个）
./LatticeCryBenchmarking
./test_powersof
./test_powersof_modswitch
./test_frd
./debug_frd
./bench_matops







#设置cmake版本
cmake_minimum_required(VERSION 3.22)

#项目名字
project(LatticeCryBenchmarking)

#设置编译版本
set(CMAKE_CXX_STANDARD 20)

#声明头文件路径
set(INC_DIR ${CMAKE_SOURCE_DIR}/include)

#声明链接库路径
set(LINK_DIR ${CMAKE_SOURCE_DIR}/lib)


#编译文件
add_executable(LatticeCryBenchmarking
        ./src/main.cpp
        ./src/mp12deltrapgen.cpp
)

# 添加头文件搜索路径
target_include_directories(LatticeCryBenchmarking PRIVATE ${INC_DIR})

# 添加库文件搜索路径
target_link_directories(LatticeCryBenchmarking PRIVATE ${LINK_DIR})


# 链接库
#target_link_libraries(LatticeCryBenchmarking PRIVATE
#        pbc
#        gmp
#        ssl
#        crypto
#)

#生成调试信息
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")