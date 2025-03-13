# Add the main build target.

set (target_name "${CMAKE_PROJECT_NAME}")
add_executable ("${target_name}" jstmap/main.cpp)
target_link_libraries ("${target_name}" PUBLIC seqan3::seqan3
                                               seqan::seqan2
                                               jstmap::global
                                               jstmap::create
                                               jstmap::index
                                               jstmap::linear
                                               jstmap::search
                                               jstmap::simulate
                                               # jstmap::view
                                               Threads::Threads
)
