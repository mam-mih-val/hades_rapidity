add_executable(rapidity rapidity.cc rapidity.h)
target_link_libraries(rapidity PUBLIC ROOT::EG at_task_main)