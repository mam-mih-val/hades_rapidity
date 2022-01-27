#add_executable(rapidity tracks_processor.cc tracks_processor.h sim_particles_processor.cc sim_particles_processor.h)
add_executable(pre_process tracks_processor.cc tracks_processor.h)
target_link_libraries(pre_process PUBLIC ROOT::EG QnTools::Base at_task_main)