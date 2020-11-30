//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>

#include <AnalysisTree/Detector.hpp>
#include <at_task/Task.h>
#include <memory>
#include <string>

class Rapidity : public UserFillTask {

public:
void Init(std::map<std::string, void *> &Map) override;
void Exec() override;
void Finish() override {}
boost::program_options::options_description GetBoostOptions() override;
void PreInit() override;
void PostFinish() override {
UserTask::PostFinish();
}

private:
int pdg_code_;

std::string tracks_branch_;

AnalysisTree::Particles *tracks_{nullptr};
AnalysisTree::Particles *rec_particles_{nullptr};
AnalysisTree::BranchConfig rec_particle_config_;

short y_cm_field_id_;

TASK_DEF(Rapidity, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
