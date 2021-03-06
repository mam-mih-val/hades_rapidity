//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>

#include <at_task/UserTask.h>
#include <at_task/Task.h>

#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Particle.hpp>
#include <AnalysisTree/BranchConfig.hpp>
#include <AnalysisTree/EventHeader.hpp>

#include <TH1F.h>
#include <TH2F.h>
#include <memory>
#include <string>

class Rapidity : public UserFillTask {

public:
  void UserInit(std::map<std::string, void *> &Map) override;
  void UserExec() override;
  void UserFinish() override {}
  boost::program_options::options_description GetBoostOptions() override;
  void PreInit() override;
  void PostFinish() override {
    UserTask::PostFinish();
  }

protected:
  virtual bool UseATI2() const override;

private:
  void ReadEfficiencyHistos();
  int pdg_code_;

  std::string tracks_branch_;
  std::string out_tracks_branch_;

  AnalysisTree::EventHeader *event_header_{nullptr};
  AnalysisTree::Particles *tracks_{nullptr};
  AnalysisTree::Particles *rec_particles_{nullptr};
  AnalysisTree::BranchConfig rec_particle_config_;

  TH1F* centrality_histo_{nullptr};

  std::string config_directory_;
  TFile* file_efficiency_protons_{nullptr};
  std::vector<TH2F*> efficiency_protons_;
  TFile* file_efficiency_pi_plus_{nullptr};
  std::vector<TH2F*> efficiency_pi_plus_;
  TFile* file_efficiency_pi_minus_{nullptr};
  std::vector<TH2F*> efficiency_pi_minus_;

TASK_DEF(Rapidity, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
