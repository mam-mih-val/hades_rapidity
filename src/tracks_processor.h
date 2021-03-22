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

class TracksProcessor : public UserFillTask {

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
//  virtual bool UseATI2() const override;

private:
  void ReadEfficiencyHistos();
  bool is_mc_=true;

  std::string out_tracks_branch_;

  ATI2::Branch *event_header_{nullptr};
  ATI2::Variable centrality_var_;

  ATI2::Branch *in_tracks_{nullptr};
  ATI2::Branch *out_tracks_{nullptr};
  ATI2::Variable out_ycm_var_;
  ATI2::Variable out_abs_ycm_var_;
  ATI2::Variable out_efficiency_var_;
  ATI2::Variable out_protons_rapidity_var_;
  ATI2::Variable out_pions_rapidity_var_;

  ATI2::Branch *in_sim_particles_{nullptr};
  ATI2::Branch *out_sim_particles_{nullptr};
  ATI2::Variable out_sim_ycm_var_;
  ATI2::Variable out_sim_abs_ycm_var_;
  ATI2::Variable out_sim_protons_rapidity_var_;
  ATI2::Variable out_sim_pions_rapidity_var_;

  std::string config_directory_;
  TFile* file_efficiency_protons_{nullptr};
  std::vector<TH2F*> efficiency_protons_;
  TFile* file_efficiency_pi_plus_{nullptr};
  std::vector<TH2F*> efficiency_pi_plus_;
  TFile* file_efficiency_pi_minus_{nullptr};
  std::vector<TH2F*> efficiency_pi_minus_;

TASK_DEF(TracksProcessor, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
