//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>

#include <AnalysisTree/Detector.hpp>
#include <TH1F.h>
#include <TH2F.h>
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

  int GetCentralityClass(int multiplicity);
  void InitCentralityHisto();
  int pdg_code_;

  std::string tracks_branch_;
  std::string out_tracks_branch_;

  AnalysisTree::Particles *tracks_{nullptr};
  AnalysisTree::Particles *rec_particles_{nullptr};
  AnalysisTree::BranchConfig rec_particle_config_;

  short in_dca_xy_id_;
  short in_dca_z_id_;
  short in_chi2_id_;
  short out_charge_id_;
  short out_dca_xy_id_;
  short out_dca_z_id_;
  short out_efficiency_id_;
  short out_chi2_id_;

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
