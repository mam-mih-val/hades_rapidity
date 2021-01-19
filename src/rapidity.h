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
#include <TH3F.h>
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
  size_t WhatSector( double phi ){
    if( 0.0 < phi && phi < M_PI/3.0 )
      return 0;
    if( M_PI/3.0 < phi && phi < 2.0*M_PI/3.0 )
      return 1;
    if( 2*M_PI/3.0 < phi && phi < M_PI )
      return 2;
    if( -M_PI < phi && phi < -2*M_PI/3.0 )
      return 3;
    if( -2*M_PI/3.0 < phi && phi < -M_PI/3.0 )
      return 4;
    if( -M_PI/3.0 < phi && phi < 0.0 )
      return 5;
    return -1;
  }
  std::array<int, 6> CalcSectorsOccupancy();
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
  short out_protons_rapidity_;
  short out_pions_rapidity_;

  TH1F* centrality_histo_{nullptr};

  std::string config_directory_;
  TFile* file_occupancy_protons_{nullptr};
  TH2F* protons_slope_{nullptr};
  TH2F* protons_offset_{nullptr};

TASK_DEF(Rapidity, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
