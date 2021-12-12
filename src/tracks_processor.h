//
// Created by mikhail on 11/30/20.
//

#ifndef HADES_RAPIDITY_SRC_RAPIDITY_H_
#define HADES_RAPIDITY_SRC_RAPIDITY_H_

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile3D.h>
#include <THnSparse.h>

#include <at_task/UserTask.h>
#include <at_task/Task.h>

#include <AnalysisTree/Matching.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Particle.hpp>
#include <AnalysisTree/BranchConfig.hpp>
#include <AnalysisTree/EventHeader.hpp>

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
  void LoopRecTracks();
  void LoopSimParticles();

private:
  void ReadEfficiencyHistos();
  double AngleDifference( double phi, double psi ){
    auto delta_phi = phi - psi;
    if( delta_phi < -M_PI )
      delta_phi+=2*M_PI;
    if( delta_phi > M_PI )
      delta_phi-=2*M_PI;
    return delta_phi;
  }
  bool is_mc_=true;

  // Input brancehs 
  ATI2::Branch *event_header_{nullptr};
  ATI2::Branch *in_tracks_{nullptr};
  ATI2::Branch *in_sim_particles_{nullptr};
  ATI2::Branch *in_sim_header_{nullptr};
  // Out event header
  ATI2::Branch *out_event_header_{nullptr};
  ATI2::Variable in_multiplicity_var_;
  ATI2::Variable out_centrality_var_;
  // Out Rec Particles 
  ATI2::Branch *out_tracks_{nullptr};
  ATI2::Variable out_theta_var_;
  ATI2::Variable out_ycm_var_;
  ATI2::Variable out_efficiency_var_;
  ATI2::Variable out_occ_weight_var_;
  // Out sim particles
  ATI2::Branch *out_sim_particles_{nullptr};
  ATI2::Variable out_sim_theta_var_;
  ATI2::Variable out_sim_ycm_var_;
  ATI2::Variable out_sim_is_charged_;
  // External centrality
  std::string str_centrality_file_;
  TFile* file_centrality_{nullptr};
  TH1F* h1_centrality_bins_{nullptr};
  // Protons efficiency
  std::string str_protons_efficiency_;
  TFile* file_efficiency_protons_{nullptr};
  std::vector<TH2F*> efficiency_protons_;
  // Pi pos efficiency
  std::string str_pi_plus_efficiency_;
  TFile* file_efficiency_pi_plus_{nullptr};
  std::vector<TH2F*> efficiency_pi_plus_;
  // Pi neg efficiency
  std::string str_pi_minus_efficiency_;
  TFile* file_efficiency_pi_minus_{nullptr};
  std::vector<TH2F*> efficiency_pi_minus_;
  // Efficiency for occupancy correction based on pairs approach
  std::string str_pair_efficiency_;
  TFile*file_pair_efficiency_{nullptr};
  TH3F*h3_2212_efficiency_centrality_phi_theta_;
  TH3F*h3_all_efficiency_centrality_phi_theta_;
  THnSparseF* hn_efficiency_pairs_centrality_phi_theta_;
TASK_DEF(TracksProcessor, 0)
};

#endif // HADES_RAPIDITY_SRC_RAPIDITY_H_
