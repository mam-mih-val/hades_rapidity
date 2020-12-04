//
// Created by mikhail on 11/30/20.
//

#include "rapidity.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include <AnalysisTree/DataHeader.hpp>

TASK_IMPL(Rapidity)

boost::program_options::options_description Rapidity::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(GetName() + " options");
  desc.add_options()
      ("tracks-branch", value(&tracks_branch_)->default_value("mdc_vtx_tracks"), "Name of branch with tracks")
      ("out-branch", value(&out_tracks_branch_)->default_value("mdc_vtx_tracks_rapidity"), "Name of output branch")
      ("pdg-code", value(&pdg_code_)->default_value(2212), "PDG-Code");
  return desc;
}

void Rapidity::PreInit() {
  SetInputBranchNames({tracks_branch_});
  SetOutputBranchName({out_tracks_branch_});
}

void Rapidity::Init(std::map<std::string, void *> &Map) {
  tracks_ = static_cast<AnalysisTree::Particles*>(Map.at(tracks_branch_));

  rec_particle_config_ = AnalysisTree::BranchConfig(out_branch_, AnalysisTree::DetType::kParticle);
  rec_particle_config_.AddField<int>("charge");
  out_charge_id_ = rec_particle_config_.GetFieldId("charge");

  out_config_->AddBranchConfig(rec_particle_config_);
  rec_particles_ = new AnalysisTree::Particles;
  out_tree_->Branch(out_branch_.c_str(), &rec_particles_);
}

void Rapidity::Exec() {

  auto &particle_config = rec_particle_config_;
  rec_particles_->ClearChannels();

  TLorentzVector momentum;
  auto mass = TDatabasePDG::Instance()->GetParticle(pdg_code_)->Mass();
  for (int i_track = 0; i_track < tracks_->GetNumberOfChannels(); ++i_track) {
    auto track = tracks_->GetChannel(i_track);
    auto pid = track.GetPid();
    auto charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge();
    if( charge == 0 )
      continue;
    auto p = track.GetP();
    auto pT = track.GetPt();
    auto phi = track.GetPhi();
    auto particle = rec_particles_->AddChannel();
    particle->Init(particle_config);
    particle->SetMomentum3(track.GetMomentum3());
    particle->SetMass(mass);
    particle->SetPid(pid);
  }
}