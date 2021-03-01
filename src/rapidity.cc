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
      ("config-directory", value(&config_directory_)->default_value("../../efficiency_files/"), "Path to directory with efficiency")
      ("pdg-code", value(&pdg_code_)->default_value(0), "PDG-Code");
  return desc;
}

void Rapidity::PreInit() {
  SetInputBranchNames({tracks_branch_, "event_header"});
  SetOutputBranchName({out_tracks_branch_});
}

void Rapidity::UserInit(std::map<std::string, void *> &Map) {
  tracks_ = static_cast<AnalysisTree::Particles*>(Map.at(tracks_branch_));
  event_header_ = static_cast<AnalysisTree::EventHeader*>(Map.at("event_header"));
  ReadEfficiencyHistos();

  rec_particle_config_ = AnalysisTree::BranchConfig(out_branch_, AnalysisTree::DetType::kParticle);
  rec_particle_config_.AddField<bool>("is_primary");
  rec_particle_config_.AddField<bool>("is_pion");
  rec_particle_config_.AddField<float>("charge");
  rec_particle_config_.AddField<float>("ycm");
  rec_particle_config_.AddField<float>("chi2");
  rec_particle_config_.AddField<float>("dca_xy");
  rec_particle_config_.AddField<float>("dca_z");
  rec_particle_config_.AddField<float>("efficiency");
  rec_particle_config_.AddField<float>("protons_rapidity");
  rec_particle_config_.AddField<float>("pions_rapidity");

  out_config_->AddBranchConfig(rec_particle_config_);
  rec_particles_ = new AnalysisTree::Particles;
  out_tree_->Branch(out_branch_.c_str(), &rec_particles_);
  out_file_->cd();
}

void Rapidity::UserExec() {
  auto &particle_config = rec_particle_config_;
  rec_particles_->ClearChannels();
  auto n_tracks = tracks_->GetNumberOfChannels();
  auto centrality = event_header_->GetField<float>(
      config_->GetBranchConfig( "event_header" ).GetFieldId("selected_tof_rpc_hits_centrality"));
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  float y_beam_2 = data_header_->GetBeamRapidity();
  size_t n_recorded=0;

  auto out_is_primary_id = rec_particle_config_.GetFieldId("is_primary");
  auto out_is_pion = rec_particle_config_.GetFieldId("is_pion");
  auto out_charge_id = rec_particle_config_.GetFieldId("charge");
  auto out_y_cm_id = rec_particle_config_.GetFieldId("ycm");
  auto out_chi2_id = rec_particle_config_.GetFieldId("chi2");
  auto out_dca_xy_id = rec_particle_config_.GetFieldId("dca_xy");
  auto out_dca_z_id = rec_particle_config_.GetFieldId("dca_z");
  auto out_efficiency_id = rec_particle_config_.GetFieldId("efficiency");
  auto out_protons_rapidity = rec_particle_config_.GetFieldId("protons_rapidity");
  auto out_pions_rapidity = rec_particle_config_.GetFieldId("pions_rapidity");

  auto in_is_primary = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("is_primary");
  auto in_geant_pid = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("geant_pid");
  auto in_chi2_id = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("chi2");
  auto in_dca_xy_id = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("dca_xy");
  auto in_dca_z_id = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("dca_z");

  for (int i_track = 0; i_track < tracks_->GetNumberOfChannels(); ++i_track) {
    auto track = tracks_->GetChannel(i_track);
    auto particle = rec_particles_->AddChannel();
    particle->Init(particle_config);

    auto pid = track.GetPid();
    auto geant_code = track.GetField<int>(in_geant_pid);
    float charge{0};
    if( pid > 1e8 ){
      charge = (int) (pid / 1E+4) % (int)1e+3;
    }
    else {
      charge = std::roundf(TDatabasePDG::Instance()->GetParticle(pid)->Charge()/3);
    }
    if( charge == 0 && pid != 0 ) {
      continue;
    }
    particle->SetField(charge, out_charge_id);
    n_recorded++;
    auto pT = track.GetPt();
    auto p = track.GetP();
    auto pz = track.GetPz();
    auto mass = track.GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
      else
        mass = track.GetMass();
    }
    float y = track.GetRapidityByMass(mass);
    float y_cm = y-y_beam_2;
    auto mass_protons = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    auto E_protons = sqrt( p*p + mass_protons*mass_protons );
    float y_protons = 0.5 * ( log( E_protons + pz ) - log( E_protons - pz ) );
    particle->SetField(y_protons, out_protons_rapidity);
    auto mass_pions = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    auto E_pions = sqrt( p*p + mass_pions*mass_pions );
    float y_pions = 0.5 * ( log( E_pions + pz ) - log( E_pions - pz ) );
    particle->SetField(y_pions, out_pions_rapidity);
    particle->SetMomentum3(track.GetMomentum3());
    particle->SetMass(mass);
    particle->SetPid(pid);
    particle->SetField(y_cm, out_y_cm_id);
    auto is_pion = abs(pid) == 211;
    particle->SetField( is_pion, out_is_pion );
    TH2F* efficiency_histogram{nullptr};
    float efficiency{1.0};
    try{
      switch (pid) {
      case 211:
        efficiency_histogram = efficiency_pi_plus_.at(centrality_class);
        break;
      case -211:
        efficiency_histogram = efficiency_pi_minus_.at(centrality_class);
        break;
      case 2212:
        efficiency_histogram = efficiency_protons_.at(centrality_class);
        break;
      case 0:
        efficiency = 0.0;
      }
      if (efficiency_histogram) {
        auto bin_y = efficiency_histogram->GetXaxis()->FindBin(y_cm);
        auto bin_pT = efficiency_histogram->GetYaxis()->FindBin(pT);
        efficiency = efficiency_histogram->GetBinContent(bin_y, bin_pT);
      }
    } catch (std::exception&) {}
    if( efficiency > 0.1 )
      particle->SetField((float) 1.0 / efficiency, out_efficiency_id);
    else
      particle->SetField(0.0f, out_efficiency_id);
    try {
      auto is_primary = track.GetField<bool>(in_is_primary);
      particle->SetField(is_primary, out_is_primary_id);
    } catch (std::exception&) {}
    try {
      auto chi2 = track.GetField<float>(in_chi2_id);
      auto dca_xy = track.GetField<float>(in_dca_xy_id);
      auto dca_z = track.GetField<float>(in_dca_z_id);
      particle->SetField(chi2, out_chi2_id);
      particle->SetField(dca_xy, out_dca_xy_id);
      particle->SetField(dca_z, out_dca_z_id);
    } catch (std::exception&) {}
  }
}

void Rapidity::ReadEfficiencyHistos(){
  auto protons_file_name = config_directory_+"/efficiency_protons_botvina.root";
  file_efficiency_protons_ = TFile::Open(protons_file_name.c_str(), "read");
  auto pi_plus_file_name = config_directory_+"/efficiency_pi_plus.root";
  file_efficiency_pi_plus_ = TFile::Open(pi_plus_file_name.c_str(), "read");
  auto pi_minus_file_name = config_directory_+"/efficiency_pi_minus.root";
  file_efficiency_pi_minus_ = TFile::Open(pi_minus_file_name.c_str(), "read");
  int p=2;
  while(p<40){
    efficiency_protons_.emplace_back();
    efficiency_pi_plus_.emplace_back();
    efficiency_pi_minus_.emplace_back();
    std::string name = "efficiency_"+std::to_string(p);
    if( file_efficiency_protons_ )
      file_efficiency_protons_->GetObject(name.c_str(), efficiency_protons_.back());
    if( file_efficiency_pi_plus_ )
      file_efficiency_pi_plus_->GetObject(name.c_str(), efficiency_pi_plus_.back());
    if(file_efficiency_pi_minus_)
      file_efficiency_pi_minus_->GetObject(name.c_str(), efficiency_pi_minus_.back());
    p+=5;
  }
}