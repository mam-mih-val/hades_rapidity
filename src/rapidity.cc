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
  SetInputBranchNames({tracks_branch_});
  SetOutputBranchName({out_tracks_branch_});
}

void Rapidity::Init(std::map<std::string, void *> &Map) {
  tracks_ = static_cast<AnalysisTree::Particles*>(Map.at(tracks_branch_));

  rec_particle_config_ = AnalysisTree::BranchConfig(out_branch_, AnalysisTree::DetType::kParticle);
  rec_particle_config_.AddField<int>("charge");
  rec_particle_config_.AddField<float>("chi2");
  rec_particle_config_.AddField<float>("dca_xy");
  rec_particle_config_.AddField<float>("dca_z");
  rec_particle_config_.AddField<float>("efficiency");
  rec_particle_config_.AddField<float>("protons_rapidity");
  rec_particle_config_.AddField<float>("pions_rapidity");

  out_charge_id_ = rec_particle_config_.GetFieldId("charge");
  out_dca_xy_id_ = rec_particle_config_.GetFieldId("dca_xy");
  out_dca_z_id_ = rec_particle_config_.GetFieldId("dca_z");
  out_efficiency_id_ = rec_particle_config_.GetFieldId("efficiency");
  out_protons_rapidity_ = rec_particle_config_.GetFieldId("protons_rapidity");
  out_pions_rapidity_ = rec_particle_config_.GetFieldId("pions_rapidity");

  in_chi2_id_ = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("chi2");
  in_dca_xy_id_ = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("dca_xy");
  in_dca_z_id_ = config_->GetBranchConfig( tracks_branch_ ).GetFieldId("dca_z");

  auto protons_file_name = config_directory_+"/efficiency_protons.root";
  file_efficiency_protons_ = TFile::Open(protons_file_name.c_str(), "read");
  auto pi_plus_file_name = config_directory_+"/efficiency_pi_plus.root";
  file_efficiency_pi_plus_ = TFile::Open(pi_plus_file_name.c_str(), "read");
  auto pi_minus_file_name = config_directory_+"/efficiency_pi_minus.root";
  file_efficiency_pi_minus_ = TFile::Open(pi_minus_file_name.c_str(), "read");
  if( file_efficiency_protons_ )
    file_efficiency_protons_->GetObject("efficiency_pT_y_n_tacks_sector", efficiency_protons_sectors_);
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

  out_config_->AddBranchConfig(rec_particle_config_);
  rec_particles_ = new AnalysisTree::Particles;
  out_tree_->Branch(out_branch_.c_str(), &rec_particles_);
  out_file_->cd();
}

void Rapidity::Exec() {

  auto &particle_config = rec_particle_config_;
  rec_particles_->ClearChannels();
  auto n_tracks = tracks_->GetNumberOfChannels();
  auto centrality_class = GetCentralityClass(n_tracks);
  TLorentzVector momentum;
  float y_beam_2{0.74};
  size_t n_recorded=0;
  for (int i_track = 0; i_track < tracks_->GetNumberOfChannels(); ++i_track) {
    auto track = tracks_->GetChannel(i_track);
    auto pid = track.GetPid();
    float charge{0.0};
    if( pid > 1e8 ){
      charge = (int)(pid / 1E+4) % (int)1e+3;
    }
    else {
      charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge()*3;
    }
    if( charge == 0 && pid != 0 ) {
      continue;
    }
    n_recorded++;
    auto pT = track.GetPt();
    auto p = track.GetP();
    auto pz = track.GetPz();
    auto y = track.GetRapidity();
    auto y_cm = y-y_beam_2;
    auto mass = track.GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
      else
        mass = track.GetMass();
    }
    auto particle = rec_particles_->AddChannel();
    particle->Init(particle_config);
    auto mass_protons = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    auto E_protons = sqrt( p*p + mass_protons*mass_protons );
    float y_protons = 0.5 * ( log( E_protons + pz ) - log( E_protons - pz ) );
    particle->SetField(y_protons, out_protons_rapidity_);
    auto mass_pions = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    auto E_pions = sqrt( p*p + mass_pions*mass_pions );
    float y_pions = 0.5 * ( log( E_pions + pz ) - log( E_pions - pz ) );
    particle->SetField(y_pions, out_pions_rapidity_);
    particle->SetMomentum3(track.GetMomentum3());
    particle->SetMass(mass);
    particle->SetPid(pid);
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
    if( efficiency > 0.5 )
      particle->SetField((float) 1.0 / efficiency, out_efficiency_id_);
    else
      particle->SetField(0.0f, out_efficiency_id_);
    try {
      auto chi2 = track.GetField<float>(in_chi2_id_);
      auto dca_xy = track.GetField<float>(in_dca_xy_id_);
      auto dca_z = track.GetField<float>(in_dca_z_id_);
      particle->SetField(chi2, out_chi2_id_);
      particle->SetField(dca_xy, out_dca_xy_id_);
      particle->SetField(dca_z, out_dca_z_id_);
    } catch (std::exception&) {}
  }
}

int Rapidity::GetCentralityClass(int multiplicity){
  if( !centrality_histo_ )
    InitCentralityHisto();
  auto bin = centrality_histo_->FindBin(multiplicity);
  return centrality_histo_->GetBinContent(bin) - 1;
}

void Rapidity::InitCentralityHisto(){
  Double_t xAxis[17] = {2, 3, 5, 7, 8, 11, 14, 17, 21, 26, 32, 38, 45, 54, 65, 102, 499};
  centrality_histo_ = new TH1F("TOFRPC_5pc_fixedCuts__1", "TOFRPC_5pc_fixedCuts", 16, xAxis);
  centrality_histo_->SetBinContent(1, 15);
  centrality_histo_->SetBinContent(2, 14);
  centrality_histo_->SetBinContent(3, 13);
  centrality_histo_->SetBinContent(4, 12);
  centrality_histo_->SetBinContent(5, 11);
  centrality_histo_->SetBinContent(6, 10);
  centrality_histo_->SetBinContent(7, 9);
  centrality_histo_->SetBinContent(8, 8);
  centrality_histo_->SetBinContent(9, 7);
  centrality_histo_->SetBinContent(10, 6);
  centrality_histo_->SetBinContent(11, 5);
  centrality_histo_->SetBinContent(12, 4);
  centrality_histo_->SetBinContent(13, 3);
  centrality_histo_->SetBinContent(14, 2);
  centrality_histo_->SetBinContent(15, 1);
  centrality_histo_->SetBinError(1, 1.8301);
  centrality_histo_->SetBinError(2, 5.29551);
  centrality_histo_->SetBinError(3, 5.07632);
  centrality_histo_->SetBinError(4, 4.68075);
  centrality_histo_->SetBinError(5, 4.83154);
  centrality_histo_->SetBinError(6, 5.30302);
  centrality_histo_->SetBinError(7, 4.84388);
  centrality_histo_->SetBinError(8, 5.20623);
  centrality_histo_->SetBinError(9, 4.82995);
  centrality_histo_->SetBinError(10, 5.22526);
  centrality_histo_->SetBinError(11, 4.99358);
  centrality_histo_->SetBinError(12, 4.83847);
  centrality_histo_->SetBinError(13, 4.94002);
  centrality_histo_->SetBinError(14, 5.02256);
  centrality_histo_->SetBinError(15, 5.08179);
  centrality_histo_->SetBinError(16, 0.00100254);
  centrality_histo_->SetEntries(16);
}