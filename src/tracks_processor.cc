//
// Created by mikhail on 11/30/20.
//

#include "tracks_processor.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"

#include <AnalysisTree/DataHeader.hpp>
#include <random>

TASK_IMPL(TracksProcessor)

boost::program_options::options_description TracksProcessor::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(GetName() + " options");
  desc.add_options()
      ("protons-efficiency", value(&str_protons_efficiency_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to file with proton's efficiency")
      ("centrality-file", value(&str_centrality_file_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to file with proton's efficiency")
      ("pair-efficiency", value(&str_pair_efficiency_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to file with proton's efficiency")
      ("pi-plus-efficiency", value(&str_pi_plus_efficiency_)->default_value("../../efficiency_files/pi_plus_efficiency.root"), "Path to file with pi-plus' efficiency")
      ("pi-minus-efficiency", value(&str_pi_minus_efficiency_)->default_value("../../efficiency_files/pi_minus_efficiency.root"), "Path to file with pi-minus' efficiency");
  return desc;
}

void TracksProcessor::PreInit() {}

void TracksProcessor::UserInit(std::map<std::string, void *> &Map) {
  auto y_cm = data_header_->GetBeamRapidity();
  in_tracks_ = GetInBranch("mdc_vtx_tracks");
  event_header_ = GetInBranch("event_header");

  out_tracks_ = NewBranch( "mdc_vtx_tracks_extra", PARTICLES );

  out_tracks_->CloneVariables(in_tracks_->GetConfig());

  out_ycm_var_ = out_tracks_->NewVariable( "ycm", FLOAT );
  out_theta_var_ = out_tracks_->NewVariable( "theta", FLOAT );
  out_efficiency_var_ = out_tracks_->NewVariable( "efficiency", FLOAT );
  out_occ_weight_var_ = out_tracks_->NewVariable( "occ_weight", FLOAT );

  out_event_header_ = NewBranch( "event_header_extra", EVENT_HEADER );
  out_centrality_var_ = out_event_header_->NewVariable("selected_tof_rpc_hits_centrality", FLOAT);

  try {
    in_sim_particles_ = GetInBranch("sim_tracks");
    in_sim_header_ = GetInBranch("sim_header");
    is_mc_=true;
  } catch (std::exception&) {
    is_mc_=false;
  }
  if( is_mc_ ) {
    out_sim_particles_ = NewBranch( "sim_particles_extra", PARTICLES );
    out_sim_particles_ ->CloneVariables( in_sim_particles_->GetConfig() );
    out_sim_theta_var_ = out_sim_particles_->NewVariable( "theta", FLOAT );
    out_sim_ycm_var_ = out_sim_particles_->NewVariable( "ycm", FLOAT );
    out_sim_is_charged_ = out_sim_particles_->NewVariable( "is_charged", BOOLEAN );
  }
  ReadEfficiencyHistos();
  out_file_->cd();
}

void TracksProcessor::UserExec() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if( h1_centrality_bins_ ){
    auto tof_rpc_hits = (*event_header_)[in_multiplicity_var_].GetInt();
    auto mult_bin = h1_centrality_bins_->GetXaxis()->FindBin(tof_rpc_hits);
    auto c_class = h1_centrality_bins_->GetBinContent(mult_bin);
    centrality = 5.0f * c_class - 2.5f;
    (*out_event_header_)[out_centrality_var_].SetVal(centrality);
  }
  this->LoopRecTracks();
  if( is_mc_ ){
    this->LoopSimParticles();
  }
}

void TracksProcessor::ReadEfficiencyHistos(){
  file_efficiency_protons_ = TFile::Open(str_protons_efficiency_.c_str(), "read");
  file_efficiency_pi_plus_ = TFile::Open(str_pi_plus_efficiency_.c_str(), "read");
  file_efficiency_pi_minus_ = TFile::Open(str_pi_minus_efficiency_.c_str(), "read");
  file_centrality_ = TFile::Open(str_centrality_file_.c_str(), "read");
  if(file_centrality_)
    file_centrality_->GetObject("Centrality/TOFRPC_5pc_fixedCuts",h1_centrality_bins_);
  // Initializing efficiency for occupancy correction
  file_pair_efficiency_ = TFile::Open(str_pair_efficiency_.c_str(), "read" );
  if(file_pair_efficiency_){
    file_pair_efficiency_->GetObject("h3_pid_efficiency_centrality_pT_theta_",
                                     h3_2212_efficiency_centrality_phi_theta_);
    file_pair_efficiency_->GetObject("h3_all_efficiency_centrality_pT_theta_",
                                     h3_all_efficiency_centrality_pT_theta_);
    file_pair_efficiency_->GetObject("hn_pair_efficiency_centrality_pT_phi_theta_",
                                     hn_efficiency_pairs_centrality_phi_theta_);
  }

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
void TracksProcessor::LoopRecTracks() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if( h1_centrality_bins_){
    auto tof_rpc_hits = (*event_header_)[GetVar("event_header/selected_tof_rpc_hits")].GetInt();
    auto mult_bin = h1_centrality_bins_->GetXaxis()->FindBin(tof_rpc_hits);
    auto c_class = h1_centrality_bins_->GetBinContent(mult_bin);
    centrality = 5.0*c_class-2.5;
  }
  auto psi_rp = 0.0;
  if( is_mc_ ){
    psi_rp = (*in_sim_header_)[GetVar( "sim_header/reaction_plane" )];
  }
  float y_beam = data_header_->GetBeamRapidity();
  out_tracks_->ClearChannels();

  for ( size_t idx1=0; idx1 < in_tracks_->size(); idx1++ ) {
    auto in_track = (*in_tracks_)[idx1];
    auto pid = in_track.DataT<Particle>()->GetPid();
    auto mass = in_track.DataT<Particle>()->GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    }
    auto mom4 = in_track.DataT<Particle>()->Get4MomentumByMass(mass);
    auto pT = mom4.Pt();
    auto p = mom4.P();
    auto pz = mom4.Pz();
    auto theta = mom4.Theta();
    float y = mom4.Rapidity();
    float y_cm = y - y_beam;
    TH2F* efficiency_histogram{nullptr};
    float efficiency{1.0};
    if( pid == 211 ){
      efficiency_histogram = centrality_class < efficiency_pi_plus_.size() ? efficiency_pi_plus_.at(centrality_class) : nullptr;
    }
    if( pid == -211 ){
      efficiency_histogram = centrality_class < efficiency_pi_minus_.size() ? efficiency_pi_minus_.at(centrality_class) : nullptr;
    }
    if( pid == 2212 ){
      efficiency_histogram = centrality_class < efficiency_protons_.size() ? efficiency_protons_.at(centrality_class) : nullptr;
    }
    if (efficiency_histogram) {
      auto bin_y = efficiency_histogram->GetXaxis()->FindBin(y_cm);
      auto bin_pT = efficiency_histogram->GetYaxis()->FindBin(pT);
      efficiency = efficiency_histogram->GetBinContent(bin_y, bin_pT);
    }
    auto occupancy_weight = 0.0;
    if( pid == 2212 ){
      double phi = mom4.Phi();
      double sector1 = floor(phi / (M_PI/3.0));
      double sector_center = sector1 *(M_PI/3.0) + M_PI/6.0;
      double sector_phi1 = AngleDifference(phi,sector_center);
      auto pT_theta_efficiency = h3_2212_efficiency_centrality_phi_theta_->Interpolate(centrality, mom4.Pt(), mom4.Theta());
      auto conditional_efficiency = pT_theta_efficiency;
      for( size_t idx2=0; idx2 < in_tracks_->size(); ++idx2){
        if(idx1==idx2)
          continue;
        auto in_track2 = (*in_tracks_)[idx2];
        auto pid2 = in_track2.DataT<Particle>()->GetPid();
        auto mass2 = in_track2.DataT<Particle>()->GetMass();
        if( pid2 != 0 ) {
          if( TDatabasePDG::Instance()->GetParticle(pid2) )
            mass2 = TDatabasePDG::Instance()->GetParticle(pid2)->Mass();
        }
        auto mom2 = in_track2.DataT<Particle>()->Get4MomentumByMass(mass2);
        double phi2 = mom2.Phi();
        double sector2 = floor(phi2 / (M_PI/3.0));
        if( fabs(sector1 - sector2) > 0.001 )
          continue;
        double sector_phi2 = AngleDifference(phi2,sector_center);
        auto pT_theta_efficiency2 =
            h3_all_efficiency_centrality_pT_theta_->Interpolate(centrality, mom2.Pt(), mom2.Theta());
        if(pT_theta_efficiency2 < 0.05 )
          continue;
        auto c_bin = hn_efficiency_pairs_centrality_phi_theta_->GetAxis(0)->FindBin(centrality);
        auto pT1_bin = hn_efficiency_pairs_centrality_phi_theta_->GetAxis(1)->FindBin(mom4.Pt());
        auto phi1_bin = hn_efficiency_pairs_centrality_phi_theta_->GetAxis(2)->FindBin(sector_phi1);
        auto theta1_bin = hn_efficiency_pairs_centrality_phi_theta_->GetAxis(3)->FindBin(mom4.Theta());
        auto phi2_bin = hn_efficiency_pairs_centrality_phi_theta_->GetAxis(4)->FindBin(sector_phi2);
        auto theta2_bin = hn_efficiency_pairs_centrality_phi_theta_->GetAxis(5)->FindBin(mom2.Theta());
        int index[] = {c_bin, pT1_bin, phi1_bin, theta1_bin, phi2_bin, theta2_bin};
        auto pair_efficiency = hn_efficiency_pairs_centrality_phi_theta_->GetBinContent( index );
        if( pair_efficiency < 0.05 )
          continue;
        conditional_efficiency*= pair_efficiency / pT_theta_efficiency2 / pT_theta_efficiency;
      }
      occupancy_weight = conditional_efficiency > 0.1 ? 1.0/conditional_efficiency : 0.0;
      if(pT_theta_efficiency < 0.05 )
        occupancy_weight = 0.0;
    }
    auto out_particle = out_tracks_->NewChannel();
    out_particle.CopyContents(in_track);
    out_particle[out_ycm_var_] = y_cm;
    out_particle[out_theta_var_] = (float) theta;
    out_particle[out_efficiency_var_] = efficiency > 0.1f ? 1.0f / efficiency : 0.0f;
    out_particle[out_occ_weight_var_] = occupancy_weight;
    out_particle.DataT<Particle>()->SetMass(mass);
  }
}
void TracksProcessor::LoopSimParticles() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if( h1_centrality_bins_){
    auto tof_rpc_hits = (*event_header_)[GetVar("event_header/selected_tof_rpc_hits")].GetInt();
    auto mult_bin = h1_centrality_bins_->GetXaxis()->FindBin(tof_rpc_hits);
    auto c_class = h1_centrality_bins_->GetBinContent(mult_bin);
    centrality = 5.0f * c_class - 2.5;
  }
  float y_beam = data_header_->GetBeamRapidity();

  out_sim_particles_->ClearChannels();
  for ( auto in_particle : in_sim_particles_->Loop() ) {
    auto out_particle = out_sim_particles_->NewChannel();
    out_particle.CopyContents(in_particle);
    auto pid = in_particle.DataT<Particle>()->GetPid();
    auto mass = in_particle.DataT<Particle>()->GetMass();
    auto mom4 = in_particle.DataT<Particle>()->Get4MomentumByMass(mass);
    auto pT = mom4.Pt();
    auto theta = mom4.Theta();
    auto p = mom4.P();
    auto pz = mom4.Pz();
    float y = mom4.Rapidity();
    float y_cm = y- y_beam;
    int charge=0;
    if( TDatabasePDG::Instance()->GetParticle(pid) )
      charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge() / 3;
    else{
      charge = (int)(pid / 1E+4) % (int)1e+3;
    }
    out_particle[out_sim_theta_var_] = (float) mom4.Theta();
    out_particle[out_sim_ycm_var_] = y_cm;
    out_particle[out_sim_is_charged_] = charge != 0;
    out_particle.DataT<Particle>()->SetMass(mass);
  }
}
