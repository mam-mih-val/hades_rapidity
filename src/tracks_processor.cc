//
// Created by mikhail on 11/30/20.
//

#include "tracks_processor.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"

#include <AnalysisTree/DataHeader.hpp>

TASK_IMPL(TracksProcessor)

boost::program_options::options_description TracksProcessor::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(GetName() + " options");
  desc.add_options()
      ("efficiency-file", value(&config_directory_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to efficiency file");
  return desc;
}

void TracksProcessor::PreInit() {}

void TracksProcessor::UserInit(std::map<std::string, void *> &Map) {
  in_tracks_ = GetInBranch("mdc_vtx_tracks");
  event_header_ = GetInBranch("event_header");
  centrality_var_ = GetVar( "event_header/selected_tof_rpc_hits_centrality" );
  out_tracks_ = NewBranch( "mdc_vtx_tracks_extra", PARTICLES );
  out_tracks_->CloneVariables(in_tracks_->GetConfig());
  out_ycm_var_ = out_tracks_->NewVariable( "ycm", FLOAT );
  out_abs_ycm_var_ = out_tracks_->NewVariable( "abs_ycm", FLOAT );
  out_efficiency_var_ = out_tracks_->NewVariable( "efficiency", FLOAT );
  out_protons_rapidity_var_ = out_tracks_->NewVariable( "protons_rapidity", FLOAT );
  out_pions_rapidity_var_ = out_tracks_->NewVariable( "pions_rapidity", FLOAT );
  try {
    in_sim_particles_ = GetInBranch("sim_tracks");
    is_mc_=true;
  } catch (std::exception&) {
    is_mc_=false;
  }
  if( is_mc_ ) {
    out_sim_particles_ = NewBranch( "sim_particles_extra", PARTICLES );
    out_sim_particles_ ->CloneVariables( in_sim_particles_->GetConfig() );
    out_sim_ycm_var_ = out_sim_particles_->NewVariable( "ycm", FLOAT );
    out_sim_abs_ycm_var_ = out_sim_particles_->NewVariable( "abs_ycm", FLOAT );
    out_sim_protons_rapidity_var_ = out_sim_particles_->NewVariable( "protons_rapidity", FLOAT );
    out_sim_pions_rapidity_var_ = out_sim_particles_->NewVariable( "pions_rapidity", FLOAT );
    out_sim_is_charged_ = out_sim_particles_->NewVariable( "is_charged", BOOLEAN );
  }
  ReadEfficiencyHistos();
  out_file_->cd();
}

void TracksProcessor::UserExec() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[centrality_var_].GetVal();
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  float y_beam = data_header_->GetBeamRapidity();
  out_tracks_->ClearChannels();
  for ( auto in_track : in_tracks_->Loop() ) {
    auto pid = in_track.DataT<Particle>()->GetPid();
    int charge=0;
    if( TDatabasePDG::Instance()->GetParticle(pid) )
      charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge() / 3;
    else{
      charge = (int)(pid / 1E+4) % (int)1e+3;
    }
    if(charge==0 && pid == 0)
      continue;
    auto mass = in_track.DataT<Particle>()->GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    }
    auto mom4 = in_track.DataT<Particle>()->Get4MomentumByMass(mass);
    auto pT = mom4.Pt();
    auto p = in_track.DataT<Particle>()->GetP();
    auto pz = in_track.DataT<Particle>()->GetPz();
    float y = mom4.Rapidity();
    float y_cm = y- y_beam;
    auto mass_protons = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    auto E_protons = sqrt( p*p + mass_protons*mass_protons );
    float y_protons = 0.5 * ( log( E_protons + pz ) - log( E_protons - pz ) ) - y_beam;
    auto mass_pions = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    auto E_pions = sqrt( p*p + mass_pions*mass_pions );
    float y_pions = 0.5 * ( log( E_pions + pz ) - log( E_pions - pz ) );
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
      efficiency = 1.0f/efficiency;
    else
      efficiency = 0.0;
    auto out_particle = out_tracks_->NewChannel();
    out_particle.CopyContents(in_track);
    out_particle[out_ycm_var_].SetVal(y_cm);
    out_particle[out_abs_ycm_var_].SetVal(fabsf(y_cm));
    out_particle[out_efficiency_var_].SetVal(efficiency);
    out_particle[out_protons_rapidity_var_].SetVal(y_protons);
    out_particle[out_pions_rapidity_var_].SetVal(y_pions);
    in_track.DataT<Particle>()->SetMass(mass);
  }
  if( is_mc_ ){
    out_sim_particles_->ClearChannels();
    for ( auto in_particle : in_sim_particles_->Loop() ) {
      auto out_particle = out_sim_particles_->NewChannel();
      out_particle.CopyContents(in_particle);
      auto pid = in_particle.DataT<Particle>()->GetPid();
      auto mass = in_particle.DataT<Particle>()->GetMass();
      auto mom4 = in_particle.DataT<Particle>()->Get4MomentumByMass(mass);
      auto pT = mom4.Pt();
      auto p = mom4.P();
      auto pz = mom4.Pz();
      float y = mom4.Rapidity();
      float y_cm = y- y_beam;
      auto mass_protons = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
      auto E_protons = sqrt( p*p + mass_protons*mass_protons );
      float y_protons = 0.5 * ( log( E_protons + pz ) - log( E_protons - pz ) ) - y_beam;
      auto mass_pions = TDatabasePDG::Instance()->GetParticle(211)->Mass();
      auto E_pions = sqrt( p*p + mass_pions*mass_pions );
      float y_pions = 0.5 * ( log( E_pions + pz ) - log( E_pions - pz ) );
      int charge=0;
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge() / 3;
      else{
        charge = (int)(pid / 1E+4) % (int)1e+3;
      }
      out_particle[out_sim_ycm_var_].SetVal(y_cm);
      out_particle[out_sim_abs_ycm_var_].SetVal(fabsf(y_cm));
      out_particle[out_sim_protons_rapidity_var_].SetVal(y_protons);
      out_particle[out_sim_pions_rapidity_var_].SetVal(y_pions);
      out_particle[out_sim_is_charged_].SetVal( charge != 0 );
      in_particle.DataT<Particle>()->SetMass(mass);
    }
  }
}

void TracksProcessor::ReadEfficiencyHistos(){
  auto protons_file_name = config_directory_;
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
//bool TracksProcessor::UseATI2() const { return false; }
