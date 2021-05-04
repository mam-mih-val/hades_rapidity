//
// Created by mikhail on 11/30/20.
//

#include "tracks_processor.h"

#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TMath.h"

#include <AnalysisTree/DataHeader.hpp>

TASK_IMPL(TracksProcessor)

boost::program_options::options_description TracksProcessor::GetBoostOptions() {
  using namespace boost::program_options;

  options_description desc(GetName() + " options");
  desc.add_options()
      ("protons-efficiency", value(&protons_efficiency_file_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to file with proton's efficiency")
      ("pi-plus-efficiency", value(&pi_plus_efficiency_file_)->default_value("../../efficiency_files/pi_plus_efficiency.root"), "Path to file with pi-plus' efficiency")
      ("pi-minus-efficiency", value(&pi_minus_efficiency_file_)->default_value("../../efficiency_files/pi_minus_efficiency.root"), "Path to file with pi-minus' efficiency")
      ("deutrons-efficiency", value(&deutrons_efficiency_file_)->default_value("../../efficiency_files/pi_minus_efficiency.root"), "Path to file with pi-minus' efficiency");
  return desc;
}

void TracksProcessor::PreInit() {}

void TracksProcessor::UserInit(std::map<std::string, void *> &Map) {
  in_tracks_ = GetInBranch("mdc_vtx_tracks");
  event_header_ = GetInBranch("event_header");
  centrality_var_ = GetVar( "event_header/selected_tof_rpc_hits_centrality" );
  geant_pid_var_ = GetVar( "mdc_vtx_tracks/geant_pid" );
  charge_var_ = GetVar( "mdc_vtx_tracks/charge" );
  out_tracks_ = NewBranch( "mdc_vtx_tracks_extra", PARTICLES );
  out_tracks_->CloneVariables(in_tracks_->GetConfig());
  out_is_positive_ = out_tracks_->NewVariable( "is_positive", BOOLEAN );
  out_is_in_protons_acceptance_ = out_tracks_->NewVariable( "is_in_protons_acceptance", BOOLEAN );
  out_ycm_var_ = out_tracks_->NewVariable( "ycm", FLOAT );
  out_abs_ycm_var_ = out_tracks_->NewVariable( "abs_ycm", FLOAT );
  out_efficiency_var_ = out_tracks_->NewVariable( "efficiency", FLOAT );
  out_protons_rapidity_var_ = out_tracks_->NewVariable( "protons_rapidity", FLOAT );
  out_protons_pT_var_ = out_tracks_->NewVariable( "protons_pT", FLOAT );
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
    out_sim_theta_var_ = out_sim_particles_->NewVariable( "theta", FLOAT );
    out_sim_ycm_var_ = out_sim_particles_->NewVariable( "ycm", FLOAT );
    out_sim_abs_ycm_var_ = out_sim_particles_->NewVariable( "abs_ycm", FLOAT );
    out_sim_protons_rapidity_var_ = out_sim_particles_->NewVariable( "protons_rapidity", FLOAT );
    out_sim_protons_pT_var_ = out_sim_particles_->NewVariable( "protons_pT", FLOAT );
    out_sim_pions_rapidity_var_ = out_sim_particles_->NewVariable( "pions_rapidity", FLOAT );
    out_sim_is_charged_ = out_sim_particles_->NewVariable( "is_charged", BOOLEAN );
    out_sim_is_positive_ = out_sim_particles_->NewVariable( "is_positive", BOOLEAN );
    out_sim_is_in_acceptance_ = out_sim_particles_->NewVariable( "is_in_acceptance", BOOLEAN );
    out_sim_is_in_protons_acceptance_ = out_sim_particles_->NewVariable( "is_in_protons_acceptance", BOOLEAN );
    out_sim_is_in_high_efficiency_ = out_sim_particles_->NewVariable( "is_in_high_efficiency", BOOLEAN );
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
    auto geant_pid = in_track[geant_pid_var_].GetInt();
    int charge=in_track[charge_var_].GetInt();
    auto mass = in_track.DataT<Particle>()->GetMass();
    if( pid != 0 ) {
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
    }
    auto mom4 = in_track.DataT<Particle>()->Get4MomentumByMass(mass);
    auto mass_protons = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    auto pT = mom4.Pt();
    auto p = mom4.P();
    auto pz = mom4.Pz();
    auto theta = mom4.Theta();
    float y = mom4.Rapidity();
    float y_cm = y- y_beam;
    float E_protons = sqrtf(p * p + mass_protons*mass_protons );
    float y_protons = 0.5f * ( logf( E_protons + pz) - logf( E_protons - pz) ) - y_beam;
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
      default:
        switch (geant_pid) {
          case 45:
            efficiency_histogram = efficiency_deutrons_.at(centrality_class);
            break;
        }
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
    double protons_efficiency=0.0;
    try{
      efficiency_histogram = efficiency_protons_.at(centrality_class);
      if (efficiency_histogram) {
        auto bin_y = efficiency_histogram->GetXaxis()->FindBin(y_protons);
        auto bin_pT = efficiency_histogram->GetYaxis()->FindBin(pT);
        protons_efficiency = efficiency_histogram->GetBinContent(bin_y, bin_pT);
      }
    } catch (std::exception&) {}
    auto out_particle = out_tracks_->NewChannel();
    out_particle.CopyContents(in_track);
    out_particle[out_ycm_var_].SetVal(y_cm);
    out_particle[out_abs_ycm_var_].SetVal(fabsf(y_cm));
    out_particle[out_efficiency_var_].SetVal(efficiency);
    out_particle[out_protons_rapidity_var_].SetVal(y_protons);
    out_particle[out_is_positive_].SetVal( charge > 0 );
    out_particle[out_is_in_protons_acceptance_].SetVal( protons_efficiency > 0.1 );
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
      auto mass_protons = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
      auto pT = mom4.Pt();
      auto theta = mom4.Theta();
      auto p = mom4.P();
      auto pz = mom4.Pz();
      float y = mom4.Rapidity();
      float y_cm = y- y_beam;
      auto E_protons = sqrt(p * p + mass_protons*mass_protons );
      float y_protons = 0.5 * ( log( E_protons + pz) - log( E_protons - pz) ) - y_beam;
      int charge=0;
      if( TDatabasePDG::Instance()->GetParticle(pid) )
        charge = TDatabasePDG::Instance()->GetParticle(pid)->Charge() / 3;
      else{
        charge = (int)(pid / 1E+4) % (int)1e+3;
      }
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
        }
        if (efficiency_histogram) {
          auto bin_y = efficiency_histogram->GetXaxis()->FindBin(y_cm);
          auto bin_pT = efficiency_histogram->GetYaxis()->FindBin(pT);
          efficiency = efficiency_histogram->GetBinContent(bin_y, bin_pT);
        }
      } catch (std::exception&) {}
      double protons_efficiency=0.0;
      try{
        efficiency_histogram = efficiency_protons_.at(centrality_class);
        if (efficiency_histogram) {
          auto bin_y = efficiency_histogram->GetXaxis()->FindBin(y_protons);
          auto bin_pT = efficiency_histogram->GetYaxis()->FindBin(pT);
          protons_efficiency = efficiency_histogram->GetBinContent(bin_y, bin_pT);
        }
      } catch (std::exception&) {}
      double p_hi = (4.95-2.5*theta)/1.2;
      double p_lo = 0.4;
      auto theta_lo = 0.3;
      auto theta_hi = 1.5;
      bool is_in_acceptance{true};
      if( p < p_lo )
        is_in_acceptance = false;
      if( p > p_hi )
        is_in_acceptance = false;
      if( theta > theta_hi )
        is_in_acceptance = false;
      if( theta < theta_lo )
        is_in_acceptance = false;
      out_particle[out_sim_theta_var_].SetVal((float)mom4.Theta());
      out_particle[out_sim_ycm_var_].SetVal(y_cm);
      out_particle[out_sim_abs_ycm_var_].SetVal(fabsf(y_cm));
      out_particle[out_sim_protons_rapidity_var_].SetVal(y_protons);
      out_particle[out_sim_is_charged_].SetVal( charge != 0 );
      out_particle[out_sim_is_positive_].SetVal( charge > 0 );
      out_particle[out_sim_is_in_acceptance_].SetVal( is_in_acceptance );
      out_particle[out_sim_is_in_high_efficiency_].SetVal( efficiency > 0.1 );
      out_particle[out_sim_is_in_protons_acceptance_].SetVal( protons_efficiency > 0.1 );
      in_particle.DataT<Particle>()->SetMass(mass);
    }
  }
}

void TracksProcessor::ReadEfficiencyHistos(){
  file_efficiency_protons_ = TFile::Open(protons_efficiency_file_.c_str(), "read");
  file_efficiency_pi_plus_ = TFile::Open(pi_plus_efficiency_file_.c_str(), "read");
  file_efficiency_pi_minus_ = TFile::Open(pi_minus_efficiency_file_.c_str(), "read");
  file_efficiency_deutrons_ = TFile::Open(deutrons_efficiency_file_.c_str(), "read");
  int p=2;
  while(p<40){
    efficiency_protons_.emplace_back();
    efficiency_pi_plus_.emplace_back();
    efficiency_pi_minus_.emplace_back();
    efficiency_deutrons_.emplace_back();
    std::string name = "efficiency_"+std::to_string(p);
    if( file_efficiency_protons_ )
      file_efficiency_protons_->GetObject(name.c_str(), efficiency_protons_.back());
    if( file_efficiency_pi_plus_ )
      file_efficiency_pi_plus_->GetObject(name.c_str(), efficiency_pi_plus_.back());
    if(file_efficiency_pi_minus_)
      file_efficiency_pi_minus_->GetObject(name.c_str(), efficiency_pi_minus_.back());
    if(file_efficiency_deutrons_)
      file_efficiency_deutrons_->GetObject(name.c_str(), efficiency_deutrons_.back());
    p+=5;
  }
}
//bool TracksProcessor::UseATI2() const { return false; }
