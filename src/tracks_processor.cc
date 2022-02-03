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
      ("q-vector-file", value(&str_qvector_file_name_), "Path to the file containing Q-vectors")
      ("q-vector-name", value(&str_qvector_name_), "Path to the file containing Q-vectors")
      ("protons-efficiency", value(&str_protons_efficiency_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to file with proton's efficiency")
      ("centrality-file", value(&str_centrality_file_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to file with proton's efficiency")
      ("efficiency-delta-phi", value(&str_efficiency_delta_phi_)->default_value("../../efficiency_files/protons_efficiency.root"), "Path to file with proton's efficiency")
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
  out_fabs_pid_var_ = out_tracks_->NewVariable( "abs_pid", INTEGER );
  out_c_eff_var_ = out_tracks_->NewVariable( "c_eff", FLOAT );

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
  if( tree_qvector_ )
    tree_qvector_->GetEntry(qvector_event_);
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if( h1_centrality_bins_ ){
    auto tof_rpc_hits = (*event_header_)[in_multiplicity_var_].GetInt();
    auto mult_bin = h1_centrality_bins_->GetXaxis()->FindBin(tof_rpc_hits);
    auto c_class = h1_centrality_bins_->GetBinContent(mult_bin);
    centrality = 5.0f * c_class - 2.5f;
    (*out_event_header_)[out_centrality_var_].SetVal(centrality);
  }
  out_tracks_->ClearChannels();
  this->LoopRecTracks();
  if( is_mc_ ){
    this->LoopSimParticles();
  }
  qvector_event_++;
}

void TracksProcessor::ReadEfficiencyHistos(){
  file_efficiency_protons_ = TFile::Open(str_protons_efficiency_.c_str(), "read");
  file_efficiency_pi_plus_ = TFile::Open(str_pi_plus_efficiency_.c_str(), "read");
  file_efficiency_pi_minus_ = TFile::Open(str_pi_minus_efficiency_.c_str(), "read");
  file_centrality_ = TFile::Open(str_centrality_file_.c_str(), "read");
  if(file_centrality_)
    file_centrality_->GetObject("Centrality/TOFRPC_5pc_fixedCuts",h1_centrality_bins_);
  // Initializing efficiency for occupancy correction
  file_efficiency_delta_phi_ = TFile::Open( str_efficiency_delta_phi_.c_str(), "read" );
  if( file_efficiency_delta_phi_ ){
    file_efficiency_delta_phi_->GetObject("h3_rec_delta_phi_theta_centrality_all",
        h3_npart_delta_phi_theta_centrality_);
  }

  file_qvector_ = TFile::Open(str_qvector_file_name_.c_str());
  if( !file_qvector_ && !str_qvector_file_name_.empty() )
    std::cout << "File "+str_qvector_file_name_+" does not exist" << std::endl;
  if( file_qvector_ ) {
    file_qvector_->GetObject("tree", tree_qvector_);
    if (!tree_qvector_)
      throw std::runtime_error("There is no tree in file " +
                               str_qvector_file_name_);
    tree_qvector_->SetBranchAddress(str_qvector_name_.c_str(), &dc_qvector_);
    if (!dc_qvector_)
      throw std::runtime_error("There is no branch " + str_qvector_name_ +
                               " in tree");
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
  out_file_->cd();

}
void TracksProcessor::LoopRecTracks() {
  using AnalysisTree::Particle;
  auto centrality = (*event_header_)[GetVar( "event_header/selected_tof_rpc_hits_centrality" )].GetVal();
  auto centrality_class = (size_t) ( (centrality-2.5)/5.0 );
  if (h1_centrality_bins_) {
    auto tof_rpc_hits =
        (*event_header_)[GetVar("event_header/selected_tof_rpc_hits")]
            .GetInt();
    auto mult_bin = h1_centrality_bins_->GetXaxis()->FindBin(tof_rpc_hits);
    auto c_class = h1_centrality_bins_->GetBinContent(mult_bin);
    centrality = 5.0 * c_class - 2.5;
  }
  auto psi_ep = 0.0;
  if (dc_qvector_) {
    auto qvec = dc_qvector_->At(0);
    auto mag = qvec.mag(1);
    if (fabs(mag) < std::numeric_limits<double>::min())
      return;
    auto x_qvec = qvec.x(1);
    auto y_qvec = qvec.y(1);
    psi_ep = atan2(y_qvec, x_qvec);
  }
  float y_beam = data_header_->GetBeamRapidity();
  double c_eff_protons = 0.0045;
  double c_eff_pi_neg = 0.05;
  double c_eff_pi_pos = 0.005;
  while( c_eff_pi_pos < 0.101 ) {
    for (auto in_track : in_tracks_->Loop()) {
      auto pid = in_track.DataT<Particle>()->GetPid();
      auto mass = in_track.DataT<Particle>()->GetMass();
      if (pid != 0) {
        if (TDatabasePDG::Instance()->GetParticle(pid))
          mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
      }
      auto mom4 = in_track.DataT<Particle>()->Get4MomentumByMass(mass);
      auto pT = mom4.Pt();
      auto p = mom4.P();
      auto pz = mom4.Pz();
      auto theta = mom4.Theta();
      float y = mom4.Rapidity();
      float y_cm = y - y_beam;
      TH2F *efficiency_histogram{nullptr};
      float efficiency{1.0};
      if (pid == 211) {
        efficiency_histogram = centrality_class < efficiency_pi_plus_.size()
                                   ? efficiency_pi_plus_.at(centrality_class)
                                   : nullptr;
      }
      if (pid == -211) {
        efficiency_histogram = centrality_class < efficiency_pi_minus_.size()
                                   ? efficiency_pi_minus_.at(centrality_class)
                                   : nullptr;
      }
      if (pid == 2212) {
        efficiency_histogram = centrality_class < efficiency_protons_.size()
                                   ? efficiency_protons_.at(centrality_class)
                                   : nullptr;
      }
      if (efficiency_histogram) {
        auto bin_y = efficiency_histogram->GetXaxis()->FindBin(y_cm);
        auto bin_pT = efficiency_histogram->GetYaxis()->FindBin(pT);
        efficiency = efficiency_histogram->GetBinContent(bin_y, bin_pT);
      }
      auto occupancy_weight = 0.0;
      auto c_eff=0.0;
      if (h3_npart_delta_phi_theta_centrality_) {
        if (pid == 2212) {
          auto delta_phi = AngleDifference(mom4.Phi(), psi_ep);
          auto dphi_bin =
              h3_npart_delta_phi_theta_centrality_->GetXaxis()->FindBin(
                  delta_phi);
          auto theta_bin =
              h3_npart_delta_phi_theta_centrality_->GetYaxis()->FindBin(
                  mom4.Theta());
          auto c_bin =
              h3_npart_delta_phi_theta_centrality_->GetZaxis()->FindBin(
                  centrality);
          auto n_tracks = h3_npart_delta_phi_theta_centrality_->GetBinContent(
              dphi_bin, theta_bin, c_bin);
          auto eff = 0.98 - c_eff_protons * pow(n_tracks, 2);
          if (eff > 0.1)
            occupancy_weight = 1.0 / eff;
          c_eff = c_eff_protons;
        }
        if (pid == -211) {
          auto delta_phi = AngleDifference(mom4.Phi(), psi_ep);
          auto dphi_bin =
              h3_npart_delta_phi_theta_centrality_->GetXaxis()->FindBin(
                  delta_phi);
          auto theta_bin =
              h3_npart_delta_phi_theta_centrality_->GetYaxis()->FindBin(
                  mom4.Theta());
          auto c_bin =
              h3_npart_delta_phi_theta_centrality_->GetZaxis()->FindBin(
                  centrality);
          auto n_tracks = h3_npart_delta_phi_theta_centrality_->GetBinContent(
              dphi_bin, theta_bin, c_bin);
          auto eff = 0.98 - c_eff_pi_neg * pow(n_tracks, 1);
          if (eff > 0.1)
            occupancy_weight = 1.0 / eff;
          c_eff = c_eff_pi_neg;
        }
        if (pid == 211) {
          auto delta_phi = AngleDifference(mom4.Phi(), psi_ep);
          auto dphi_bin =
              h3_npart_delta_phi_theta_centrality_->GetXaxis()->FindBin(
                  delta_phi);
          auto theta_bin =
              h3_npart_delta_phi_theta_centrality_->GetYaxis()->FindBin(
                  mom4.Theta());
          auto c_bin =
              h3_npart_delta_phi_theta_centrality_->GetZaxis()->FindBin(
                  centrality);
          auto n_tracks = h3_npart_delta_phi_theta_centrality_->GetBinContent(
              dphi_bin, theta_bin, c_bin);
          auto eff = 0.98 - c_eff_pi_pos * pow(n_tracks, 1);
          if (eff > 0.1)
            occupancy_weight = 1.0 / eff;
          c_eff = c_eff_pi_pos;
        }
      }
      auto out_particle = out_tracks_->NewChannel();
      out_particle.CopyContents(in_track);
      out_particle[out_ycm_var_] = y_cm;
      out_particle[out_theta_var_] = (float)theta;
      out_particle[out_efficiency_var_] =
          efficiency > 0.3f ? 1.0f / efficiency : 0.0f;
      out_particle[out_occ_weight_var_] = occupancy_weight;
      out_particle[out_fabs_pid_var_] = (int)abs(pid);
      out_particle[out_c_eff_var_] = (float)c_eff;
      out_particle.DataT<Particle>()->SetMass(mass);
    }
    c_eff_pi_pos+=0.01;
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