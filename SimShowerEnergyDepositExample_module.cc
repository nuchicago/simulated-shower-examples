//////////////////////////////////////////////////////////////
// Name:      SimShowerEnergyDepositExample
// Date:      15 October 2018
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef SimShowerEnergyDepositExample_Module
#define SimShowerEnergyDepositExample_Module

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h" 
//#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// NuTools includes
#include "nusimdata/SimulationBase/MCParticle.h"
//#include "nusimdata/SimulationBase/MCTrajectory.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
//#include "larsim/MCCheater/BackTrackerService.h"
//#include "larsim/Simulation/LArG4Parameters.h"

// ROOT includes
#include "TTree.h"

// C++ includes
#include <map>
#include <set>
#include <string>
#include <vector>

// type definition of particle maps
typedef std::map< int, art::Ptr< simb::MCParticle > > ParticleMap;
typedef std::map< int, ParticleMap > ShowerParticleMaps;

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class SimShowerEnergyDepositExample : public art::EDAnalyzer
{

 public:

  // standard constructor and destructor for an art module
  explicit SimShowerEnergyDepositExample(fhicl::ParameterSet const& pset);
  virtual ~SimShowerEnergyDepositExample();

  // this method is called once, at the start of the job
  virtual void beginJob() override;

  // this method is called once, at the start of each run
  virtual void beginRun(const art::Run & run) override;

  // this method is called once, at the start of each subrun
  virtual void beginSubRun(const art::SubRun & subrun) override;

  // this method is called once, at the end of the job
  virtual void endJob() override;

  // this method is called once, at the end of each run
  virtual void endRun(const art::Run & run) override;

  // this method is called once, at the end of each subrun
  virtual void endSubRun(const art::SubRun & subrun) override;

  // this method reads in any parameters from the .fcl files
  void reconfigure(fhicl::ParameterSet const& pset) override;

  // the analyze routine, called once per event
  void analyze(art::Event const & event) override;

 private:

  // parameters read from FHiCL (.fcl) file
  std::string simulation_producer_label_;
  std::string mcshower_producer_label_;

  // pointer to geometry provider
  geo::GeometryCore const* geometry_;

  // pointer to detector properties
  //detinfo::DetectorProperties const* detector_;

  // reset once per event
  void reset_();

  // get particle map
  ParticleMap particle_map_(
      std::vector< art::Ptr< simb::MCParticle > > const& particle_vector);

  // get shower particle maps
  ShowerParticleMaps shower_particle_maps_(
      ParticleMap                              const& particle_map,
      std::vector< art::Ptr< sim::MCShower > > const& mcshower_vector);

  // fill shower particle map
  void fill_shower_particle_map_(
      int         const  track_id,
      ParticleMap const& particle_map,
      ParticleMap      & shower_particle_map);

  // get energy deposited
  double sim_energy_deposited_(
      std::set< int >                            const& track_ids,
      std::vector< art::Ptr< sim::SimChannel > > const& simchannel_vector);

  // pointers to TTree object
  TTree * ttree_;

  // variables that will go into the TTree objects
  int event_;     // number of the event being processed
  int run_;       // number of the run being processed
  int subrun_;    // number of the sub-run being processed
};

//-----------------------------------------------------------------------
// constructor
SimShowerEnergyDepositExample::SimShowerEnergyDepositExample(fhicl::ParameterSet const& pset)
  :
  EDAnalyzer(pset)
{
  // reconfigure parameters
  this->reconfigure(pset);

  // get a pointer to the geometry service provider
  geometry_ = &*(art::ServiceHandle<geo::Geometry>());

  // get a pointer to the detector properties provider
  //detector_ = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

//-----------------------------------------------------------------------
// destructor
SimShowerEnergyDepositExample::~SimShowerEnergyDepositExample() {}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::beginJob()
{
  // access art's TFileService
  art::ServiceHandle< art::TFileService > tfs;

  ttree_ = tfs->make<TTree>("simshowerenergydepositexample", "simshowerenergydepositexample");

  ttree_->Branch("event",  &event_,  "event/I");
  ttree_->Branch("run",    &run_,    "run/I");
  ttree_->Branch("subrun", &subrun_, "subrun/I");
}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::beginRun(const art::Run & /*run*/)
{}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::beginSubRun(const art::SubRun & /*subrun*/)
{}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::endJob()
{}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::endRun(const art::Run & /*run*/)
{}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::endSubRun(const art::SubRun & /*subrun*/)
{}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::reconfigure(fhicl::ParameterSet const& pset)
{
  // read parameters from the .fcl file
  simulation_producer_label_ = pset.get< std::string >("SimulationLabel", "largeant");
  mcshower_producer_label_   = pset.get< std::string >("MCShowerLabel",   "mcreco");
}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::analyze(art::Event const & event)
{
  //-------------------------------------------------------------------
  // reset once per event
  //-------------------------------------------------------------------

  //this->reset_();

  //-------------------------------------------------------------------
  // get event, run, and subrun numbers
  //-------------------------------------------------------------------

  event_  = event.id().event();
  run_    = event.run();
  subrun_ = event.subRun();

  //-------------------------------------------------------------------
  // get data products
  //-------------------------------------------------------------------

  // get art::ValidHandle< std::vector< sim::SimChannel > >
  auto simchannel_handle = event.getValidHandle< std::vector< sim::SimChannel > >
      (simulation_producer_label_);

  // fill vector of SimChannel objects
  std::vector< art::Ptr< sim::SimChannel > > simchannel_vector;
  art::fill_ptr_vector(simchannel_vector, simchannel_handle);

  // get art::ValidHandle< std::vector< simb::MCParticle > >
  auto particle_handle = event.getValidHandle< std::vector< simb::MCParticle > >
      (simulation_producer_label_);

  // fill vector of MCParticle objects
  std::vector< art::Ptr< simb::MCParticle > > particle_vector;
  art::fill_ptr_vector(particle_vector, particle_handle);

  // get art::ValidHandle< std::vector< sim::MCShower > >
  auto mcshower_handle = event.getValidHandle< std::vector< sim::MCShower > >
      (mcshower_producer_label_);

  // fill vector of MCShower objects
  std::vector< art::Ptr< sim::MCShower > > mcshower_vector;
  art::fill_ptr_vector(mcshower_vector, mcshower_handle);

  //-------------------------------------------------------------------
  // fill particle maps
  //-------------------------------------------------------------------

  // ParticleMap
  auto const particle_map = this->particle_map_(particle_vector);

  // ShowerParticleMaps
  auto const shower_particle_maps = this->shower_particle_maps_(
      particle_map, mcshower_vector);

  //-------------------------------------------------------------------
  // get energy deposited by shower particles
  //-------------------------------------------------------------------

  // loop over map of particle maps
  for (auto const& shower_particle_maps_iter : shower_particle_maps)
  {
    // get ancestor shower particle track ID
    int const ancestor_track_id = shower_particle_maps_iter.first;

    // get ancestor shower particle
    // art::Ptr< simb::MCParticle > const ancestor_particle
    //     = particle_map.at(ancestor_track_id);

    // set of track IDs of descendent particles
    std::set< int > descendent_track_ids;

    // loop over particle maps for descendent particles
    for (auto const& particle_map_iter : shower_particle_maps_iter.second)
    {
      // get track ID of descendent particle
      int const descendent_track_id = particle_map_iter.first;

      // insert into set of descendent track IDs
      descendent_track_ids.insert(descendent_track_id);
    } // end loop over particle maps

    //-------------------------------------------------------
    // this is where you can do something with the total
    // visible energy deposited on the TPC wires by the
    // electromagnetic shower
    //-------------------------------------------------------

    // get total visible energy deposited by all descendent particles
    double const energy_deposited
        = this->sim_energy_deposited_(descendent_track_ids, simchannel_vector);

  } // end loop over map of particle maps

  //-------------------------------------------------------------------
  // fill TTree object
  //-------------------------------------------------------------------

  ttree_->Fill();
}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::reset_()
{}

//-----------------------------------------------------------------------
ParticleMap SimShowerEnergyDepositExample::particle_map_(
    std::vector< art::Ptr< simb::MCParticle > > const& particle_vector)
{
  // initialize particle map
  ParticleMap particle_map;

  // loop over simb::MCParticle objects and fill particle map
  for (auto const& particle : particle_vector)
  {
    // add the address of the MCParticle to the map, with the
    // track ID as the key
    particle_map[particle->TrackId()] = particle;
  } // end loop over simb::MCParticle objects

  return particle_map;
}

//-----------------------------------------------------------------------
ShowerParticleMaps SimShowerEnergyDepositExample::shower_particle_maps_(
    ParticleMap                              const& particle_map,
    std::vector< art::Ptr< sim::MCShower > > const& mcshower_vector)
{
  // initialize shower particle maps
  ShowerParticleMaps shower_particle_maps;

  // loop over MCShower objects
  for (auto const& mcshower : mcshower_vector)
  {
    // get track ID and particle
    int const track_id = mcshower->TrackID();
    auto const& particle = particle_map.at(track_id);

    // fill shower particle maps
    this->fill_shower_particle_map_(
        particle->TrackId(), particle_map,
        shower_particle_maps[particle->TrackId()]);
  } // end loop over MCShower objects

  return shower_particle_maps;
}

//-----------------------------------------------------------------------
void SimShowerEnergyDepositExample::fill_shower_particle_map_(
    int         const  track_id,
    ParticleMap const& particle_map,
    ParticleMap      & shower_particle_map)
{
  // create iterator for particle map
  ParticleMap::const_iterator particle_iter = particle_map.find(track_id);

  // return if the end of the particle map is reached
  if (particle_iter == particle_map.end()) return;

  // pointer to MCParticle object
  art::Ptr< simb::MCParticle > particle = (*particle_iter).second;

  // get number of daughters of particle
  int number_daughters = particle->NumberDaughters();

  // add pointer of MCParticle object to the shower particle map
  shower_particle_map[track_id] = particle;

  // loop over daughter particles and recur
  for (int idx = 0; idx < number_daughters; ++idx)
  {
    // get daughter track ID
    int daughter_track_id = particle->Daughter(idx);

    // recur
    this->fill_shower_particle_map_(
        daughter_track_id, particle_map, shower_particle_map);
  } // end loop over daughter particles
}

//-----------------------------------------------------------------------
double SimShowerEnergyDepositExample::sim_energy_deposited_(
    std::set< int >                            const& track_ids,
    std::vector< art::Ptr< sim::SimChannel > > const& simchannel_vector)
{
  double energy_deposited = 0;

  // look at all the energy deposited on the TPC wires
  // sim::SimChannel
  for (auto const& channel : simchannel_vector)
  {
    // get the numeric ID associated with this channel
    // raw::ChannelID_t
    auto channel_number = channel->Channel();

    // only look at the energy on the collection plane
    if (geometry_->SignalType(channel_number) != geo::kCollection)
        continue;

    // each channel has a map inside it that connects
    // a time slice to energy deposits in the detector
    // std::map< unsigned short, std::vector< sim::IDE > >
    auto const& time_slices = channel->TDCIDEMap();

    // for every time slice in this channel:
    for (auto const& time_slice : time_slices)
    {
      // get energy deposits from time slice
      // std::vector< sim::IDE >
      auto const& energy_deposits = time_slice.second;

       // loop through energy deposits
       // sim::IDE
       for (auto const& energy_deposit : energy_deposits)
       {
         if (track_ids.find(energy_deposit.trackID) != track_ids.end() and
             energy_deposit.trackID                 != sim::NoParticleId)
         {
           // add to energy deposited on TPC wires
           energy_deposited += energy_deposit.energy;
         } // if energy deposited by particle
       } // for each energy deposit
    } // for each time slice
  } // for each SimChannel

  return energy_deposited;
}


DEFINE_ART_MODULE(SimShowerEnergyDepositExample)

#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop

#endif // SimShowerEnergyDepositExample_module
