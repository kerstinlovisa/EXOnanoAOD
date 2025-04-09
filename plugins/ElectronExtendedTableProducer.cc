// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class ElectronExtendedTableProducer : public edm::global::EDProducer<> {
  public:
    explicit ElectronExtendedTableProducer(const edm::ParameterSet &iConfig) :
      name_(iConfig.getParameter<std::string>("name")),
      electronTag_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons")))
    {
      produces<nanoaod::FlatTable>();
    }

    ~ElectronExtendedTableProducer() override {};

    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("electrons")->setComment("input electron collection");
      desc.add<std::string>("name")->setComment("name of the electron nanoaod::FlatTable we are extending");
      descriptions.add("electronTable", desc);
    }

  private:
    void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

    std::string name_;
    edm::EDGetTokenT<std::vector<pat::Electron>> electronTag_;

};

void ElectronExtendedTableProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const 
{
  edm::Handle<std::vector<pat::Electron>> electrons;
  iEvent.getByToken(electronTag_, electrons);

  unsigned int nElectrons = electrons->size();
  std::vector<float> idx;
  for (unsigned int i = 0; i < nElectrons; i++) {
    idx.push_back(i);
  }

  auto tab  = std::make_unique<nanoaod::FlatTable>(nElectrons, name_, false, true);
  tab->addColumn<float>("idx", idx, "LLPnanoAOD electron index");
  iEvent.put(std::move(tab));
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(ElectronExtendedTableProducer);
