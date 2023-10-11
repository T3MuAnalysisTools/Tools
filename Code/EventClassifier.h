#ifndef EventClassifier_h
#define EventClassifier_h
#include "Ntuple_Controller.h"
#include "PDG_Var.h"
#include "PDGInfo.h"


class EventClassifier {

 public:
  EventClassifier(Ntuple_Controller *ntp);
  ~EventClassifier();

  enum DecayCategory{AmBmm=1, mAmm, mmAm, AmBmCm, mAmBm, mmm, XmYmm, XmYmZm};


  bool isDecayInFlight(unsigned int i);
  bool isKpiFake(unsigned int i);

  std::vector<unsigned int> findDuplicates(std::vector<unsigned int> vec);
  std::vector<unsigned int> findSameAncestors(std::vector<unsigned int> vec);
  std::vector<int>  triplet_has_mother(std::vector<unsigned int> vec);

  int EventType(std::vector<unsigned int> vec);
  int ClassifyTypeI(std::vector<unsigned int> vec);
  int ClassifyTypeII(std::vector<unsigned int> vec);
  int TypeI1_I2_pair_parent(std::vector<unsigned int> vec);
  std::vector<unsigned int> AllParentsOfParticle(int index,std::vector<unsigned int> out = {});
  std::vector<unsigned int> Intersection(std::vector<unsigned int> v1, std::vector<unsigned int> v2);


 private:
  Ntuple_Controller *Ntp;

};
#endif
