#include "EventClassifier.h"


EventClassifier::EventClassifier(Ntuple_Controller *ntp)
{
  Ntp=ntp;
}

EventClassifier::~EventClassifier()
{
}


bool
EventClassifier::isKpiFake(unsigned int i) //   so far fake only means K/pi
{

  if(isDecayInFlight(i)) return false;

  if(abs(Ntp->MCParticle_pdgid(i))  ==  PDGInfo::mu_minus) return false;

  if( abs(Ntp->MCParticle_pdgid(i)) ==  PDGInfo::pi_plus   ||
      abs(Ntp->MCParticle_pdgid(i)) ==  PDGInfo::K_plus ) return true;
  
  
  return false;
}      



bool
EventClassifier::isDecayInFlight(unsigned int i)
{

  if(Ntp->MCParticle_hasMother(i))
    if( abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(i))) ==  PDGInfo::pi_plus   ||
        abs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(i))) ==  PDGInfo::K_plus )
      {

	std::vector<unsigned int> decay_chain;
        for(auto j : Ntp->MCParticle_childidx(Ntp->MCParticle_midx(i)))
          decay_chain.push_back(abs(Ntp->MCParticle_pdgid(j)));
        if(decay_chain.size() == 2)
          if(std::find(decay_chain.begin(), decay_chain.end(), PDGInfo::mu_minus) != decay_chain.end()&&
             std::find(decay_chain.begin(), decay_chain.end(), PDGInfo::nu_mu)    != decay_chain.end() )
            return true;
      }

  if( abs(Ntp->MCParticle_pdgid(i)) ==  PDGInfo::pi_plus   ||
      abs(Ntp->MCParticle_pdgid(i)) ==  PDGInfo::K_plus )
    {
      std::vector<unsigned int> decay_chain;
      for(auto j : Ntp->MCParticle_childidx(i))
        decay_chain.push_back(abs(Ntp->MCParticle_pdgid(j)));
      if(decay_chain.size() == 2)
        if(std::find(decay_chain.begin(), decay_chain.end(), PDGInfo::mu_minus) != decay_chain.end()&&
           std::find(decay_chain.begin(), decay_chain.end(), PDGInfo::nu_mu)    != decay_chain.end() )
          return true;

    }

  return false;

}




std::vector<int> 
EventClassifier::triplet_has_mother(std::vector<unsigned int> vec)
{
  std::vector<int> out;
  for(auto i : vec)
    out.push_back(Ntp->MCParticle_hasMother(i));
  return out;
}



std::vector<unsigned int> 
EventClassifier::findDuplicates(std::vector<unsigned int> vec)
{
  std::vector<unsigned int> duplicates;
  if(vec.size()!=0)
    for (unsigned int i=0; i < vec.size() - 1 ; ++i)
      for(unsigned int j=i+1; j < vec.size(); ++j)
	if(vec.at(i) == vec.at(j)) 
	  duplicates.push_back(vec.at(i));

  return duplicates;
}



std::vector<unsigned int> 
EventClassifier::AllParentsOfParticle(int index, std::vector<unsigned int> out)
{

  if(Ntp->MCParticle_hasMother(index))
    if(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(index)) != 2212)
      {
	
	out.push_back(Ntp->MCParticle_midx(index));
	return AllParentsOfParticle(Ntp->MCParticle_midx(index), out);

      }
  return out;
}


std::vector<unsigned int> EventClassifier::Intersection(std::vector<unsigned int> v1, std::vector<unsigned int> v2)
{
  std::vector<unsigned int> out(v1.size() + v2.size());
  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());
  std::vector<unsigned int>::iterator it;
  it=std::set_intersection (v1.begin(), v1.end(), v2.begin(), v2.end(), out.begin());
  out.resize(it - out.begin());

  return out;

}

int 
EventClassifier::EventType(std::vector<unsigned int> vec)
{
  ///////////////////////////////////////////
  // X -> A(mu)         B (mumu)  return  1
  // X ->   mu          B (mumu)  return  2
  // X ->   mu      mu    B (mu)  return  3
  // X -> A(mu)   B(mu)    C(mu)  return  4
  // X ->   mu    A(mu)    B(mu)  return  5
  // X ->   mu      mu       mu   return  6
  ////////////////////////////////////////////
  ///////////////////////////////////////////
  // X -> mu  Y -> mumu           return  7
  // X -> mu  Y -> mu Z -> mu     return  8
  // none above                   return -1
  ////////////////////////////////////////////



  int type(-1);

  std::vector<unsigned int> ancestors;
  ancestors = findSameAncestors(vec);

  if(ancestors.size() == 3) 
    {
      return ClassifyTypeI(vec);
    }
  else 
    {
      return ClassifyTypeII(vec);
    }

  return type;
}


int 
EventClassifier::ClassifyTypeII(std::vector<unsigned int> vec)
{

  ///////////////////////////////////////////
  // X -> mu  Y -> mumu           return  7
  // X -> mu  Y -> mu Z -> mu     return  8
  // none above                   return -1
  ////////////////////////////////////////////

  int type(-1);

  std::vector<unsigned int> ancestors;
  ancestors = findSameAncestors(vec);
  if(ancestors.size() == 1 || ancestors.size() == 2)
    {

      std::vector<unsigned int> A  = AllParentsOfParticle(vec.at(0));
      std::vector<unsigned int> B  = AllParentsOfParticle(vec.at(1));
      std::vector<unsigned int> C  = AllParentsOfParticle(vec.at(2));

      if(  Intersection(A,B).size() +  Intersection(A,C).size() +  Intersection(B,C).size() != 0 ) type = 7; // at least on of the intersectio in non-zero
      if(  Intersection(A,B).size() +  Intersection(A,C).size() +  Intersection(B,C).size() == 0 ) type = 8; // at least on of the intersectio in non-zero

    }


  return type;
}



int
EventClassifier::ClassifyTypeI(std::vector<unsigned int> vec)
{
  ///////////////////////////////////////////
  // X -> A(mu)         B (mumu)  return  1
  // X ->   mu          B (mumu)  return  2
  // X ->   mu      mu    B (mu)  return  3
  // X -> A(mu)   B(mu)    C(mu)  return  4
  // X ->   mu    A(mu)    B(mu)  return  5
  // X ->   mu      mu       mu   return  6
  // none above                   return -1  
  ////////////////////////////////////////////

  int type(-1);


  std::vector<unsigned int> parents;
  std::vector<unsigned int> ancestors;
  std::vector<unsigned int> parents_duplicates;

  ancestors = findSameAncestors(vec);
  if(ancestors.size() != 3) return type;


  for(auto i  :  vec)
    if(Ntp->MCParticle_hasMother(i))
      parents.push_back(Ntp->MCParticle_midx(i));

  parents_duplicates = findDuplicates(parents);

  if( parents_duplicates.size() == 1 )
    {
      if(  std::count(parents.begin(), parents.end(), ancestors.front()) == 0  )
	{
	  type = 1;
	  std::sort(parents.begin(), parents.end());
          if(parents.size() == 3)
	    if(  Ntp->MCParticle_hasMother( parents.at(2))  )
	      {
		if(Ntp->MCParticle_midx( parents.at(2) )  == (int)parents.at(0))      type = 2;
		
		std::vector<int> parents_triplet_has_mothers = triplet_has_mother(parents);
		if(std::count(parents_triplet_has_mothers.begin(), parents_triplet_has_mothers.end(), 1 ) == 3 )
		  if(Ntp->MCParticle_hasMother( Ntp->MCParticle_midx(parents.at(2))) && Ntp->MCParticle_hasMother( Ntp->MCParticle_midx(parents.at(1))))
		    if(parents.at(2) == parents.at(1) && Ntp->MCParticle_midx(Ntp->MCParticle_midx(parents.at(2))) == (int)parents.at(0) )type = 2;
		  
	      }

	}
      if(  std::count(parents.begin(), parents.end(), ancestors.front()) == 1  ) type = 2;
      if(  std::count(parents.begin(), parents.end(), ancestors.front()) == 2  ) type = 3;
      
    }


  if( parents_duplicates.size() == 0 )
    {

      if(  std::count(parents.begin(), parents.end(), ancestors.front()) == 0  ) 
	{
	  type = 4;

	  ///  Type 4 might be misclassified as type 1, check this case here.
	  std::sort(parents.begin(), parents.end()); 
	  if(parents.size() == 3)
	    if(  Ntp->MCParticle_hasMother( parents.at(2))  )
	      {
		if( Ntp->MCParticle_hasMother( Ntp->MCParticle_midx( parents.at(2) ) )  )
		  {
		    if(Ntp->MCParticle_midx( parents.at(2) )  == (int)parents.at(1))      type = 1;    
		    if( Ntp->MCParticle_hasMother( Ntp->MCParticle_midx( Ntp->MCParticle_midx( parents.at(2) ) ) ) )
		      if(   Ntp->MCParticle_midx( Ntp->MCParticle_midx( parents.at(2) ) )  ==                       (int)parents.at(1) ||
			    Ntp->MCParticle_midx( Ntp->MCParticle_midx( parents.at(2) ) )  ==  Ntp->MCParticle_midx( parents.at(1)) )    type = 1;
		    
		    if(Ntp->MCParticle_midx( parents.at(2) )                       == (int)parents.at(1) && 
		       Ntp->MCParticle_midx(Ntp->MCParticle_midx(parents.at(2) ) ) == (int)parents.at(0) && 
		       Ntp->MCParticle_midx(parents.at(1) )                        == (int)parents.at(0)) type = 2;
		        
		  }

		if(Ntp->MCParticle_hasMother( parents.at(1)))
		  if(Ntp->MCParticle_midx( parents.at(2)) == Ntp->MCParticle_midx( parents.at(1)) )                                  type = 1;
	      }

	  if(parents.size() == 2)// Rare case, but anyways change type to 1
	    if(  Ntp->MCParticle_hasMother( parents.at(1))  )
	      if( Ntp->MCParticle_hasMother( Ntp->MCParticle_midx( parents.at(1) ) )  )
		if( Ntp->MCParticle_hasMother( Ntp->MCParticle_midx( Ntp->MCParticle_midx( parents.at(1) ) ) ) )
		  if(   Ntp->MCParticle_midx( Ntp->MCParticle_midx( parents.at(1) ) )  == (int)parents.at(0)) type =1; 
	    

	}
      

      if(  std::count(parents.begin(), parents.end(), ancestors.front()) == 1  ) 
	{
	  type = 5;
	  std::sort(parents.begin(), parents.end());
          if(parents.size() == 3)
	    if(Ntp->MCParticle_hasMother( parents.at(2)) )
	      if(Ntp->MCParticle_midx( parents.at(2) ) ==  (int)parents.at(1) ) type = 2;
	}


    }

  if(parents_duplicates.size() == 3) type = 6; // This is typycally DecaysInFlights D-> K pi pi

  return type;
}


int
EventClassifier::LeParentPremier(std::vector<unsigned int> vec)
{

  int index(-1);
  std::vector<unsigned int> ancestors;
  ancestors = findSameAncestors(vec);

  if(ancestors.size() == 3)
    {
      if(std::count(ancestors.begin(), ancestors.end(), ancestors.front()) == (int)ancestors.size()) 
	index = ancestors.front();
    }

  return index;
}

std::vector<unsigned int> 
EventClassifier::findSameAncestors(std::vector<unsigned int> vec)// <--- It returns vector of size 3 with same lement - MC index of the first parent
{
  std::vector<unsigned int> out;
  int count_protons_childs(0);
  for(auto i  :  vec)
    if(Ntp->MCParticle_hasMother(i))
      {
	if(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(i)) == 2212)
	  {
	    count_protons_childs++;
	    if(  count_protons_childs == 3 && findDuplicates(vec).size()  != 0 )
	      {
		//		std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!  requires further classification " << std::endl;
		return out;  //if more that one particle is coming from proton break it;
	      }

	    if(count_protons_childs == 2 && findDuplicates(vec).size()  == 0 )
		return out;


	    ////////////////////////////////////////////
	    ///////////  these are going to be Type II 
	    std::vector<int> TripleHasMother = triplet_has_mother(vec);
	    if(count_protons_childs == 2 && findDuplicates(vec).size() ==1 && std::count(TripleHasMother.begin(), TripleHasMother.end(), 0 )==1 && 
	       std::count(TripleHasMother.begin(), TripleHasMother.end(), 1  ) == 2 ) return out;  

	    if(count_protons_childs == 1 && findDuplicates(vec).size()  == 1 && std::count(TripleHasMother.begin(), TripleHasMother.end(), 0 )==2 && 
	       std::count(TripleHasMother.begin(), TripleHasMother.end(), 1 )==1 ) return out;

	    out.push_back(i);
	  }
	else
	  {
	    out.push_back(Ntp->MCParticle_midx(i));
	  }
      }
    else if( Ntp->MCParticle_status(i) == 2) // keep the particle in the stack if it has not mother -> candidate to be to treee top
      {
	out.push_back(i);
      }
  

  int allThreeHave_mother(0);
  for(auto o : out) 
    {
      if(Ntp->MCParticle_hasMother(o))
	allThreeHave_mother+=1;
    }
  
  if(allThreeHave_mother == 0)
    {
      if(std::count(out.begin(), out.end(), out.front()) != (int)out.size()) 
	return findDuplicates(out);

      return out;
    }



  if( std::count(out.begin(), out.end(), out.front()) == (int)out.size() )
      return out;

  return findSameAncestors(out);
}


int 
EventClassifier::TypeI1_I2_pair_parent(std::vector<unsigned int> vec)
{

  if(ClassifyTypeI(vec) == 1 ||          // X -> A(mu)         B (mumu)
     ClassifyTypeI(vec) == 2    )        // X ->   mu          B (mumu) 
    {

      std::vector<unsigned int> parents;
      std::vector<unsigned int> parents_duplicates;

      for(auto i  :  vec)
	if(Ntp->MCParticle_hasMother(i))
	  parents.push_back(Ntp->MCParticle_midx(i));

      parents_duplicates = findDuplicates(parents);

      if(parents_duplicates.size() == 1)
	return parents_duplicates.at(0);
    }

  return -1;

}





//bool EventClassifier::CheckDoubleEvents(int run, int event) {
//  //static std::set<std::pair<int, int> > fRunEventPair;
//   //check whether pair is already in set and processed, otherwise insert it in set
//  return fRunEventPair.insert(std::make_pair(run,event)).second;
//} 
