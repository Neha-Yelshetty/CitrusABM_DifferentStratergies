#ifndef BIOLOGICALMODEL_HPP
#define BIOLOGICALMODEL_HPP
#include "grove.hpp"
#include "math.h"
#include "planningFunc.hpp"
#include "behavior.hpp"
#include "bioABM.h"
#include<vector>

using namespace std;
class Behavior;

class BiologicalModel {
public:
     std::vector<Behavior*> behaviorPatterns;
     Grove* grove;
     BiologicalModel() : grove(nullptr) {}

     void biologicalmodel(int farmid,Grove* g); 
     double biohlbseverity(int i,int j,Grove* g);
     bool bioisTreeAlive(int i,int j,Grove* g);
    //void initialize(Grove* g);
};

#endif