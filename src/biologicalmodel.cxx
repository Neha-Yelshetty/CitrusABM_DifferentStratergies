#include "biologicalmodel.hpp"
#include "../headers/behavior.hpp"
#include "bioABM.h"

//Behavior* behaviorPatterns[3];

void BiologicalModel::biologicalmodel(int farmid,Grove* g) {
    bioABM::advanceBiologicalModel(farmid,g);  
}

double BiologicalModel::biohlbseverity(int i, int j,Grove* g){
  return bioABM::getSeverityAtgrove(i, j,g);
}

bool BiologicalModel::bioisTreeAlive (int i, int j,Grove* g){
   return bioABM::isTreeAliveAtgrove(i,j,g);
}

/*void BiologicalModel::initialize(Grove* g){
     grove = g; 
    bioABM :: initalize(grove);
 }*/