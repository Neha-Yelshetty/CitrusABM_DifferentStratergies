#ifndef BIO_ABM
#define BIO_ABM
/********************************
 * INCLUDES
 * ******************************/
#include <iostream>
#include <vector>
#include <queue>
#include <array>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/assign.hpp>
#include <memory>
#include <list>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <cereal/archives/json.hpp>
#include <cereal/types/memory.hpp>
#include <ctime>
#include <chrono>
#include <atomic>
#include <random>
#include <numeric>
#include <stdlib.h>
#include <iostream>
#include "../headers/grove.hpp"

namespace bioABM {
    int getModelDay();
    int getModelDuration();

    int getNumRows();
    int getRowLength();
    
    // Returns age at coordinates, with an optional age differential in days
    int getAgeAt(int i, int j, int differential = 0);

    //returns an array of flush starting days (spring, summer, fall)
    int getFallStart();
    int getSummerStart();
    int getSpringStart();
    int countPsyllids();
    //Returns severity at coordinates
    //double getSeverityAt(int i, int j);

    //Returns true if tree is symptomatic
    bool isSymptomatic(int i, int j,Grove* agent);

    //Returns psyllids at coordinates
    double getPsyllidsAt(int i, int j);

    typedef boost::tuple<int, int> coord;
    //Functions to interface w/ from economic model
    bool rogueTreeAt(int i, int j,Grove* agent);

    void sprayTrees(double, std::vector<coord>,Grove* agent);

    bool isValidCoordinate(coord);

    bool isValidCoordinateAtgrove(coord,Grove* agent);

    void parseParameterFile(std::string);

    void advanceBiologicalModel(int farmid,Grove* agent);

    //void initalize(Grove* agent);
    
    //void advanceBiologicalModel_parallel();

    //bool isTreeAlive(int i, int j);

    void finishRun();

    void bioTestSuite();
    
    void setExperimentID(int);

    double getSeverityAtgrove(int i, int j,Grove* agent);

    bool isTreeAliveAtgrove(int i, int j,Grove* agent);


    /*******************************
 * FLUSH PATCH
 * The core grid unit, holds flush
 * shoots. The term flush patch
 * is used interchangebly with
 * tree in our discussions/paper
 * *****************************/
struct FlushPatch {
    int hlbseverityon = false;
    int age = 0;
    bool alive = true;
    int oldInfectedShoots = 0;
    int oldUninfectedShoots = 0;
    int numPsyllids_male = 0;
    int numPsyllids_female = 0;
    int numInfectedPsyllids_male = 0;
    int numInfectedPsyllids_female = 0;
    bool symptomatic = false;
    int daysInfected = 0;
    bool infectious = false;
    std::array<int,17> numNymphs;
    std::array<int,17> numInfectedNymphs;
    std::array<int, 30> numShoots;
    std::array<int, 30> numInfectedShoots;
    Grove* grove;


    void setGrove(Grove* g) {
        grove = g;
    }

    Grove* getGrove() {
        return grove;
    }

    //Constructor
    FlushPatch() {      
        for (int i = 0; i < 17; i++) {
            numNymphs[i] = 0;
            numInfectedNymphs[i] = 0;
        }
        for (int i = 0; i < 30; i++) {
            numShoots[i] = 0;
            numInfectedShoots[i] = 0;
        }
    }

    //Used on tree death
    void kill() {
        oldInfectedShoots = 0;
        oldUninfectedShoots = 0;
        numShoots.fill(0);
        numInfectedShoots.fill(0);
        clearPsyllids();
        alive = false;
        symptomatic = false;
    }

    //Debugging helper
    bool validate() {
        for (int i = 0; i < 17; i++) {
            if (numNymphs[i] < 0 || numInfectedNymphs[i] < 0) {
                return false;
            }
        }
        for (int i = 0; i < 30; i++) {
            if (numShoots[i] < 0 || numInfectedShoots[i] < 0) {
                return false;
            }
        }
        return true;
    }

    //Age getter
    int getAge() {
        return age;
    }

    //Psyllid getter
    int getNumPsyllids() {
        int numPsyllids = 0;
        numPsyllids += numPsyllids_male;
        numPsyllids += numPsyllids_female;
        numPsyllids += numInfectedPsyllids_male;
        numPsyllids += numInfectedPsyllids_female;
        //numPsyllids += accumulate(numNymphs.begin(), numNymphs.end(), 0);
        //numPsyllids += accumulate(numInfectedNymphs.begin(), numInfectedNymphs.end(), 0);
        return numPsyllids;
    }

    //Calculate HLB severity, which is measured as the proportion of infected flushes
    double getHLBSeverity() {
        
        int uninfected = std::accumulate(numShoots.begin(), numShoots.end(), 0);
        int infected = std::accumulate(numInfectedShoots.begin(), numInfectedShoots.end(), 0);
        double hlbNum = (double)infected + (double)oldInfectedShoots;
        double hlbDenom = (double)uninfected + (double)oldUninfectedShoots + (double)hlbNum;
       
      
            if (hlbDenom == 0 || hlbseverityon) {
            
                return 0;
            }
            else {
                assert((hlbNum / hlbDenom) >= 0 && (hlbNum / hlbDenom) <= 1);
                double severity = hlbNum / hlbDenom;
                return severity;
            }
        

    }

    //Remove all psyllids
    void clearPsyllids() {
        numPsyllids_female = 0;
        numPsyllids_male = 0;
        numInfectedPsyllids_female = 0;
        numInfectedPsyllids_male = 0;
        numNymphs.fill(0);
        numInfectedNymphs.fill(0);
        return;
    }

    //Wrapper used somewhere, unsure why this is here
    int getTotalPsyllids() {
        return getNumPsyllids();
    }

    //Helper to add a psyllid based on its characteristics
    void placePsyllid(bool female, bool adult, bool infected, int age = -1) {
        if (adult && female && infected) {
            numInfectedPsyllids_female++;
        }
        else if (adult && female && !infected) {
            numPsyllids_female++;
        }
        else if (adult && !female && infected) {
            numInfectedPsyllids_male++;
        }
        else if (adult && !female && !infected) {
            numPsyllids_male++;
        }
        else if (!adult && !infected) {
            numNymphs[age]++;
        }
        else if (!adult && infected) {
            numInfectedNymphs[age]++;
        }
    }
};


/********************************
 * FLUSH SHOOT
 * A single flush shoot, is part
 * of a greater flush patch
 * ******************************/
struct FlushShoot {
    int age = 0;
    bool infected = false;
    bool symptomatic = false;
    int numEggs = 0;
    int daysAsymptomatic = 0;
    bool alive = true;
    bool bark = false;
    int numPsyllids = 0;
    int numInfectedPsyllids = 0;
};

}

















#endif