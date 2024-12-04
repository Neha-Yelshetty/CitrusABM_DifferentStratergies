#include "../headers/grove.hpp"
#include "../headers/bioABM.h"
/*
Grove::Grove() {
    //Please dont use this
    setAge(0);
}*/



Grove::Grove() : bioModel(nullptr) {}

Grove::Grove(Commodity crop, bool agency, int i_lb, int i_ub, int j_lb, int j_ub,int farmid)
: crop(crop), agency(agency), farmid(farmid) {
    this->crop = crop;
    setAgency(agency);
    ibounds[0] = i_lb;
    ibounds[1] = i_ub;
    jbounds[0] = j_lb;
    jbounds[1] = j_ub;
    this->farmid = farmid;

    bioModel = new BiologicalModel();
    initializeLattice(i_ub, j_ub);
}

int* Grove::getIBounds() {
    return this->ibounds;
}

int* Grove::getJBounds() {
    return this->jbounds;
}

//set agency
void Grove::setAgency(bool agency) {
    this->agency = agency;
}

void Grove::initializeLattice(int numRows, int rowLength) {
    lattice.resize(numRows);  // Resize to the number of rows
    for (int i = 0; i < numRows; i++) {
        std::vector<bioABM::FlushPatch> row(rowLength);  // Create a row with pre-allocated patches
        for (int j = 0; j < rowLength; j++) {
            row[j] = bioABM::FlushPatch();  // Explicitly initialize each patch
        }
        lattice[i] = row;  // Assign the row to the lattice
    }
}

