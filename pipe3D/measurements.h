#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "../palabos/src/palabos3D.h"
#include "../palabos/src/palabos3D.hh"   // include full template code

#include <string>
#include <vector>
#include <unordered_map>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics.h>

#include "../palabos/tools/tools.h"
#include "../palabos/tools/tools.hh"
#include "../palabos/tools/tools_plb.h"
#include "../palabos/tools/tools_plb.hh"
#include "../palabos/tools/tools_gls.h"
#include "../palabos/tools/tools_gls.hh"


using namespace plb;
using namespace plb::descriptors;
using namespace std;


template<typename T>
class Measurements3D {
    public:
        Measurements3D(MultiBlockLattice3D<T, NSDESCRIPTOR> &lattice_,
                MultiScalarField3D<int>& topology_, plint resolution_,
                double smoothFactor_, bool flatMode_, bool doHist_, double area=0);
        void load(std::vector<T> recentData);
        T    calculateDiameter(MultiScalarField3D<int>& topo, double area=0);
        void calculateWallStress(MultiScalarField3D<int>& topo);
        void calculateWallStressHist();
        void takeMeasurements(plint iT, T time, Box3D inlet, bool trans);
        T getErosThreshold(T fraction);
        T getFractionEroded();
        T getSediThreshold(T fraction, plint compensate=0);
        plint getFluidCellsIncrease();
        const string getDataHeaderFlow();
        const string getDataHeaderSeros();
        const string getDataFlow();
        const string getDataSeros();
        bool checkBestSeros();
        void resetAvg();
        const string getHistogramData();
        const string getDistributionData();
        const T getKineticEnergyAvgSmoothGrad();
        const T getRhobarGlobalAvgSmoothGradSmooth();
        const T getRhobarInletAvgSmoothGradSmooth();
        const int getSerosCount();
        void increaseSerosCount();
        void calculateSmoothingAndGradient();
        MultiScalarField3D<T>& getDensityField();
        MultiScalarField3D<T>& getQcriterionField();
        MultiTensorField3D<T,3>& getVelocityField();
        MultiTensorField3D<T,3>& getVelocityTransField();
        MultiTensorField3D<T,6>& getStrainRateField();
        MultiTensorField3D<T,6>& getStrainRateTransField();
        MultiTensorField3D<T,3>& getVorticityField();
        MultiTensorField3D<T,3>& getVorticityTransField();
        MultiTensorField3D<T,6>& getShearStressField();
        MultiTensorField3D<T,6>& getShearStressTransField();
        MultiScalarField3D<T>& getShearStressNormField();
    private:
        MultiBlockLattice3D<T, NSDESCRIPTOR> &lattice;
        const plint resolution;
        const plint Nx, Ny, Nz;
        MultiScalarField3D<int> &topology;
        MultiScalarField3D<T>   kineticEnergyField;
        MultiScalarField3D<T>   densityField;
        MultiScalarField3D<T>   qcriterionField;
        MultiTensorField3D<T,3> velocityField;
        MultiTensorField3D<T,3> velocityTransField;
        MultiTensorField3D<T,6> strainRateField;
        MultiTensorField3D<T,6> strainRateTransField;
        MultiTensorField3D<T,3> vorticityField;
        MultiTensorField3D<T,3> vorticityTransField;
        MultiTensorField3D<T,6> shearStressField;
        MultiTensorField3D<T,6> shearStressTransField;
        MultiScalarField3D<T>   shearStressNormField;

        unique_ptr<gsl_histogram> hist;
        double smoothFactor;

        std::unordered_map<std::string, T> data;
        std::vector<std::string> dataNames;
        std::vector<std::string> dataNamesFlow;
        std::vector<std::string> dataNamesSeros;
        //std::vector<T> wallStressVec;
        std::vector<T> wallStressVec;
        T densitySum;
        T kineticEnergySum;
        plint avgCount, transCount;
        inRange<int> fluidRange;
        inRange<int> wallRange;
        inRange<int> serosRange;
        plint fluidCellsInit;
        T volumeToSurfaceRatioInit;

        T erosThreshold;
        T sediThreshold;
        plint fluidCellsIncrease;
        const bool flatMode, doHist;
        int bestSeros;
        T minRhoBar;
};


#endif
