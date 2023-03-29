#ifndef MEASUREMENTS_HH
#define MEASUREMENTS_HH

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
#include "./flags.h"
#include "./measurements.h"

using namespace plb;
using namespace plb::descriptors;
using namespace std;


template<typename T>
Measurements3D<T>::Measurements3D(MultiBlockLattice3D<T, NSDESCRIPTOR>& lattice_,
        MultiScalarField3D<int>& topology_, plint resolution_, double smoothFactor_,
        bool flatMode_, bool doHist_, double areaFromSTL)
    : lattice(lattice_),
      resolution(resolution_),
      Nx(lattice.getNx()),
      Ny(lattice.getNy()),
      Nz(lattice.getNz()),
      topology(topology_),
      kineticEnergyField(        MultiScalarField3D<T>(lattice)),
      densityField(              MultiScalarField3D<T>(lattice)),
      qcriterionField(           MultiScalarField3D<T>(lattice)),
      velocityField(             MultiTensorField3D<T,3>(lattice)),
      velocityTransField(        MultiTensorField3D<T,3>(lattice)),
      strainRateField(           MultiTensorField3D<T,6>(lattice)),
      strainRateTransField(      MultiTensorField3D<T,6>(lattice)),
      vorticityField(            MultiTensorField3D<T,3>(lattice)),
      vorticityTransField(       MultiTensorField3D<T,3>(lattice)),
      shearStressField(          MultiTensorField3D<T,6>(lattice)),
      shearStressTransField(     MultiTensorField3D<T,6>(lattice)),
      shearStressNormField(      MultiScalarField3D<T>(lattice)),
      hist( gsl_histogram_alloc(100)),
      smoothFactor(smoothFactor_),
      densitySum(0.),
      kineticEnergySum(0.),
      avgCount(0),
      transCount(0),
      fluidRange( inRange<int>(flagBB,  flagBuffer)      ),
      wallRange(  inRange<int>(flagBB,  flagSedimenting) ),
      serosRange( inRange<int>(flagWet, flagWet)         ),
      fluidCellsInit(0),
      volumeToSurfaceRatioInit(0.),
      erosThreshold(0.),
      sediThreshold(0.),
      fluidCellsIncrease(0),
      flatMode(flatMode_),
      doHist(doHist_),
      bestSeros(0),
      minRhoBar(0)
{

    std::vector<std::string> xAxis, flow, seros;
    xAxis.push_back("iT");
    xAxis.push_back("time");
    xAxis.push_back("seros"); // number of seros steps performed

    flow.push_back("rhobarInlet");
    flow.push_back("rhobarGlobal");
    flow.push_back("kineticEnergy");
    flow.push_back("shearStressMin");
    flow.push_back("shearStressMax");
    flow.push_back("shearStressDiff");
    flow.push_back("histMean");
    flow.push_back("histSigma");
    flow.push_back("histSkew");
    flow.push_back("histKurt");

    seros.push_back("rhobarInletAvg");
    seros.push_back("rhobarInletAvgSmooth");
    seros.push_back("rhobarInletAvgSmoothGrad");
    seros.push_back("rhobarInletAvgSmoothGradSmooth");
    seros.push_back("rhobarGlobalAvg");
    seros.push_back("rhobarGlobalAvgSmooth");
    seros.push_back("rhobarGlobalAvgSmoothGrad");
    seros.push_back("rhobarGlobalAvgSmoothGradSmooth");
    seros.push_back("kineticEnergyAvg");
    seros.push_back("kineticEnergyAvgSmooth");
    seros.push_back("kineticEnergyAvgSmoothGrad");
    seros.push_back("diameter");
    seros.push_back("fluidCells");
    seros.push_back("wallCells");
    seros.push_back("volumeToSurfaceRatio");

    dataNamesFlow.reserve( xAxis.size() + flow.size() );
    dataNamesFlow.insert( dataNamesFlow.end(), xAxis.begin(), xAxis.end() );
    dataNamesFlow.insert( dataNamesFlow.end(), flow.begin(),  flow.end() );

    dataNamesSeros.reserve( xAxis.size() + seros.size() );
    dataNamesSeros.insert( dataNamesSeros.end(), xAxis.begin(),  xAxis.end() );
    dataNamesSeros.insert( dataNamesSeros.end(), seros.begin(),  seros.end() );

    dataNames.reserve( dataNamesFlow.size() + dataNamesSeros.size() );
    dataNames.insert( dataNames.end(), dataNamesFlow.begin(),  dataNamesFlow.end() );
    dataNames.insert( dataNames.end(), dataNamesSeros.begin(), dataNamesSeros.end() );

    for (unsigned int i=0; i<dataNames.size(); i++) {
        data.emplace(dataNames[i], 0.);
    }

    data["diameter"] = 1.;
    data["fluidCells"] = count(topology, fluidRange);
    fluidCellsInit = data["fluidCells"];
    data["wallCells"]  = count(topology, wallRange);
    if (areaFromSTL != 0) {
        data["volumeToSurfaceRatio"] = data["fluidCells"] / areaFromSTL;
    } else
        data["volumeToSurfaceRatio"] = data["fluidCells"] / data["wallCells"];
    volumeToSurfaceRatioInit = data["volumeToSurfaceRatio"];
}

template<typename T>
const string Measurements3D<T>::getDataHeaderFlow() {
    std::string dataHeader = createStringFromVector(dataNamesFlow, " ");
    return dataHeader;
}

template<typename T>
const string Measurements3D<T>::getDataHeaderSeros() {
    std::string dataHeader = createStringFromVector(dataNamesSeros, " ");
    return dataHeader;
}

template<typename T>
const string Measurements3D<T>::getDataFlow() {
    std::string recentData = createStringFromMap(data, dataNamesFlow, " ");
    return recentData;
}

template<typename T>
const string Measurements3D<T>::getDataSeros() {
    std::string recentData = createStringFromMap(data, dataNamesSeros, " ");
    return recentData;
}

template<typename T>
bool Measurements3D<T>::checkBestSeros() {
    bool newBestSeros = false;
    if (data["rhobarInletAvg"] < minRhoBar || data["seros"] == 0) {
        minRhoBar = data["rhobarInletAvg"];
        bestSeros = data["seros"];
        newBestSeros = true;
    }
    return newBestSeros;
}

template<typename T>
void Measurements3D<T>::resetAvg() {
    data["rhobarInletAvg"]=0;
    data["rhobarGlobalAvg"]=0;
    data["kineticEnergyAvg"]=0;
    velocityTransField = MultiTensorField3D<T,3>(lattice);
    strainRateTransField = MultiTensorField3D<T,6>(lattice);
    vorticityTransField = MultiTensorField3D<T,3>(lattice);
    shearStressTransField = MultiTensorField3D<T,6>(lattice);
    avgCount=0;
    transCount=0;
}

template<typename T>
const string Measurements3D<T>::getHistogramData() {
    std::stringstream histSS;
    histSS << "max value: " << gsl_histogram_max(hist.get()) << endl;
    for (unsigned int i=0; i<gsl_histogram_bins(hist.get()); i++)
        histSS << gsl_histogram_get(hist.get(), i) << endl;
    return histSS.str();
}

template<typename T>
const string Measurements3D<T>::getDistributionData() {
    std::stringstream distrSS;
    for (unsigned int i=0; i<wallStressVec.size(); i++)
        distrSS << wallStressVec.at(i) << endl;
    return distrSS.str();
}

template<typename T>
void Measurements3D<T>::load(std::vector<T> recentData) {
    if (recentData.size() == 0) return;
    changeElementsInMap(data, dataNames, recentData);
}

template<typename T>
const T Measurements3D<T>::getRhobarGlobalAvgSmoothGradSmooth() {
    return data["rhobarGlobalAvgSmoothGradSmooth"]; }

template<typename T>
const T Measurements3D<T>::getRhobarInletAvgSmoothGradSmooth() {
    return data["rhobarInletAvgSmoothGradSmooth"]; }

template<typename T>
const T Measurements3D<T>::getKineticEnergyAvgSmoothGrad() {
    return data["kineticEnergyAvgSmoothGrad"]; }

template<typename T>
MultiScalarField3D<T>& Measurements3D<T>::getDensityField() {
    return densityField; }

template<typename T>
MultiScalarField3D<T>& Measurements3D<T>::getQcriterionField() {
    qcriterionField = *computeQcriterion(vorticityTransField, strainRateTransField);
    return qcriterionField; }

template<typename T>
MultiTensorField3D<T,3>& Measurements3D<T>::getVelocityField() {
    return velocityField; }

template<typename T>
MultiTensorField3D<T,3>& Measurements3D<T>::getVelocityTransField() {
    return velocityTransField; }

template<typename T>
MultiTensorField3D<T,6>& Measurements3D<T>::getStrainRateField() {
    return strainRateField; }

template<typename T>
MultiTensorField3D<T,6>& Measurements3D<T>::getStrainRateTransField() {
    return strainRateTransField; }

template<typename T>
MultiTensorField3D<T,3>& Measurements3D<T>::getVorticityField() {
    return vorticityField; }

template<typename T>
MultiTensorField3D<T,3>& Measurements3D<T>::getVorticityTransField() {
    return vorticityTransField; }

template<typename T>
MultiTensorField3D<T,6>& Measurements3D<T>::getShearStressField() {
    return shearStressField; }

template<typename T>
MultiTensorField3D<T,6>& Measurements3D<T>::getShearStressTransField() {
    return shearStressTransField; }

template<typename T>
MultiScalarField3D<T>& Measurements3D<T>::getShearStressNormField() {
    return shearStressNormField; }

template<typename T>
T Measurements3D<T>::calculateDiameter(MultiScalarField3D<int>& topo, double areaFromSTL) {
    data["fluidCells"] = count(topo, fluidRange);
    data["wallCells"]  = count(topo, wallRange);
    if (areaFromSTL != 0) {
        data["volumeToSurfaceRatio"] = data["fluidCells"] / areaFromSTL;
    } else
        data["volumeToSurfaceRatio"] = data["fluidCells"] / data["wallCells"];
    data["diameter"] = data["volumeToSurfaceRatio"] / volumeToSurfaceRatioInit;
    return data["diameter"];
}

template<typename T>
plint Measurements3D<T>::getFluidCellsIncrease() {
    plint fluidCellsIncrease = data["fluidCells"] - fluidCellsInit;
    if (flatMode)
        fluidCellsIncrease /= Nz; // -> per layer
    return fluidCellsIncrease;
}

template<typename T>
void Measurements3D<T>::calculateWallStress(MultiScalarField3D<int>& topo) {

    wallStressVec.clear();
    MaskedBoxScalarListFunctional3D<T, inRange<int>> shearStressListFunctional(serosRange);
    if (flatMode) {
        applyProcessingFunctional(
            shearStressListFunctional, layer(Nx,Ny,0), shearStressNormField, topology);
    } else {
        applyProcessingFunctional(
            shearStressListFunctional, lattice.getBoundingBox(), shearStressNormField, topology);
    }
    wallStressVec    = shearStressListFunctional.getListScalar();
    sort(wallStressVec.begin(), wallStressVec.end());
}

template<typename T>
void Measurements3D<T>::calculateWallStressHist() {
    gsl_histogram_reset(hist.get());
    if (wallStressVec.back() > 0.) {
        T tmp_maxWallStress = wallStressVec.back() * 1.001;
        histFromVect(hist.get(), 0., tmp_maxWallStress, wallStressVec);

        data["histMean"]  = gsl_histogram_mean(hist.get());
        data["histSigma"] = gsl_histogram_sigma(hist.get());
        data["histSkew"]  = gsl_stats_skew_m_sd(&wallStressVec[0], 1, wallStressVec.size(), data["histMean"], data["histSigma"]);
        data["histKurt"]  = gsl_stats_kurtosis_m_sd(&wallStressVec[0], 1, wallStressVec.size(), data["histMean"], data["histSigma"]);
    }
}

template<typename T>
void Measurements3D<T>::takeMeasurements(plint iT, T time, Box3D inlet, bool trans) {

    data["iT"]  = iT;
    data["time"]  = time;

    // fields:
    densityField       = *computeDensity(lattice);
    velocityField      = *computeVelocity(lattice);
    strainRateField    = *computeStrainRate(velocityField);
    vorticityField     = *computeVorticity(velocityField);
    shearStressField   = *computeShearStress(lattice);
    kineticEnergyField = *computeKineticEnergy(lattice);

    // density (equivalent to pressure)
    T densityInletBulk = computeAverage(densityField, topology, flagBulk, inlet);
    T densityInletBuff = computeAverage(densityField, topology, flagBuffer, inlet);
    if (densityInletBulk != 0. && densityInletBuff != 0.)
        data["rhobarInlet"] = (densityInletBulk + densityInletBuff) / 2 -1;
    else
        data["rhobarInlet"] =  densityInletBulk + densityInletBuff -1;
    densitySum               = computeAverage(densityField, topology, flagWet,  lattice.getBoundingBox());
    densitySum              += computeAverage(densityField, topology, flagBulk, lattice.getBoundingBox());
    data["rhobarGlobal"]     = densitySum / 2 -1;

    // kinetic energy
    kineticEnergySum         = computeAverage(kineticEnergyField, topology, flagWet, lattice.getBoundingBox());
    kineticEnergySum        += computeAverage(kineticEnergyField, topology, flagBulk, lattice.getBoundingBox());
    data["kineticEnergy"] = kineticEnergySum / 2;

    // add to averaging between seros steps:
    avgCount++;
    addToAverage(data["rhobarInletAvg"],  data["rhobarInlet"],  avgCount);
    addToAverage(data["rhobarGlobalAvg"], data["rhobarGlobal"], avgCount);
    addToAverage(data["kineticEnergyAvg"], data["kineticEnergy"], avgCount);

    if (trans) {
        transCount++;
        applyProcessingFunctional(
            new UpdateAveTensorTransientStatistics3D<T,3>(transCount),
            lattice.getBoundingBox(), velocityField, velocityTransField);
        applyProcessingFunctional(
            new UpdateAveTensorTransientStatistics3D<T,6>(transCount),
            lattice.getBoundingBox(), strainRateField, strainRateTransField);
        applyProcessingFunctional(
            new UpdateAveTensorTransientStatistics3D<T,3>(transCount),
            lattice.getBoundingBox(), vorticityField, vorticityTransField);
        applyProcessingFunctional(
            new UpdateAveTensorTransientStatistics3D<T,6>(transCount),
            lattice.getBoundingBox(), shearStressField, shearStressTransField);
        shearStressNormField = *computeSymmetricTensorNorm(shearStressTransField);
    } else {
        velocityTransField = velocityField;
        strainRateTransField = strainRateField;
        vorticityTransField = vorticityField;
        shearStressTransField = shearStressField;
        shearStressNormField = *computeSymmetricTensorNorm(shearStressField);
    }

    // shear stress related:
    calculateWallStress(topology);
    if (doHist)
        calculateWallStressHist();
    data["shearStressMax"]   = wallStressVec.back();
    data["shearStressMin"]   = wallStressVec.front();
    data["shearStressDiff"]  = data["shearStressMax"]-data["shearStressMin"];
}


template<typename T>
void Measurements3D<T>::calculateSmoothingAndGradient() {

    // prepare gradient calculation (save last values in gradient variables):
    if (data["seros"] > 0) {
        data["rhobarInletAvgSmoothGrad"]  = -1 * data["rhobarInletAvgSmooth"];
        data["rhobarGlobalAvgSmoothGrad"]  = -1 * data["rhobarGlobalAvgSmooth"];
        data["kineticEnergyAvgSmoothGrad"] = -1 * data["kineticEnergyAvgSmooth"];
    }

    // exponential smoothing for current values
    // http://www.inventoryops.com/articles/exponential_smoothing.htm
    if (data["seros"] == 0) {
        data["rhobarInletAvgSmooth"]  = data["rhobarInletAvg"];
        data["rhobarGlobalAvgSmooth"]  = data["rhobarGlobalAvg"];
        data["rhobarInletAvgSmoothGradSmooth"]  = data["rhobarInletAvgSmoothGrad"];
        data["rhobarGlobalAvgSmoothGradSmooth"]  = data["rhobarGlobalAvgSmoothGrad"];
        data["kineticEnergyAvgSmooth"] = data["kineticEnergyAvg"];
    } else {
        data["rhobarInletAvgSmooth"]  = (data["rhobarInletAvg"]*smoothFactor) + (data["rhobarInletAvgSmooth"]*(1-smoothFactor));
        data["rhobarGlobalAvgSmooth"]  = (data["rhobarGlobalAvg"]*smoothFactor) + (data["rhobarGlobalAvgSmooth"]*(1-smoothFactor));
        data["kineticEnergyAvgSmooth"] = (data["kineticEnergyAvg"] * smoothFactor) + (data["kineticEnergyAvgSmooth"] * (1-smoothFactor));
    }

    // finish gradient calculation:
    if (data["seros"] > 0) {
        data["rhobarInletAvgSmoothGrad"]  += data["rhobarInletAvgSmooth"];
        data["rhobarGlobalAvgSmoothGrad"]  += data["rhobarGlobalAvgSmooth"];
        data["kineticEnergyAvgSmoothGrad"] += data["kineticEnergyAvgSmooth"];
    }

    // now apply smoothing to newly calculated gradients:
    if (data["seros"] != 0) {
        data["rhobarInletAvgSmoothGradSmooth"]  = (data["rhobarInletAvgSmoothGrad"]*smoothFactor) + (data["rhobarInletAvgSmoothGradSmooth"]*(1-smoothFactor));
        data["rhobarGlobalAvgSmoothGradSmooth"]  = (data["rhobarGlobalAvgSmoothGrad"]*smoothFactor) + (data["rhobarGlobalAvgSmoothGradSmooth"]*(1-smoothFactor));
    }
}

template<typename T>
const int Measurements3D<T>::getSerosCount() {
    return data["seros"];
}

template<typename T>
void Measurements3D<T>::increaseSerosCount() {
    data["seros"]++;
}

template<typename T>
T Measurements3D<T>::getErosThreshold(T fraction) {

    if (fraction <= 0.) {
        erosThreshold = wallStressVec.back()+1.;
        return erosThreshold;
    }

    if (fraction >= 1.) {
        erosThreshold = wallStressVec.front();
        return erosThreshold;
    }

    plint numberOfErosValues = (T)fraction * (T)wallStressVec.size();
    if (numberOfErosValues == 0)
        erosThreshold = wallStressVec.back();
    else
        erosThreshold = wallStressVec[wallStressVec.size()-numberOfErosValues];

    //pcout << "numberOfErosValues = " << numberOfErosValues << endl;
    return erosThreshold;
}

template<typename T>
T Measurements3D<T>::getFractionEroded() {

    plint fluidCellsAfterEros = count(topology, fluidRange);
    plint numberOfErodedCells = fluidCellsAfterEros - data["fluidCells"];
    if (flatMode)
        numberOfErodedCells /= Nz; // -> = numberOfErodedCells per layer
    T fractionEroded = (T)numberOfErodedCells / wallStressVec.size();
    //pcout << "numberOfErodedCells = " << numberOfErodedCells << endl;
    return fractionEroded;
}

template<typename T>
T Measurements3D<T>::getSediThreshold(T fraction, plint compensate) {

    plint numberOfSediValues = 0;
    numberOfSediValues += compensate;

    numberOfSediValues += (T)fraction * (T)wallStressVec.size();
    if (numberOfSediValues > 0) {
        if ((unsigned int) numberOfSediValues < wallStressVec.size())
            sediThreshold = wallStressVec[numberOfSediValues-1];
        else
            sediThreshold = wallStressVec.back();
    } else
        sediThreshold = 0.;

    //pcout << "numberOfSediValues = " << numberOfSediValues << endl;
    return sediThreshold;
}


#endif
