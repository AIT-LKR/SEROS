#define NSDESCRIPTOR D3Q19Descriptor

#include "../palabos/src/palabos3D.h"
#include "../palabos/src/palabos3D.hh"   // include full template code
#include "../palabos/tools/tools.h"
#include "../palabos/tools/tools.hh"
#include "../palabos/tools/tools_files.h"
#include "../palabos/tools/tools_files.hh"
#include "../palabos/tools/tools_plb.h"
#include "../palabos/tools/tools_plb.hh"
#include "./simulationSettings.h"
#include "./simulationSettings.hh"
#include "./measurements.h"
#include "./measurements.hh"
#include "./geometry.h"
#include "./geometry.hh"
#include "./flags.h"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;


bool correctArguments(const int argc, char const* const* argv) {
    if (argc <  2)  return false;
    if (argc == 3)  return false;
    if (argc >  4)  return false;
    return true;
}

class ppmWriter {
public:
    ppmWriter(Measurements3D<T>& measure_,
              const T rhoMax_, const T tauMax_, const T velMax_,
              const bool showCurrent_, const bool tauScaled_)
        : measure(measure_),
          rhoMax(rhoMax_),
          tauMax(tauMax_),
          velMax(velMax_),
          showCurrent(showCurrent_),
          tauScaled(tauScaled_),
          imageWriterEarth(  ImageWriter<T>("earth") ),
          imageWriterLeeloo( ImageWriter<T>("leeloo") )
    { }
    void operator() (plint iter, plint iterSeros)
    {
        ScalarField2D<T> density          = extractMiddleLayer(measure.getDensityField());
        //ScalarField2D<T> qcriterion       = extractMiddleLayer(measure.getQcriterionField());
        ScalarField2D<T> velocity         = extractMiddleLayer(*computeNorm(measure.getVelocityField()));
        ScalarField2D<T> velocityTrans    = extractMiddleLayer(*computeNorm(measure.getVelocityTransField()));
        //ScalarField2D<T> strainRate       = extractMiddleLayer(*computeSymmetricTensorNorm(measure.getStrainRateField()));
        //ScalarField2D<T> strainRateTrans  = extractMiddleLayer(*computeSymmetricTensorNorm(measure.getStrainRateTransField()));
        //ScalarField2D<T> vorticity        = extractMiddleLayer(*computeNorm(measure.getVorticityField()));
        ScalarField2D<T> vorticityTrans   = extractMiddleLayer(*computeNorm(measure.getVorticityTransField()));
        ScalarField2D<T> shearStress      = extractMiddleLayer(*computeSymmetricTensorNorm(measure.getShearStressField()));
        ScalarField2D<T> shearStressTrans = extractMiddleLayer(*computeSymmetricTensorNorm(measure.getShearStressTransField()));

        imageWriterEarth.writePpm(createLongFileName("rho", iter, 8, "seros", iterSeros),
                density, 0.99, rhoMax);
        //imageWriterEarth.writeScaledPpm(createLongFileName("q", iter, 8, "seros", iterSeros),
        //        qcriterion);
        imageWriterLeeloo.writePpm(createLongFileName("uTrans", iter, 8, "seros", iterSeros),
                velocityTrans, 0, velMax);
        //imageWriterLeeloo.writeScaledPpm(createLongFileName("strainTransScaled", iter, 8, "seros", iterSeros),
        //        strainRateTrans);
        imageWriterLeeloo.writeScaledPpm(createLongFileName("vortTransScaled", iter, 8, "seros", iterSeros),
                vorticityTrans);
        if (tauScaled) {
            imageWriterLeeloo.writeScaledPpm(createLongFileName("tauTransScaled", iter, 8, "seros", iterSeros),
                    shearStressTrans);
        } else {
            imageWriterLeeloo.writePpm(createLongFileName("tauTrans", iter, 8, "seros", iterSeros),
                    shearStressTrans, 0, tauMax);
        }

        if (showCurrent) {
            imageWriterLeeloo.writePpm(createLongFileName("u", iter, 8, "seros", iterSeros),
                    velocity, 0, velMax);
            //imageWriterLeeloo.writeScaledPpm(createLongFileName("strainScaled", iter, 8, "seros", iterSeros),
            //        strainRate);
            //imageWriterLeeloo.writeScaledPpm(createLongFileName("vortScaled", iter, 8, "seros", iterSeros),
            //        vorticity);
            if (tauScaled) {
                imageWriterLeeloo.writeScaledPpm(createLongFileName("tauScaled", iter, 8, "seros", iterSeros),
                        shearStress);
            } else {
                imageWriterLeeloo.writePpm(createLongFileName("tau", iter, 8, "seros", iterSeros),
                        shearStress, 0, tauMax);
            }
        }
        return;
    }
private:
    Measurements3D<T> &measure;
    T rhoMax, tauMax, velMax;
    bool showCurrent, tauScaled;
    ImageWriter<T> imageWriterEarth;
    ImageWriter<T> imageWriterLeeloo;
};

void updateBoundaryVelocity(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice,
                            SimulationSettings<T>& sims, plint iter)
{
    T uNow;
    if (sims.inverseFlow) {
        uNow = -1*rampVelocity(sims.uInlet, iter, sims.rampT);
        if (sims.flatMode) {
            setBoundaryVelocity(NSlattice, sims.getOutletWindow(),
                    PoiseuilleVelocity<plint,T>(-1*uNow, sims.uDev,
                    sims.directionOfOutlet, sims.getOutletCentre(), sims.outletRadius+1));
            if (sims.shape == "T")
                setBoundaryVelocity(NSlattice, sims.getOutletTwoWindow(),
                        PoiseuilleVelocity<plint,T>(-1*uNow, sims.uDev,
                        sims.directionOfOutlet, sims.getOutletCentre(), sims.outletRadius+1));
            if (sims.shape == "Y")
                setBoundaryVelocity(NSlattice, sims.getOutletTwoWindow(),
                        PoiseuilleVelocity<plint,T>(-1*uNow, sims.uDev,
                        sims.directionOfOutlet, sims.getOutletCentre(false), sims.outletRadius+1));
        } else {
            setBoundaryVelocity(NSlattice, sims.getOutletWindow(),
                    PoiseuilleVelocity<plint,T>(-1*uNow, sims.uDev,
                    sims.directionOfOutlet, sims.getOutletCentre(), sims.Nz/2, sims.outletRadius+1));
            if (sims.shape == "T")
                setBoundaryVelocity(NSlattice, sims.getOutletTwoWindow(),
                        PoiseuilleVelocity<plint,T>(-1*uNow, sims.uDev,
                        sims.directionOfOutlet, sims.getOutletCentre(), sims.Nz/2, sims.outletRadius+1));
            if (sims.shape == "Y")
                setBoundaryVelocity(NSlattice, sims.getOutletTwoWindow(),
                        PoiseuilleVelocity<plint,T>(-1*uNow, sims.uDev,
                        sims.directionOfOutlet, sims.getOutletCentre(false), sims.Nz/2, sims.outletRadius+1));
        }
    } else {
        uNow = rampVelocity(sims.uInlet, iter, sims.rampT);
        if (sims.flatMode) {
            setBoundaryVelocity(NSlattice, sims.getInletWindow(),
                    PoiseuilleVelocity<plint,T>(uNow, sims.uDev,
                    sims.directionOfInlet, sims.inletCentre, sims.inletRadius+1));
            if (sims.shape == "Ti")
                setBoundaryVelocity(NSlattice, sims.getInletTwoWindow(),
                        PoiseuilleVelocity<plint,T>(-1*uNow, sims.uDev,
                        sims.directionOfInlet, sims.inletCentre, sims.inletRadius+1));
        } else {
            setBoundaryVelocity(NSlattice, sims.getInletWindow(),
                    PoiseuilleVelocity<plint,T>(uNow, sims.uDev,
                    sims.directionOfInlet, sims.inletCentre, sims.Nz/2, sims.outletRadius+1, sims.inletHeight/2+1));
            if (sims.shape == "Ti")
                setBoundaryVelocity(NSlattice, sims.getInletTwoWindow(),
                        PoiseuilleVelocity<plint,T>(uNow, sims.uDev,
                        sims.directionOfInlet, sims.inletCentre, sims.Nz/2, sims.inletRadius+1));
        }
    }
}

void writeTopoVTK(Geometry<T>& geo, plint iter, plint iterSeros)
{
    VtkImageOutput3D<int> vtkOutInt(createLongFileName("vtkInt", iter, 8, "seros", iterSeros));
    vtkOutInt.writeData<int>(geo.getTopology(), "topology");
}

void writeFlowVTK(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice,
                  Measurements3D<T>& measure, plint iter, plint iterSeros)
{
    VtkImageOutput3D<T> vtkOut(createLongFileName("vtk", iter, 8, "seros", iterSeros));
    vtkOut.writeData<float>(measure.getDensityField(), "density");
    vtkOut.writeData<float>(measure.getQcriterionField(), "qcriterion");
    vtkOut.writeData<3,float>(measure.getVelocityTransField(), "velocity");
    //vtkOut.writeData<3,float>(measure.getVelocityTransField(), "velocity", dx/dt);
    vtkOut.writeData<6,float>(measure.getStrainRateTransField(), "strainRate");
    vtkOut.writeData<3,float>(measure.getVorticityTransField(), "vorticity");
    vtkOut.writeData<6,float>(measure.getShearStressTransField(), "shearStress");
    //vtkOut.writeData<float>(*computeKineticEnergy(NSlattice), "kineticEnergy");
    vtkOut.writeData<float>(*computeOmega(NSlattice), "viscosity");
}

void loadOrCreateTopology(Geometry<T>& geo, SimulationSettings<T>& sims)
{
    std::string topoFileName;

    // overload with custom topology PPM from own dir?
    topoFileName = sims.outDir+"overloadTopo.ppm";
    if (sims.flatMode && sims.startT==0 && doesFileExist(topoFileName)) {
        pcout << messageLoading(topoFileName);
        geo.loadTopologyFromPPM(topoFileName);
        geo.saveTopology(sims.topoNameStem);
    } else {

        // is there a file in own dir?
        topoFileName = sims.outDir+sims.topoName(sims.startT);
        if (doesFileExist(topoFileName)) {
            pcout << messageLoading(topoFileName);
            geo.loadTopology(topoFileName);
        } else {

            // is there a file in parent dir?
            topoFileName = sims.topoName(sims.startT);
            if (doesFileExist(topoFileName)) {
                pcout << messageLoading(topoFileName);
                geo.loadTopology(topoFileName);
            } else {

                pcout << " create new topology..." << endl;
                geo.createTopology(
                        sims.shape, sims.inletCentre, sims.inletRadius, sims.outletCentre, sims.outletRadius,
                        sims.bufferShiftIn, sims.bufferShiftOut,
                        sims.cross, sims.freezeInnerSide, sims.openFlow,
                        sims.inletHeight, sims.outletHeight);
                topoFileName = sims.outDir+sims.topoNameStem;
                pcout << messageSaving(topoFileName);
                geo.saveTopology(topoFileName);
            }
        }
    }
    geo.updateTopology(geo.getTopology());
}

void writeTopoPPM(MultiScalarField3D<int>& topo, plint iter, plint iterSeros, bool color)
{
    TopoWriter<int> topologyExport(256);
    topologyExport.writePpm(createLongFileName("topoZ", iter, 8, "seros", iterSeros), color, topo);
    //topologyExport.writePpm(createLongFileName("topoX", iter, 8, "seros", iterSeros), color, topo, 0);
    //topologyExport.writePpm(createLongFileName("topoY", iter, 8, "seros", iterSeros), color, topo, 1);
}

void writeTopoBinary(MultiScalarField3D<int>& topo, SimulationSettings<T>& sims, plint iter)
{
    TopoWriter<int> topologyExport(256);
    topologyExport.writeBinary(sims.topoName(iter), topo);
}

void updateLattice(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice, Geometry<T>& geo)
{
    defineDynamics(NSlattice, geo.getTopology(),
            new BounceBack<T,NSDESCRIPTOR>(1.), flagBB);
    defineDynamics(NSlattice, geo.getTopology(),
            new BounceBack<T,NSDESCRIPTOR>(1.), flagConstraintBB);
    defineDynamics(NSlattice, geo.getTopology(),
            new BounceBack<T,NSDESCRIPTOR>(0.), flagSolid);
    defineDynamics(NSlattice, geo.getTopology(),
            new BounceBack<T,NSDESCRIPTOR>(0.), flagConstraintSolid);
}

void setupBoundaryConditions(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice,
        SimulationSettings<T>& sims)
{
    unique_ptr<OnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>>
        NSboundaryCondition(createLocalBoundaryCondition3D<T,NSDESCRIPTOR>());

    if (sims.freeSlip) {
        NSboundaryCondition->setVelocityConditionOnBlockBoundaries(
                NSlattice, sims.getBottomLayer(), boundary::freeslip);
        NSboundaryCondition->setVelocityConditionOnBlockBoundaries(
                NSlattice, sims.getTopLayer(), boundary::freeslip);
    }

    if (sims.inverseFlow) {
        NSboundaryCondition->setVelocityConditionOnBlockBoundaries(
                NSlattice, sims.getOutletWindow(), boundary::dirichlet);
        if (sims.shape == "T" || sims.shape == "Y")
            NSboundaryCondition->setVelocityConditionOnBlockBoundaries(
                    NSlattice, sims.getOutletTwoWindow(), boundary::dirichlet);

        if (sims.shape == "U") integrateProcessingFunctional(
                new FluidPressureOutlet3D<T,NSDESCRIPTOR,1,1>(), sims.getInletWindow(), NSlattice, 0);
        else                   integrateProcessingFunctional(
                new FluidPressureOutlet3D<T,NSDESCRIPTOR,1,1>(), sims.getInletWindow(), NSlattice, 0);
        if (sims.shape == "Ti") integrateProcessingFunctional(
                new FluidPressureOutlet3D<T,NSDESCRIPTOR,1,-1>(), sims.getInletTwoWindow(), NSlattice, 0);
    } else {
        NSboundaryCondition->setVelocityConditionOnBlockBoundaries(
                NSlattice, sims.getInletWindow(), boundary::dirichlet);
        if (sims.shape == "Ti")
            NSboundaryCondition->setVelocityConditionOnBlockBoundaries(
                    NSlattice, sims.getInletTwoWindow(), boundary::dirichlet);

        if (sims.shape == "Y") {
            integrateProcessingFunctional(
                new FluidPressureOutlet3D<T,NSDESCRIPTOR,1,-1>(), sims.getOutletWindow(), NSlattice, 0);
            integrateProcessingFunctional(
                new FluidPressureOutlet3D<T,NSDESCRIPTOR,1,-1>(), sims.getOutletTwoWindow(), NSlattice, 0);
        } else {

            if (sims.shape == "U") integrateProcessingFunctional(
                    new FluidPressureOutlet3D<T,NSDESCRIPTOR,1,1>(), sims.getOutletWindow(), NSlattice, 0);
            else                   integrateProcessingFunctional(
                    new FluidPressureOutlet3D<T,NSDESCRIPTOR,0,1>(), sims.getOutletWindow(), NSlattice, 0);
            if (sims.shape == "T")
                integrateProcessingFunctional(
                    new FluidPressureOutlet3D<T,NSDESCRIPTOR,0,-1>(), sims.getOutletTwoWindow(), NSlattice, 0);
        }
    }

    NSlattice.periodicity().toggleAll(false);
    if (sims.flatMode)
        NSlattice.periodicity().toggle(2, true);
}


bool loadOrInitPopulation(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice,
        SimulationSettings<T>& sims)
{
    std::string popuFileName;

    // is there a file in own dir?
    popuFileName = sims.outDir+sims.popuName(sims.startT);
    if (doesFileExist(popuFileName) == false) {

        // is there a file in parent dir?
        popuFileName = sims.popuName(sims.startT);
        if (doesFileExist(popuFileName) == false) {

            if (sims.startT == 0) {
                pcout << " initialize lattice population..." << endl;
                initializeAtEquilibrium(NSlattice, NSlattice.getBoundingBox(),
                        (T) 1.0, Array<T,3>((T)0, (T)0, (T)0));
                NSlattice.initialize();
                return false;
            } else {
                pcout << messageLost(popuFileName); // file is not in own directory
                pcout << " -> terminate program." << endl << endl;
                exit(1);
            }
        }
    }

    if (sims.startT == 0) pcout << " initialized fluid flow population exists." << endl;
    // Load saved data from a previous simulation
    pcout << messageLoading(popuFileName);
    loadBinaryBlock(NSlattice, popuFileName);
    // set inlet velocity:
    updateBoundaryVelocity(NSlattice, sims, sims.rampT);
    // e.g. population_seros3D_L_Re500_x50_iT194000_time38.8.dat
    pcout << "Resuming iteration after iTotal=" << sims.startT << endl;

    return true;
}

void simulateFlow(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice, plint& iter)
{
    global::timer("collideAndStream").start();
    NSlattice.collideAndStream();
    global::timer("collideAndStream").stop();
    iter++;
}

void analyseFlow(Measurements3D<T>& measure, SimulationSettings<T>& sims,
                 plint iter, bool isInitialized)
{
    global::timer("measurements").start();
    measure.takeMeasurements(iter, sims.param->stepToDimless(iter), sims.getInletRegion(), sims.trans);

    std::stringstream histFileName, distFileName;
    distFileName << "distribution" << iter;
    histFileName << "histogram" << iter;
    if (isInitialized) {
        appendToFile(measure.getDataFlow(), sims.outDir+sims.iterFileName);
            rewriteFile(measure.getDistributionData(), sims.distDir+"/"+distFileName.str());
        if (sims.doHist)
            rewriteFile(measure.getHistogramData(), sims.histDir+"/"+histFileName.str());
    } else {
        appendToFile(measure.getDataFlow(), sims.outDir+sims.iterInitializeFileName);
        rewriteFile(measure.getDistributionData(), sims.distDirInit+"/"+distFileName.str());
        if (sims.doHist)
            rewriteFile(measure.getHistogramData(), sims.histDirInit+"/"+histFileName.str());
    }
    global::timer("measurements").stop();
}

bool analyseFlowAvg(Measurements3D<T>& measure, SimulationSettings<T>& sims, bool writeFile=true)
{
    //get averaged data since last seros and reset averaging:
    measure.calculateSmoothingAndGradient();
    if (writeFile)
        appendToFile(measure.getDataSeros(), sims.outDir+sims.serosFileName);
    return measure.checkBestSeros();
}

bool loadOrCreateMeasurementsFile(Measurements3D<T>& measure, SimulationSettings<T>& sims, bool isInitialized)
{
    std::string dataFlowNewFile;
    std::string dataSerosNewFile;
    if (sims.startT == 0) {
        if (isInitialized) {
            dataFlowNewFile  = sims.outDir+sims.iterFileName;
            renameAllExisting(sims.outDir+sims.iterStemFileName);
        } else {
            dataFlowNewFile  = sims.outDir+sims.iterInitializeFileName;
        }
        dataSerosNewFile = sims.outDir+sims.serosFileName;
        renameAllExisting(sims.outDir+sims.serosStemFileName);

        pcout << " setup new measurements files..." << endl;
        rewriteFile(measure.getDataHeaderFlow(), dataFlowNewFile);
        rewriteFile(measure.getDataHeaderSeros(), dataSerosNewFile);

        pcout << " take first measurements..." << endl << endl;
        measure.takeMeasurements(sims.startT, sims.param->stepToDimless(sims.startT), sims.getInletRegion(), false);
        pcout << " append first measurements to files..." << endl << endl;
        appendToFile(measure.getDataFlow(),  dataFlowNewFile);
        analyseFlowAvg(measure, sims, false);
        measure.resetAvg();
        return false;

    } else {
        pcout << " continue from existing measurements files..." << endl << endl;
        std::string lastLine = readLastLineInFile(sims.outDir+sims.iterFileName) +
                               readLastLineInFile(sims.outDir+sims.serosFileName);
        std::vector<string> lastLineElements = splitString(lastLine);
        std::vector<T> recentData = convertStringToDouble(lastLineElements);
        measure.load(recentData);
        return true;
    }
}

void createCheckpoint(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice, Geometry<T>& geo,
                      Measurements3D<T>& measure, SimulationSettings<T>& sims,
                      bool savePopu, plint iter)
{
    if (savePopu) {
        std::string popuSaveName = sims.popuNameStem;
        pcout << messageSaving(popuSaveName);
        saveBinaryBlock(NSlattice, sims.outDir+popuSaveName);
    }
    std::string topoSaveName = sims.topoName(iter);
    pcout << messageSaving(topoSaveName);
    geo.saveTopology(sims.outDir+topoSaveName);
    pcout << endl;
}

void logger(Measurements3D<T>& measure, SimulationSettings<T>& sims,
            bool showSeros=false)
{
    if (showSeros) {
        pcout << "seros: "
              << measure.getSerosCount() << endl;
        pcout << "duration of seros iteration: "
              << global::timer("seros").getTime() << endl;
    }
    pcout << "duration of collideAndStream cycle: "
          << global::timer("collideAndStream").getTime() << endl
          << "duration of measurement cycle: "
          << global::timer("measurements").getTime() << endl;
    global::timer("collideAndStream").reset();
    global::timer("measurements").reset();
    global::timer("seros").reset();
}

double getAreaFromSTL(SimulationSettings<T>& sims, bool newBestSeros=false)
{
    double area;
    if (global::mpi().isMainProcessor()) {
        system(isoSurfaceCommand(sims, newBestSeros).c_str());
        area = numberInString(ReadNthLine(sims.isoSurface, 2));
    }
    global::mpi().bCast(&area, 1);
    return area;
}

void reshape(MultiBlockLattice3D<T, NSDESCRIPTOR>& NSlattice,
             Geometry<T>& geo, Measurements3D<T>& measure,
             SimulationSettings<T>& sims)
{
    global::timer("seros").start();
    T fractionSedi = sims.factorSeros;
    if (sims.doEros) {
        geo.erode(measure.getShearStressNormField(), measure.getErosThreshold(sims.factorSeros));
        if (sims.cSma == 0.) {
            defineDynamics(NSlattice, geo.getTopology(),
                new RegularizedBGKdynamics<T,NSDESCRIPTOR>( sims.omega ), flagEroding );
        } else {
            defineDynamics(NSlattice, geo.getTopology(),
                new SmagorinskyRegularizedDynamics<T,NSDESCRIPTOR>( sims.omega, 0.14 ), flagEroding );
        }
        fractionSedi = measure.getFractionEroded();
    }

    if (sims.serosConstraint != 'V' && sims.serosConstraint != 'D' && sims.doSedi) {
        // apply change to topo (with compensation if necessary)
        geo.sediment(geo.getTopology(), measure.getShearStressNormField(),
                     measure.getSediThreshold(fractionSedi, 0));
    }

    if (sims.serosConstraint == 'V' && sims.doSedi) {
        MultiScalarField3D<int> topoTemp = geo.getTopology();
        // calculate new geometry without applying change to topo:
        geo.sediment(topoTemp, measure.getShearStressNormField(),
                     measure.getSediThreshold(fractionSedi));
        geo.updateTopology(topoTemp, false);
        measure.calculateDiameter(topoTemp);
        // calculate unwanted change of fluid volume
        plint fluidCellsIncrease = measure.getFluidCellsIncrease();
        // apply change to topo (with compensation if necessary)
        geo.sediment(geo.getTopology(), measure.getShearStressNormField(),
                     measure.getSediThreshold(fractionSedi, fluidCellsIncrease));
    }

    if (sims.serosConstraint == 'D' && sims.doSedi) {
        MultiScalarField3D<int>& topo = geo.getTopology();
        geo.updateTopology(topo);
        T diameter, areaFromSTL=0;
        if (sims.doSTL) {
            writeTopoBinary(geo.getTopologyForSTL(), sims, -1);
            areaFromSTL = getAreaFromSTL(sims);
        }
        diameter = measure.calculateDiameter(topo, areaFromSTL);
        //pcout << " -> diameter = " << diameter << "(0)" << endl;
        T lastFractionSedi = fractionSedi;
        for(int i=1; diameter>=sims.targetDiameter && i<30; i++) {
            geo.markErodingCellsInTopology(topo);
            if (sims.flatMode)
                fractionSedi = (diameter*diameter)-1.;
            else
                fractionSedi = (diameter*diameter*diameter)-1.;
            //pcout << "fractionSedi = " << fractionSedi;
            if (fractionSedi == lastFractionSedi) break;
            measure.calculateWallStress(topo);
            geo.sediment(topo, measure.getShearStressNormField(),
                         measure.getSediThreshold(fractionSedi));
            geo.updateTopology(topo);
            if (sims.doSTL) {
                writeTopoBinary(geo.getTopologyForSTL(), sims, -1);
                areaFromSTL = getAreaFromSTL(sims);
            }
            diameter = measure.calculateDiameter(topo, areaFromSTL);
            //pcout << " -> diameter = " << diameter << "(" << i << ")" << endl;
        }
    }

    geo.updateTopology(geo.getTopology());
    double areaFromSTL = 0;
    if (sims.doSTL) {
        writeTopoBinary(geo.getTopologyForSTL(), sims, -1);
        areaFromSTL = getAreaFromSTL(sims);
    }
    measure.calculateDiameter(geo.getTopology(), areaFromSTL);

    global::timer("seros").stop();
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);

    if (!correctArguments(argc, argv)) {
        pcout << "Usage: " << argv[0] << " xmlConfigFile (opional: ReynoldsNumber resolution)" << endl;
        exit(1);
    }

    SimulationSettings<T> sims;
    if (argc == 4) {
        try { global::argv(2).read(sims.Re); }
        catch(PlbIOException& exception) {
            // Print the corresponding error message.
            pcout << exception.what() << endl;
            exit(1);
        }
        try { global::argv(3).read(sims.resolution); }
        catch(PlbIOException& exception) {
            // Print the corresponding error message.
            pcout << exception.what() << endl;
            exit(1);
        }
    } else sims.resolution = sims.Re = 0;

    string xmlConfigFile;
    try { global::argv(1).read(xmlConfigFile); }
    catch(PlbIOException& exception) {
        // Print the corresponding error message.
        pcout << exception.what() << endl;
        exit(1);
    }

    sims.read(xmlConfigFile);
    checkCreateDir(sims.outDir);
    if (sims.doHist) {
        checkCreateDir(sims.histDir);
        checkCreateDir(sims.histDirInit);
    }
    checkCreateDir(sims.distDir);
    checkCreateDir(sims.distDirInit);
    pcout << " copy xml file to directory..." << endl;
    copyFile(xmlConfigFile, sims.outDir+xmlConfigFile);
    pcout << " copy analysis script to directory..." << endl;
    copyFile("./analysis.gnu", sims.outDir+"analysis.gnu");


    pcout << " create lattice..." << endl;
    MultiBlockLattice3D<T,NSDESCRIPTOR> NSlattice(
        sims.Nx,sims.Ny,sims.Nz, new RegularizedBGKdynamics<T,NSDESCRIPTOR>(sims.omega));

    if (sims.cSma == 0.) {
        sims.printWriteLogFile("seros3D - RegularizedBGKDynamics");
    } else {
        NSlattice = MultiBlockLattice3D<T,NSDESCRIPTOR>(
            sims.Nx,sims.Ny,sims.Nz, new SmagorinskyRegularizedDynamics<T,NSDESCRIPTOR>(sims.omega, sims.cSma));
        sims.printWriteLogFile("seros3D - SmagorinskyRegularizedDynamics");
    }

    pcout << " create stationary checker..." << endl;
    belowThresholdCounter checkKineticEnergy(sims.kineticEnergyAvgSmoothGradThreshold, sims.rampT/sims.analysisT*2);
    belowThresholdCounter checkDensityDiff(sims.densityAvgDiffSmoothGradSmoothThreshold, sims.rampT/sims.analysisT*2);
    bool isInitialized = false;
    bool isOptimised = false;

    pcout << " create geometry..." << endl;
    Geometry<T> geo(NSlattice, sims.flatMode);
    pcout << " load or create topology..." << endl;
    loadOrCreateTopology(geo, sims);
    pcout << " write topology to disk as PPM..." << endl;
    writeTopoPPM(geo.getTopology(), sims.startT, 0, true);
    //writeTopoPPM(geo.getTopologyForSTL(), sims.startT, 0, false);
    pcout << "                    and as VTK..." << endl;
    writeTopoVTK(geo, sims.startT, 0);
    pcout << "                    and as binary..." << endl;
    writeTopoBinary(geo.getTopologyForSTL(), sims, sims.startT);

    pcout << " update lattice to resemble topology..." << endl;
    updateLattice(NSlattice, geo);
    pcout << " setup lattice boundary conditions..." << endl;
    setupBoundaryConditions(NSlattice, sims);
    pcout << " load or init lattice population..." << endl;
    isInitialized = loadOrInitPopulation(NSlattice, sims);

    double areaFromSTL=0;
    if (sims.doSTL) {
        writeTopoBinary(geo.getTopologyForSTL(), sims, -1);
        areaFromSTL = getAreaFromSTL(sims);
        pcout << " save initial STL" << endl;
        copyFile(sims.isoSurface, sims.isoSurface_init);
    }
    pcout << " create measurements..." << endl;
    Measurements3D<T> measure(NSlattice, geo.getTopology(), sims.resolution,
                              sims.smoothFactor, sims.flatMode, sims.doHist, areaFromSTL);

    pcout << " create ppm writing device..." << endl;
    ppmWriter writeFlowPPMs(measure, sims.rhoMax, sims.tauMax, sims.uLB*1.2, sims.trans, sims.tauScaled);

    pcout << endl;

    T ReNow;
    plint iTotal = 0;
    plint iT = 0;
    // Initialization process:
    if (isInitialized == false) {
        ReNow = 0;
        pcout << "Start with initializing fluid flow from Re=0" << endl;
        loadOrCreateMeasurementsFile(measure, sims, isInitialized);
        pcout << " write initial PPMs..." << endl;
        writeFlowPPMs(0, 0);
        pcout << " initialization process starts now..." << endl;

        do {
            iT++;
            if (ReNow < sims.Re)
                ReNow = T(iT)/sims.rampT * sims.Re;

            if (iT%sims.logT==0) {
                writeFlowPPMs(iT, 0);
                pcout << "Time " << sims.param->stepToDimless(iT) << " (step " << iT << ")"
                      << ": Current Re=" << ReNow << endl;
            }

            updateBoundaryVelocity(NSlattice, sims, iT);
            simulateFlow(NSlattice, iTotal);
            if (iT%sims.analysisT==0)
                analyseFlow(measure, sims, iTotal, isInitialized);

        } while (iT<sims.rampT*1.5 || (!sims.trans && checkKineticEnergy.repeatedlyBelowThrehold(measure.getKineticEnergyAvgSmoothGrad())));

        logger(measure, sims, false);
        pcout << "Fluid velocity field initialized!" << endl << endl;
        isInitialized = true;
        checkKineticEnergy.reset();
        iTotal = 0;

        pcout << " write VTK file of flow field..." << endl;
        writeFlowVTK(NSlattice, measure, iTotal, measure.getSerosCount());
        createCheckpoint(NSlattice, geo, measure, sims, isInitialized, iTotal);
    } else {
        ReNow = sims.Re;
    }

    pcout << "Ready for reshaping of flow geometry!" << endl;
    loadOrCreateMeasurementsFile(measure, sims, isInitialized);
    pcout << " write initial PPMs..." << endl;
    writeFlowPPMs(0, 0);
    pcout << " reshaping process starts now..." << endl;

    plint serosCount = -1;
    // Reshaping loop:
    while (serosCount<sims.maxSeros && (sims.doEros || sims.doSedi)) {

        if (serosCount >= 0) {
            reshape(NSlattice, geo, measure, sims);
            measure.increaseSerosCount();
        }
        serosCount++;
        //writeTopoPPM(geo.getTopology(), iT, serosCount, true);
        //writeTopoPPM(geo.getTopologyForSTL(), iT, serosCount, false);
        measure.resetAvg();

        /*
        if (serosCount < sims.duringRampS && serosCount%betweenRampIntervall==0) {
            pcout << "Increase Re... " << endl;

            for (int iRamp=1; ReNow<sims.ReStartSeros; iRamp++) {
            for (int iRamp=1; iRamp%betweenSerosT!=0; iRamp++) {
                ReNow = T(iRamp)/sims.rampInitT * sims.ReStartSeros;
                if (iRamp%sims.logT==0)
                    pcout << "step " << iRamp << ": Current Re=" << ReNow << endl;
                updateBoundaryVelocity(NSlattice, sims, iRamp);
                simulateFlow(NSlattice, iTotal);
                if (iRamp%sims.analysisT==0)
                    analyseFlow(measure, sims, iTotal, isInitialized);
            }

            pcout << "Reached Re=" << ReNow << endl;
        }
        if (iTotal+sims.rampBeforeSerosT == sims.rampT) pcout << endl << "Ramp finished!!" << endl << endl;
        */

        //pcout << " let flow adapt to new geometry..." << endl << endl;
        iT=0;
        do {
            iT++;
            simulateFlow(NSlattice, iTotal);
            if (iT > sims.delayAnalysisT && iT%sims.analysisT==0)
                analyseFlow(measure, sims, iTotal, isInitialized);
        } while (iT%sims.serosT!=0);
        logger(measure, sims, !isOptimised);

        //pcout << " write PPM files..." << endl;
        writeFlowPPMs(iTotal, measure.getSerosCount());
        if (serosCount%sims.vtkFlowIntervall==0) {
            pcout << " write VTK file of flow field..." << endl;
            writeFlowVTK(NSlattice, measure, iTotal, measure.getSerosCount());
        }
        if (serosCount%sims.vtkTopoIntervall==0) {
            pcout << " write VTK file of topology ..." << endl;
            writeTopoVTK(geo, iTotal, measure.getSerosCount());
        }

        bool newBestSeros = analyseFlowAvg(measure, sims);
        if (newBestSeros) {
            stringstream bestSerosCount;
            bestSerosCount << measure.getSerosCount();
            rewriteFile(bestSerosCount.str(), sims.outDir+"bestSeros");
            if (sims.doSTL) {
                writeTopoBinary(geo.getTopologyForSTL(), sims, -2);
                getAreaFromSTL(sims, true);
            }
            writeFlowVTK(NSlattice, measure, -1, -1);
            createCheckpoint(NSlattice, geo, measure, sims, false, -2);
        }

        if (sims.checkIntervall!=0 && serosCount%sims.checkIntervall==0)
            createCheckpoint(NSlattice, geo, measure, sims, false, iTotal);
    }

    createCheckpoint(NSlattice, geo, measure, sims, false, iTotal);
    return 0;
}
