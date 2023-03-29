#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "../palabos/src/palabos3D.h"
#include "../palabos/src/palabos3D.hh"   // include full template code
#include "../palabos/externalLibraries/easyppm/easyppm.h"
#include "../palabos/externalLibraries/easyppm/easyppm.c"

#include <string>
#include <vector>

#include "../palabos/tools/tools_geometric.h"
#include "../palabos/tools/tools_geometric.hh"
#include "../palabos/tools/tools_plb.h"
#include "../palabos/tools/tools_plb.hh"


using namespace plb;
using namespace plb::descriptors;
using namespace std;


template<typename T>
class Geometry {
    public:
        Geometry(MultiBlockLattice3D<T, NSDESCRIPTOR>& lattice_,
                bool flatMode_);
        void updateTopology(MultiScalarField3D<int>& topo, bool updateDynamics = true);
        void createTopology(std::string shape, plint inletCentre, plint inletRadius,
                plint outletCentre, plint outletRadius,
                plint bufferShiftIn, plint bufferShiftOut,
                plint cross, bool freeze, bool openFlow, plint inletHeight, plint outletHeight);
        void loadTopology(std::string fileName);
        void loadTopologyFromPPM(std::string ppmName);
        void saveTopology(std::string fileName, bool binary=false);
        void markErodingCellsInTopology(MultiScalarField3D<int> &topo);
        void erode(MultiScalarField3D<T> &shearStressNormField, T erosThreshold);
        void sediment(MultiScalarField3D<int> &topo,
                      MultiScalarField3D<T> &shearStressNormField,
                      T sediThreshold);
        MultiScalarField3D<int>& getTopology();
        MultiScalarField3D<int>& getTopologyForSTL();

    private:
        MultiBlockLattice3D<T, NSDESCRIPTOR> &lattice;
        const plint Nx, Ny, Nz;
        MultiScalarField3D<int> topology;
        MultiScalarField3D<int> topologyForSTL;
        MultiScalarField3D<int> eroding;
        const bool flatMode; // true: only bottom layer of topology is optimised -> copied to all other layers
};


/* ******** CreateTopologyFunctional3D ************************************** */
class CreateTopologyFunctional3D : public plb::BoxProcessingFunctional3D_S<int>
{
public:
    CreateTopologyFunctional3D(
        std::string shape_, plint Nx, plint Ny, plint Nz,
        plint inletCentre,  plint inletRadius, plint outletCentre, plint outletRadius,
        plint bufferShiftIn_, plint bufferShiftOut_, plint cross_, bool openFlow_, plint inletHeight, plint outletHeight);
    virtual void process(Box3D domain, ScalarField3D<int>& topology);
    virtual CreateTopologyFunctional3D* clone() const {
        return new CreateTopologyFunctional3D(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
private:
    std::string shape;
    plint X1;
    plint Y1;
    plint Z1;
    plint InC;
    plint InR;
    plint InRSq;
    plint InR2;
    plint InR2Sq;
    plint OutC;
    plint OutR;
    plint OutRSq;
    plint OutR2;
    plint OutR2Sq;
    plint bufferShiftIn, bufferShiftOut;
    plint cross;
    bool openFlow;
};

/* ******** PPM2TopologyFunctional3D ************************************** */
template<typename T>
class PPM2TopologyFunctional3D : public plb::BoxProcessingFunctional3D_S<T>
{
public:
    PPM2TopologyFunctional3D(std::string PPMname_);
    virtual void process(Box3D domain, ScalarField3D<T>& topology);
    virtual PPM2TopologyFunctional3D<T>* clone() const {
        return new PPM2TopologyFunctional3D<T>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
private:
    std::string PPMname;
    PPM iPPM;
};

/* ******** UpdateTopologyFunctional3D ************************************** */
template<typename T>
class UpdateTopologyFunctional3D : public plb::BoxProcessingFunctional3D_S<T>
{
public:
    UpdateTopologyFunctional3D() {};
    virtual void process(Box3D domain, ScalarField3D<T>& topology);
    virtual UpdateTopologyFunctional3D<T>* clone() const {
        return new UpdateTopologyFunctional3D<T>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
};

/* ******** FreezeTopologyFunctional3D ************************************** */
template<typename T>
class FreezeTopologyFunctional3D : public plb::BoxProcessingFunctional3D_S<T>
{
public:
    FreezeTopologyFunctional3D() {}
    virtual void process(Box3D domain, ScalarField3D<T>& topology);
    virtual FreezeTopologyFunctional3D<T>* clone() const {
        return new FreezeTopologyFunctional3D<T>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
};

/* ******** SedimentFunctional3D ********************************************** */
template<typename T1, typename T2>
class SedimentFunctional3D : public plb::BoxProcessingFunctional3D_SS<T1, T2>
{
public:
    SedimentFunctional3D(T2 sediThreshold_);
    virtual void process( Box3D domain,
                          ScalarField3D<T1>& topology,
                          ScalarField3D<T2>& shearStressField );
    virtual SedimentFunctional3D<T1,T2>* clone() const {
        return new SedimentFunctional3D<T1,T2>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    T2 sediThreshold;
};

/* ******** ErodeFunctional3D ********************************************** */
template<typename T1, typename T2>
class ErodeFunctional3D : public plb::BoxProcessingFunctional3D_SS<T1, T2>
{
public:
    ErodeFunctional3D(T2 erosThreshold_) : erosThreshold(erosThreshold_) {}
    virtual void process( Box3D domain,
                          ScalarField3D<T1>& topology,
                          ScalarField3D<T2>& shearStressField );
    virtual ErodeFunctional3D<T1,T2>* clone() const {
        return new ErodeFunctional3D<T1,T2>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    T2 erosThreshold;
};

// inspired by src/io/imageWriter.h
namespace plb {
template<typename T>
class TopoWriter {
public:
    TopoWriter(plint valueRange_);
    void writePpm(
        std::string const& fName, bool colored,
        MultiScalarField2D<T>& field) const;
    void writePpm (
        std::string const& fName, bool colored,
        MultiScalarField3D<T>& field, int dir=2) const;
    void writeBinary(
        std::string const& fName,
        MultiScalarField3D<T>& field) const;

private:
    void writePpmImplementation (
        std::string const& fName, bool colored,
        ScalarField2D<T>& localField) const;
    void writeBinaryImplementation (
        std::string const& fName,
        ScalarField3D<T>& localField) const;
private:
    plint valueRange;
};
}


#endif
