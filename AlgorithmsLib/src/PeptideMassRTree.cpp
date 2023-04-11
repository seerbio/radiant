//
// Created by anichols on 4/10/23.
//

#include "PeptideMassRTree.h"

#include "BiophysicalCalcs.h"
#include "ErrorUtils.h"
#include "PeptidesLibraryTron.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include <boost/geometry.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>
#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


class Q_DECL_HIDDEN PeptideMassRTree::Private
{
    using rTreeCoor = bg::model::point<double, 1, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, int> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(
            const QVector<PeptideStringWithMods> &peptideStrings,
            const AminoAcids &aminoAcids
            );

private:

    AminoAcids m_aminoAcids;
    RTree *m_rTree;

    int m_maxElements;

    QVector<PeptideStringWithMods> m_peptideStringWithMods;

};

PeptideMassRTree::Private::Private()
: m_rTree(nullptr)
, m_maxElements(16) {}

PeptideMassRTree::Private::~Private() {
    delete m_rTree;
}

namespace {

    Err decodePeptideStringWithMods(
            const QVector<PeptideStringWithMods> &peptideStringsWithMods,
            QVector<QPair<PeptideString,QHash<ResidueIndex, ModificationMass>>> *output
            ) {

        ERR_INIT

        for (const PeptideStringWithMods &pswm : peptideStringsWithMods) {

            QString peptideSequence;
            QString modString;
            e = Peptide::decodePeptideStringWithMods(
                    pswm,
                    &peptideSequence,
                    &modString
                    ); ree;

            QHash<ResidueIndex, ModificationMass> mods;
            e = Peptide::modificationFromString(
                    modString,
                    &mods
                    ); ree;

            output->push_back({peptideSequence, mods});
        }

        ERR_RETURN
    }

}//namespace
Err PeptideMassRTree::Private::init(
        const QVector<PeptideStringWithMods> &peptideStringsWithMods,
        const AminoAcids &aminoAcids
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(peptideStringsWithMods); ree;
    m_aminoAcids = aminoAcids;
    m_peptideStringWithMods = peptideStringsWithMods;

    QVector<QPair<PeptideString, QHash<ResidueIndex, ModificationMass>>> decodedPeptideStrings;
    e = decodePeptideStringWithMods(
            m_peptideStringWithMods,
            &decodedPeptideStrings
            ); ree;

    const QVector<QPair<PeptideString, double>> peptideMasses = BiophysicalCalcs::calculatePeptideMasses(
            decodedPeptideStrings,
            m_aminoAcids
            );

    e = ErrorUtils::isEqual(peptideMasses.size(), m_peptideStringWithMods.size()); ree;

    std::vector<rTreePoint> cloudLoader;
    for (int i = 0; i < m_peptideStringWithMods.size(); i++) {

        const double mass = peptideMasses.at(i).second;

        rTreeCoor coor(mass);
        cloudLoader.emplace_back(coor, i);
    }

    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(m_maxElements));

    ERR_RETURN
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

PeptideMassRTree::PeptideMassRTree() : d_ptr(new Private()) {}

PeptideMassRTree::~PeptideMassRTree() {}

Err PeptideMassRTree::init(
        const QVector<PeptideStringWithMods> &peptideStringsWithMods,
        const AminoAcids &aminoAcids
        ) {
    ERR_INIT
    e = d_ptr->init(
            peptideStringsWithMods,
            aminoAcids
            ); ree;
    ERR_RETURN
}
