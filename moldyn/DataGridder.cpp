/*
 * DataGridder.cpp
 *
 * Copyright (C) 2009 by Universitaet Stuttgart (VISUS). 
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "DataGridder.h"
#include <climits>
#include "MultiParticleDataCall.h"
#include "ParticleGridDataCall.h"
#include "param/IntParam.h"
#include "vislib/Array.h"
#include "vislib/Cuboid.h"
#include "vislib/Pair.h"
#include "vislib/PtrArray.h"
#include "vislib/RawStorageWriter.h"

using namespace megamol::core;


/*
 * moldyn::DataGridder::DataGridder
 */
moldyn::DataGridder::DataGridder(void) : Module(),
        inDataSlot("indata", "Slot to fetch flat data"),
        outDataSlot("outdata", "Slot to publicate gridded data"),
        gridSizeXSlot("gridsizex", "The grid size in x direction"),
        gridSizeYSlot("gridsizey", "The grid size in y direction"),
        gridSizeZSlot("gridsizez", "The grid size in z direction"),
        datahash(0), frameID(UINT_MAX), gridSizeX(0), gridSizeY(0),
        gridSizeZ(0), types(), grid(), vertData(), colData() {

    this->inDataSlot.SetCompatibleCall<moldyn::MultiParticleDataCallDescription>();
    this->MakeSlotAvailable(&this->inDataSlot);

    this->outDataSlot.SetCallback(moldyn::ParticleGridDataCall::ClassName(),
        moldyn::ParticleGridDataCall::FunctionName(0), &DataGridder::getData);
    this->outDataSlot.SetCallback(moldyn::ParticleGridDataCall::ClassName(),
        moldyn::ParticleGridDataCall::FunctionName(1), &DataGridder::getExtend);
    this->MakeSlotAvailable(&this->outDataSlot);

    this->gridSizeXSlot << new param::IntParam(5, 1);
    this->MakeSlotAvailable(&this->gridSizeXSlot);

    this->gridSizeYSlot << new param::IntParam(5, 1);
    this->MakeSlotAvailable(&this->gridSizeYSlot);

    this->gridSizeZSlot << new param::IntParam(5, 1);
    this->MakeSlotAvailable(&this->gridSizeZSlot);

}


/*
 * moldyn::DataGridder::~DataGridder
 */
moldyn::DataGridder::~DataGridder(void) {
    this->Release(); // implicitly calls 'release'
}


/*
 * moldyn::DataGridder::create
 */
bool moldyn::DataGridder::create(void) {
    this->types.Clear();
    this->grid.Clear();
    this->vertData.Clear();
    this->colData.Clear();
    this->gridSizeX = this->gridSizeY = this->gridSizeZ = 0;
    return true;
}


/*
 * moldyn::DataGridder::release
 */
void moldyn::DataGridder::release(void) {
    this->types.Clear();
    this->grid.Clear();
    this->vertData.Clear();
    this->colData.Clear();
    this->gridSizeX = this->gridSizeY = this->gridSizeZ = 0;
}


/*
 * moldyn::DataGridder::doSort
 */
int moldyn::DataGridder::doSort(const vislib::Pair<SIZE_T, unsigned int>& lhs,
        const vislib::Pair<SIZE_T, unsigned int>& rhs) {
    return lhs.Second() - rhs.Second();
}


/*
 * moldyn::DataGridder::getData
 */
bool moldyn::DataGridder::getData(Call& call) {
    ParticleGridDataCall *pgdc = dynamic_cast<ParticleGridDataCall *>(&call);
    MultiParticleDataCall *mpdc = this->inDataSlot.CallAs<MultiParticleDataCall>();
    if ((pgdc == NULL) || (mpdc == NULL)) return false;
    vislib::math::Cuboid<float> bbox(0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);

    *static_cast<AbstractGetData3DCall*>(mpdc) = *pgdc;
    if ((*mpdc)(1)) {
        bbox = mpdc->AccessBoundingBoxes().ClipBox();
    }

    *static_cast<AbstractGetData3DCall*>(mpdc) = *pgdc;
    if (!(*mpdc)(0)) return false; // unable to get data

    if ((mpdc->DataHash() == 0) // data not hashable
            || (mpdc->DataHash() != this->datahash) // new data
            || this->gridSizeXSlot.IsDirty() // grid changed
            || this->gridSizeYSlot.IsDirty()
            || this->gridSizeZSlot.IsDirty()
            || (mpdc->FrameID() != this->frameID)) { // new frame
        // renew the grid

        // reset all dirty flags
        this->datahash = mpdc->DataHash();
        this->frameID = mpdc->FrameID();
        this->gridSizeXSlot.ResetDirty();
        this->gridSizeYSlot.ResetDirty();
        this->gridSizeZSlot.ResetDirty();

        // allocate new grid
        this->gridSizeX = this->gridSizeXSlot.Param<param::IntParam>()->Value();
        this->gridSizeY = this->gridSizeYSlot.Param<param::IntParam>()->Value();
        this->gridSizeZ = this->gridSizeZSlot.Param<param::IntParam>()->Value();
        SIZE_T gridSize = this->gridSizeX * this->gridSizeY * this->gridSizeZ;
        if (gridSize == 0) {
            this->grid.Clear();
            return false;
        }
        this->grid.SetCount(gridSize);

        // copy types
        unsigned int typeCnt = mpdc->GetParticleListCount();
        this->types.SetCount(typeCnt);
        for (unsigned int i = 0; i < typeCnt; i++) {
            ParticleGridDataCall::ParticleType &t = this->types[i];
            MultiParticleDataCall::Particles &p = mpdc->AccessParticles(i);
           
            t.SetColourDataType(p.GetColourDataType());
            t.SetColourMapIndexValues(p.GetMinColourIndexValue(), p.GetMaxColourIndexValue());
            t.SetGlobalColour(p.GetGlobalColour());
            t.SetGlobalRadius(p.GetGlobalRadius());
            t.SetVertexDataType(p.GetVertexDataType());
        }

        // index buffer action blub
        vislib::Array<vislib::Pair<SIZE_T, unsigned int> > indexBuffer;
        SIZE_T oaCnt = 0;
        for (unsigned int i = 0; i < typeCnt; i++) {
            oaCnt += static_cast<SIZE_T>(mpdc->AccessParticles(i).GetCount());
        }
        indexBuffer.SetCount(oaCnt);
        oaCnt = 0;
        for (unsigned int i = 0; i < typeCnt; i++) {
            MultiParticleDataCall::Particles &p = mpdc->AccessParticles(i);
            const unsigned char *vertPtr = static_cast<const unsigned char*>(p.GetVertexData());
            unsigned int vertStep = p.GetVertexDataStride();
            SIZE_T c = static_cast<SIZE_T>(p.GetCount());

            switch (p.GetVertexDataType()) {
                case MultiParticleDataCall::Particles::VERTDATA_NONE:
                    continue; // done with that type already!
                case MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZ:
                    vertStep = vislib::math::Max(vertStep, 12U);
                    break;
                case MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZR:
                    vertStep = vislib::math::Max(vertStep, 16U);
                    break;
                default:
                    continue; // there is something wrong
            }

            for (SIZE_T j = 0; j < c; j++, oaCnt++, vertPtr += vertStep) {
                const float *v = reinterpret_cast<const float*>(vertPtr);

                int x = static_cast<int>((v[0] - bbox.Left()) * static_cast<float>(this->gridSizeX) / bbox.Width());
                if (x < 0) x = 0; else if (static_cast<unsigned int>(x) >= this->gridSizeX) x = this->gridSizeX - 1;
                int y = static_cast<int>((v[1] - bbox.Bottom()) * static_cast<float>(this->gridSizeY) / bbox.Height());
                if (y < 0) y = 0; else if (static_cast<unsigned int>(y) >= this->gridSizeY) y = this->gridSizeY - 1;
                int z = static_cast<int>((v[2] - bbox.Back()) * static_cast<float>(this->gridSizeZ) / bbox.Depth());
                if (z < 0) z = 0; else if (static_cast<unsigned int>(z) >= this->gridSizeZ) z = this->gridSizeZ - 1;

                indexBuffer[oaCnt].First() = j;
                indexBuffer[oaCnt].Second() = (x + (y + z * this->gridSizeY) * this->gridSizeX) * typeCnt + i;
            }
        }
        ASSERT(indexBuffer.Count() == oaCnt);
        indexBuffer.Sort(&DataGridder::doSort);

        // data buffer setup
        this->vertData.SetCount(gridSize * typeCnt);
        this->colData.SetCount(gridSize * typeCnt);

        // copy data to grid cells
        SIZE_T start = 0, pos = 0, cnt = indexBuffer.Count();
        for (unsigned int j = 0; j < gridSize; j++) {
            for (unsigned int i = 0; i < typeCnt; i++) {
                MultiParticleDataCall::Particles &p = mpdc->AccessParticles(i);
                const unsigned char *colPtr = static_cast<const unsigned char*>(p.GetColourData());
                const unsigned char *vertPtr = static_cast<const unsigned char*>(p.GetVertexData());
                unsigned int colSize = 0;
                unsigned int colStep = p.GetColourDataStride();
                unsigned int vertSize = 0;
                unsigned int vertStep = p.GetVertexDataStride();
                unsigned int ij = j * typeCnt + i;

                switch (p.GetColourDataType()) {
                    case MultiParticleDataCall::Particles::COLDATA_NONE:
                        colSize = 0;
                        break;
                    case MultiParticleDataCall::Particles::COLDATA_FLOAT_I:
                        colSize = 4;
                        break;
                    case MultiParticleDataCall::Particles::COLDATA_FLOAT_RGB:
                        colSize = 12;
                        break;
                    case MultiParticleDataCall::Particles::COLDATA_FLOAT_RGBA:
                        colSize = 16;
                        break;
                    case MultiParticleDataCall::Particles::COLDATA_UINT8_RGB:
                        colSize = 3;
                        break;
                    case MultiParticleDataCall::Particles::COLDATA_UINT8_RGBA:
                        colSize = 4;
                        break;
                    default:
                        colSize = 0;
                        break;
                }
                if (colStep < colSize) {
                    colStep = colSize;
                }

                switch (p.GetVertexDataType()) {
                    case MultiParticleDataCall::Particles::VERTDATA_NONE:
                        continue; // done with that type already!
                    case MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZ:
                        vertSize = 12;
                        break;
                    case MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZR:
                        vertSize = 16;
                        break;
                    default:
                        continue; // there is something wrong
                }
                if (vertStep < vertSize) {
                    vertStep = vertSize;
                }

                ASSERT(start == pos);
                for (; (pos < cnt) && (indexBuffer[pos].Second() == ij); pos++) ;

                if (i == 0) {
                    this->grid[j].AllocateParticleLists(typeCnt);
                }
                ParticleGridDataCall::Particles& parts = this->grid[j].AccessParticleLists()[i];
                parts.SetCount(pos - start);
                this->vertData[ij].EnforceSize((pos - start) * vertSize);
                this->colData[ij].EnforceSize((pos - start) * colSize);
                for (SIZE_T k = 0; start < pos; k++, start++) {
                    ::memcpy(this->vertData[ij].At(k * vertSize), vertPtr + (indexBuffer[start].First() * vertStep), vertSize);
                    ::memcpy(this->colData[ij].At(k * colSize), colPtr + (indexBuffer[start].First() * colStep), colSize);
                }
                float maxRad = 0.0f;
                if (vertSize == 12) {
                    maxRad = this->types[i].GetGlobalRadius();
                } else if (vertSize == 16) {
                    for (SIZE_T p = 12; p < this->vertData[ij].GetSize(); p += 16) {
                        float r = *this->vertData[ij].AsAt<float>(p);
                        if (r > maxRad) {
                            maxRad = r;
                        }
                    }
                }
                parts.SetMaxRadius(maxRad);
                parts.SetColourData(this->colData[ij]);
                parts.SetVertexData(this->vertData[ij]);
            }
        }
        ASSERT(indexBuffer.Count() == pos);

        // calc grid bounding boxes
        for (unsigned int i = 0; i < gridSize; i++) {
#ifdef _WIN32
            float minX, minY, minZ, maxX, maxY, maxZ;
#else
            float minX = 0.0f, minY = 0.0f, minZ = 0.0f, maxX = 0.0f, maxY = 0.0f, maxZ = 0.0f;
#endif
            bool first = true;

            for (unsigned int j = 0; j < typeCnt; j++) {
                unsigned int vertSize;
                switch (this->types[j].GetVertexDataType()) {
                    case MultiParticleDataCall::Particles::VERTDATA_NONE:
                        continue;
                    case MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZ:
                        vertSize = 3;
                        break;
                    case MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZR:
                        vertSize = 4;
                        break;
                    default:
                        continue;
                }
                const float gr = this->types[j].GetGlobalRadius();
                const float *verts = static_cast<const float*>(
                    this->grid[i].AccessParticleLists()[j].GetVertexData());
                // because, I know, that stride is zero
                if (first && (this->grid[i].AccessParticleLists()[j].GetCount() > 0)) {
                    minX = maxX = verts[0];
                    minY = maxY = verts[1];
                    minZ = maxZ = verts[2];
                    // radius will be fixed in loop
                    first = false;
                }
                for (SIZE_T k = 0; k < this->grid[i].AccessParticleLists()[j].GetCount(); k++, verts += vertSize) {
                    float rad = (vertSize == 4) ? verts[3] : gr;
                    if (verts[0] - rad < minX) minX = verts[0] - rad;
                    if (verts[0] + rad > maxX) maxX = verts[0] + rad;
                    if (verts[1] - rad < minY) minY = verts[1] - rad;
                    if (verts[1] + rad > maxY) maxY = verts[1] + rad;
                    if (verts[2] - rad < minZ) minZ = verts[2] - rad;
                    if (verts[2] + rad > maxZ) maxZ = verts[2] + rad;
                }
            }

           if (first) {
                minX = minY = minZ = 0.0f;
                maxX = maxY = maxZ = vislib::math::FLOAT_EPSILON;
            }

            this->grid[i].SetBoundingBox(vislib::math::Cuboid<float>(minX, minY, minZ, maxX, maxY, maxZ));
        }

        // data grid complete
    }

    pgdc->SetDataHash(this->datahash);
    pgdc->SetFrameID(this->frameID);
    pgdc->SetUnlocker(NULL);
    pgdc->SetGridDataRef(this->gridSizeX, this->gridSizeY, this->gridSizeZ, this->grid.PeekElements());
    pgdc->SetTypeDataRef(static_cast<unsigned int>(this->types.Count()), this->types.PeekElements());

    return ((this->gridSizeX * this->gridSizeY * this->gridSizeZ) > 0);
}


/*
 * moldyn::DataGridder::getExtend
 */
bool moldyn::DataGridder::getExtend(Call& call) {
    ParticleGridDataCall *pgdc = dynamic_cast<ParticleGridDataCall *>(&call);
    if (pgdc == NULL) return false;

    MultiParticleDataCall *mpdc = this->inDataSlot.CallAs<MultiParticleDataCall>();
    if (mpdc != NULL) {
        *static_cast<AbstractGetData3DCall*>(mpdc) = *pgdc;
        if ((*mpdc)(1)) {
            *static_cast<AbstractGetData3DCall*>(pgdc) = *mpdc;
            return true;
        }
    }

    return false;
}
