
#ifndef _CPPELEMENTFACTORY_H_
#define _CPPELEMENTFACTORY_H_
#include "cppelement.h"
#include "cppmarker.h"
#include "cppdrift.h"
#include "cppdipole.h"
#include "cppquadrupole.h"
#include "cppsextupole.h"
#include "cppoctupole.h"
#include "cpptuning.h"


class CppElementFactory{
public:
    // 根据元件类型创建对应的元件对象
    CppElement *CreateElement(CppElement* elem)
    {

        switch (elem->kind)
        {
        case MARKER:
            return new CppMarker(*dynamic_cast<CppMarker*>(elem));
            break;
        case DRIFT:
            return new CppDrift(*(dynamic_cast<CppDrift*>(elem)) );
            break;
        case DIPOLE:
            return new CppDipole(*dynamic_cast<CppDipole*>(elem));
            break;
        case QUADRUPOLE:
            return new CppQuadrupole(*dynamic_cast<CppQuadrupole*>(elem));
            break;
        case SEXTUPOLE:
            return new CppSextupole(*dynamic_cast<CppSextupole*>(elem));
            break;
        case OCTUPOLE:
            return new CppOctupole(*dynamic_cast<CppOctupole*>(elem));
            break;
        case TUNING:
            return new CppTuning(*dynamic_cast<CppTuning*>(elem));
            break;
        case EXACTDRIFT:
            return new CppExactDrift(*(dynamic_cast<CppExactDrift*>(elem)) );
            break;
        default:
            return nullptr;
            break;
        }
    }
};

#endif