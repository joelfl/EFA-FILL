/*******************************************************************************
    Copyright (C) 2018 Kevin Sahr

    This file is part of DGGRID.

    DGGRID is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DGGRID is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/
////////////////////////////////////////////////////////////////////////////////
//
// DgHexC1Grid2D.h: DgHexC1Grid2D class definitions
//
////////////////////////////////////////////////////////////////////////////////

#ifndef DGHEXC1GRID2D_H 
#define DGHEXC1GRID2D_H

#include <cmath>

#include "DgDiscRF2D.h"
#include "DgDVec2D.h"
#include "DgIVec2D.h"

class DgPolygon;

////////////////////////////////////////////////////////////////////////////////
class DgHexC1Grid2D : public DgDiscRF2D {

   public:

      DgHexC1Grid2D (DgRFNetwork& networkIn, 
                     const DgRF<DgDVec2D, long double>& ccFrameIn,
                     const string& nameIn = "HexC12D")
         : DgDiscRF2D (networkIn, ccFrameIn, nameIn, 1.0 / sqrt(3.0), 
                       1.0 / sqrt(3.0), sqrt(3.0) / 2.0, 1.0) 
           { area_ = c(); }

      DgHexC1Grid2D (const DgHexC1Grid2D& grd) : DgDiscRF2D (grd) {}

      DgHexC1Grid2D& operator= (const DgHexC1Grid2D& grd)
           { DgDiscRF2D::operator=(grd); return *this; }

      virtual long long int dist (const DgIVec2D& add1, const DgIVec2D& add2) const
           { DgIVec2D diff = add2 - add1;
             if ((diff.i() >= 0 && diff.j() <= 0) ||
                 (diff.i() <= 0 && diff.j() >= 0))   /* different signs */
             {
                return abs(diff.i()) + abs(diff.j());
             }
             else /* same signs */
             {
                diff.setI(abs(diff.i()));
                diff.setJ(abs(diff.j()));
                return (diff.i() > diff.j()) ? diff.i() : diff.j();
             }}

   protected:

      static const long double sin60_;

      virtual void setAddVertices (const DgIVec2D& add, DgPolygon& vec) const;
      virtual void setAddNeighbors 
                                  (const DgIVec2D& add, DgLocVector& vec) const;

      virtual DgIVec2D quantify (const DgDVec2D& point) const;

      virtual DgDVec2D invQuantify (const DgIVec2D& add) const
               { return DgDVec2D(add.i() - 0.5 * add.j(), 
                                 add.j() * 1.5 * r()); }

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#endif
