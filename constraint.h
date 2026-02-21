/***************************************************************************
 *   Copyright (C) 2005-2024 by T. C. Raymond                              *
 *   tcraymond@inductivereasoning.com                                      *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  *
 *                                                                         *
 ***************************************************************************/

#pragma once

 /// This class represents a generic constraint. 

#include "Mesh/mesh.h"
//#include "field.h"
#include "fe_space.hpp"

template <class T>
class Constraint {

public:

	Constraint(CellView* Element_in, FESpaceBase<T>* fe_space_in) : entity(Element_in), fe_space(fe_space_in) { };

	virtual ~Constraint() = default;

	virtual void apply() = 0;

protected:

	CellView* entity;
	FESpaceBase<T>* fe_space;

private:

};
