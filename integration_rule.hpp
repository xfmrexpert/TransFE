/***************************************************************************
 *   Copyright (C) 2005-2025 by T. C. Raymond                              *
 *   tcraymond@inductivereasoning.com                                      *
 *                                                                         *
 *   Use of this source code is governed by an MIT-style                   *
 *   license that can be found in the LICENSE.txt file or at               *
 *   https://opensource.org/licenses/MIT.                                  *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  *
 *                                                                         *
 ***************************************************************************/

#pragma once
#include <vector>
#include "Mesh/point.h"
#include "typedefs.h"

namespace TFEM
{
	class IntegrationRule {
	public:
		virtual const std::vector<Point>& IntPts() const = 0;
		virtual const Vector<double>& Weights() const = 0;
		virtual int numIntPts() const = 0;

	};
}