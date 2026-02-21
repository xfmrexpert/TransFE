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

#include <list>
#include "Mesh/mesh.h"
#include "form_integrator.hpp"
#include "fe_space.hpp"

namespace TFEM
{
	template <typename T>
	class PrincipalForm
	{
	public:
		PrincipalForm(FiniteElementSpace* fe_space) : fe_space(fe_space) {}

		void AddDomainIntegrator(std::unique_ptr<FormIntegrator<T>> integrator)
		{
			integrators.push_back(std::move(integrator));
		}

		void Assemble()
		{
			for (const auto& cell : fe_space->getMesh()->Cells())
			{
				FiniteElementBase fe_space->getFiniteElement();
				for (auto& integrator : integrators)
				{
					integrator->evaluate(*entity);
				}
			}
		}

	protected:
		FiniteElementSpace* fe_space;
		std::list<std::unique_ptr<FormIntegrator<T>> integrators;

	};
}