// Copyright (C) 2015 by Félix Darricau, AIST, CNRS, EPITA
//                       Grégoire Duchemin, AIST, CNRS, EPITA
//
// This file is part of the roboptim.
//
// roboptim is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim.  If not, see <http://www.gnu.org/licenses/>.

#include <roboptim/core/manifold-map/decorator/manifold-merger.hh>

namespace roboptim
{
  ManifoldMerger::ManifoldMerger()
    : mergedManifold_(new mnf::CartesianProduct())
  {}

  void ManifoldMerger::addManifold(const mnf::Manifold& instanceManifold)
  {
    std::function<void(const mnf::Manifold&)> addElementaries =
      [this, &addElementaries]
      (const mnf::Manifold& manifold)
      {
	if (manifold.isElementary())
	  {
	    if(!contains(manifold))
	      {
		std::cout << "ADDED " << manifold.name() << std::endl;
	      }

	    this->elementaryInstanceManifolds_[manifold.getInstanceId()] = &manifold;

	  }
	else
	  {
	    for (size_t i = 0; i < manifold.numberOfSubmanifolds(); ++i)
	      {
		addElementaries(manifold(i));
	      }
	  }
      };

    addElementaries(instanceManifold);
  }

  bool ManifoldMerger::contains(const mnf::Manifold& instanceManifold)
  {
    return elementaryInstanceManifolds_.find(instanceManifold.getInstanceId()) != elementaryInstanceManifolds_.end();
  }

  void ManifoldMerger::addManifolds(const ManifoldMerger& other)
  {
    for (auto& ite : other.elementaryInstanceManifolds_)
      {
	if (ite.second == nullptr)
	  {
	    continue;
	  }

	addManifold(*(ite.second));
      }
  }

  void ManifoldMerger::clear()
  {
    mergedManifold_.reset(new mnf::CartesianProduct());
    elementaryInstanceManifolds_.clear();
  }

  mnf::CartesianProduct* ManifoldMerger::getManifold()
  {
    for (auto ite = this->elementaryInstanceManifolds_.begin();
	 ite != elementaryInstanceManifolds_.end();
	 ++ite)
      {
	mergedManifold_->multiply(*(ite->second));
      }

    if (mergedManifold_->representationDim() <= 0)
      {
	throw std::runtime_error("The problem should not be empty.");
      }

    return &(*mergedManifold_);
  }

}
