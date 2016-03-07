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
#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MERGER_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MERGER_HH

# include <map>
# include <vector>
# include <iostream>

# include <roboptim/core/manifold/config.hh>

# include <manifolds/Manifold.h>
# include <manifolds/CartesianProduct.h>

namespace roboptim
{
  /// \addtogroup roboptim_manifolds
  /// @{

  class ROBOPTIM_CORE_MANIFOLD_DLLAPI ManifoldMerger
  {
  public:
    ManifoldMerger();

    void addManifold(const mnf::Manifold& instanceManifold);

    bool contains(const mnf::Manifold& instanceManifold) const;

    void addManifolds(const ManifoldMerger& other);

    void clear();

    // TODO: deprecate this getter
    const mnf::CartesianProduct* getManifold() const;

    std::shared_ptr<const mnf::CartesianProduct> mergedManifold() const;

  private:
    std::map<long, const mnf::Manifold*> elementaryInstanceManifolds_;
    std::shared_ptr<mnf::CartesianProduct> mergedManifold_;
  };

  /// @}

} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_DECORATOR_MANIFOLD_MERGER_HH
