// Copyright (C) 2015 by Benjamin Chrétien, CNRS-LIRMM.
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

#ifndef ROBOPTIM_CORE_MANIFOLD_MAP_UTIL_HH
# define ROBOPTIM_CORE_MANIFOLD_MAP_UTIL_HH

# include <memory>

# include <boost/shared_ptr.hpp>

namespace roboptim
{
  /// \brief Dummy deleter doing nothing.
  /// This can be used (with great care) to wrap references into smart
  /// pointers, and should be removed as soon as memory is properly handled
  /// in this project...
  /// \tparam U object type.
  template<typename U>
  struct NoopDeleter
  {
    inline void operator() (const U*) const {}
  };

  /// \brief Convert a Boost shared pointer to a STL shared pointer.
  ///
  /// Note: does not work with weak references to the source pointer.
  /// See: http://stackoverflow.com/a/12315035/1043187
  ///
  /// \tparam T object type.
  /// \param ptr Boost shared pointer.
  /// \return STL shared pointer.
  template <typename T>
  std::shared_ptr<T> toStd (boost::shared_ptr<T>& ptr)
  {
    return std::shared_ptr<T>
      (ptr.get(), [ptr](T*) mutable { ptr.reset(); });
  }

  /// \brief Convert a STL shared pointer to a Boost shared pointer.
  ///
  /// Note: does not work with weak references to the source pointer.
  /// See: http://stackoverflow.com/a/12315035/1043187
  ///
  /// \tparam T object type.
  /// \param ptr STL shared pointer.
  /// \return Boost shared pointer.
  template <typename T>
  boost::shared_ptr<T> toBoost (std::shared_ptr<T>& ptr)
  {
    return boost::shared_ptr<T>
      (ptr.get(), [ptr](T*) mutable { ptr.reset(); });
  }
} // end of namespace roboptim

#endif //! ROBOPTIM_CORE_MANIFOLD_MAP_UTIL_HH
