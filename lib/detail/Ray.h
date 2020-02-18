/*
 * \file Header for Ray class
 */

#pragma once

#include "Vector.h"

namespace nut {

/**
 * \class Ray: a point and a direction
 */
struct Ray {

  Vector const & position() const {
    return position_;
  }

  Vector const & direction() const {
    return direction_;
  }

  // private:
  Vector position_;
  Vector direction_;
}; // struct Ray

/**
 * \class Ray: a point and a direction
 */
struct Ray1 {

  Vector1 const & position() const {
    return position_;
  }

  Vector1 const & direction() const {
    return direction_;
  }

  // private:
  Vector1 position_;
  Vector1 direction_;
}; // struct Ray

} // namespace nut

// End of file
