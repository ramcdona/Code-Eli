/*********************************************************************************
* Copyright (c) 2024 Rob McDonald <rob.a.mcdonald@gmail.com>
*
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the Eclipse Public License v1.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v10.html
*
* Contributors:
*    Rob McDonald
********************************************************************************/

#ifndef eli_util_clamp_hpp
#define eli_util_clamp_hpp

#include "eli/code_eli.hpp"

// Should use clamp from C++17, but this will do until Code-Eli is set up for that.

namespace eli
{
  namespace util
  {
    template<class T>
    T clamp(const T& v, const T& lo, const T& hi)
    {
      return v < lo ? lo : hi < v ? hi : v;
    }
  }
}
#endif
