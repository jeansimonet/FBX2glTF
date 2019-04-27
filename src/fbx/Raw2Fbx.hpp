/**
 * Copyright (c) 2019-present, Adobe Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#pragma once

#include "raw/RawModel.hpp"

bool Raw2FBX(const std::string& outputFolder, const RawModel& raw, const GltfOptions& options);
