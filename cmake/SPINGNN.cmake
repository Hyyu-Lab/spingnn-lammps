# Copyright (2026) Yuanmai Zheng
# SPDX-License-Identifier: GPL-2.0-only
# Spin-aware ML pair styles use the native SPIN atom/pair interfaces.
if(NOT PKG_SPIN)
  message(FATAL_ERROR "The SPINGNN package requires PKG_SPIN=ON")
endif()
