# Copyright 2020 Daniel Runcie
# Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
# that can be found in the LICENSE file and available at
# https://polyformproject.org/licenses/noncommercial/1.0.0/


reinstall_MegaLMM = function(local_library = '~/R/x86_64-pc-linux-gnu-library/3.3/')
{
  require(devtools)
  require(withr)
  MegaLMM_path = 'https://github.com/deruncie/SparseFactorMixedModel'
  withr::with_libpaths(local_library,devtools::install_git(MegaLMM_path,branch = 'develop',subdir = 'MegaLMM'))
}
