reinstall_MegaLMM = function(local_library = '~/R/x86_64-pc-linux-gnu-library/3.3/')
{
  require(devtools)
  require(withr)
  MegaLMM_path = 'https://github.com/deruncie/SparseFactorMixedModel'
  withr::with_libpaths(local_library,devtools::install_git(MegaLMM_path,branch = 'develop',subdir = 'MegaLMM'))
}
