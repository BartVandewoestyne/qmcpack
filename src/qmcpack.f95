! Generic module, to allow a user to only 'use qmcpack'.

module qmcpack

  use cranley_patterson
  use f_chi
  use f_gajda
  use f_genz
  use f_mathe_wei
  use f_nuyens
  use f_ronald
  use f_simple
  use f_testfunctions
  use f_timps
  use f_warnock
  use f_effective_dimension
  use f_unbounded_variation
  use f_vanbeylen
  use f_worstcase
  use mod_asian_option 
  use mod_adaptive_cube_qmc_rd
  use mod_assert
  use mod_bernoulli
  use mod_bootstrap_filter
  use mod_braaten_weller
  use mod_cli
  use mod_constants
  use mod_debug
  use mod_discrepancy
  use mod_distance
  use mod_diaphony
  use mod_extensible_lattice
  use mod_faure
  use mod_filewriter
  use mod_file_utils
  use mod_function
  use mod_haber
  use mod_haselgrove
  use mod_haselgrove_niederreiter
  use mod_halton
  use mod_hammersley
  use mod_integration
  use mod_kaino
  use mod_lattice_criteria
  use mod_lhs
  use mod_mathewei
  use mod_measurement_model
  use mod_monte_carlo
  use mod_niederreiter
  use mod_niederreiter_xing
  use mod_number_theory
  use mod_scrambling
  use mod_pdf
  use mod_periodize
  use mod_permutation
  use mod_primes
  use mod_primpoly
  use mod_radical_inverse
  use mod_rank1_lattice
  use mod_sobol
  use mod_sort
  use mod_special_functions
  use mod_sphere
  use mod_statistics
  use mod_stringutil
  use mod_sugihara_murota
  use mod_system_model
  use mod_testing
  use mod_transform
  use mod_utilities
  use mod_weighted_weyl
  use mod_weightfunction
  use mod_weyl
  use numeric_kinds

  public

end module qmcpack
