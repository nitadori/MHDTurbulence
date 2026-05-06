#cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_AMD_ZEN5=ON
cmake -B build_openmp -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_ARMV9_GRACE=ON
cmake -B build_cuda -DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_HOPPER90=ON
