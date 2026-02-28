# Verify C++ Compilation
source("R/interface/cpp_permutation_interface.R")

cat("Attempting to compile C++ engine...\n")
success <- compile_cpp_engine()

if (success) {
    cat("Compilation SUCCESS!\n")
} else {
    cat("Compilation FAILED!\n")
    quit(status = 1)
}
