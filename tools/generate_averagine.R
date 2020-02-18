library(OrgMassSpecR)

atom_names <- c("C", "H", "N", "O", "S")
atom_numbers <- c(4.9384, 7.7583, 1.3577, 1.4773, 0.0417) / 111.1254
atom_mass <- c(12.011, 1.008, 14.007, 15.999, 32.066)

min_mz <- 100
max_mz <- 4000
by <- 10
mz_range <- seq(min_mz, max_mz, by)
lines <- c()
for (base_mz in mz_range) {
    print(base_mz)

    # Calculate averagine atoms.
    n_atoms <- round(atom_numbers * base_mz)
    rounded_mass <- sum(atom_mass * n_atoms)
    mass_diff <- round(base_mz - rounded_mass)
    if (mass_diff < 0) {
        n_atoms <- floor(atom_numbers * base_mz)
        rounded_mass <- sum(atom_mass * n_atoms)
        mass_diff <- round(base_mz - rounded_mass)
    }
    n_atoms[2] <- n_atoms[2] + mass_diff
    rounded_mass <- sum(atom_mass * n_atoms)
    n_atoms <- as.list(n_atoms)
    n_atoms <- setNames(n_atoms, atom_names)

    # Generate theoretical distribution.
    isotopes <- IsotopicDistribution(formula = n_atoms)
    perc <- isotopes$percent[isotopes$percent > 0.01]

    # Build the output string.
    perc <- paste(perc, collapse=", ")
    lines <- c(lines, paste("{", rounded_mass, ", {", perc, "}},"))
}

fileConn <- file("tools/averagine.txt")
writeLines(lines, fileConn)
close(fileConn)
