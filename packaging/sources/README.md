# Pinned sources

The source archives supply all 202 packages in renv.lock and the native MAGeCK source. SHA256SUMS records the checksum of each archive. Populate `r/` from the [release source archive](../../docs/DEVELOPMENT.md#locked-r-sources). This directory contains the CRAN and Bioconductor source distributions and the GitHub sources at their locked commits. The restore script verifies these files before installation.

- MAGeCK 0.5.9.5 is the original v0.5.9.5 source archive, SHA-256 b06a18036da63959cd7751911a46727aefe2fb1d8dd79d95043c3e3bdaf1d93a. Its accompanying BSD license is in MAGeCK-LICENSE. [Upstream source](https://bitbucket.org/liulab/mageck/get/v0.5.9.5.tar.gz)
- GseaVis 0.1.1 was rebuilt from commit 1678cfbafe4388b40d9111299e5644b9454b95d4. In R/funcs-from-others.R the build imports get_organism from GOSemSim and removes the unused DOSE build_Anno binding. DESCRIPTION pins ggsankey to b675d0d5144b1b5758d3b2b41e86ceee66a1e071. These changes reproduce the plotting implementation used by the desktop release. [Upstream commit](https://github.com/xpf10/GseaVis/commit/1678cfbafe4388b40d9111299e5644b9454b95d4)
- renv 1.2.1 bootstraps installation into a dedicated library. [CRAN archive](https://cran.r-project.org/src/contrib/Archive/renv/)
