# Apple silicon runtime

Copy the runtime directories from the installed app by following the [developer instructions](../../../docs/DEVELOPMENT.md). The desktop build reads the R framework and its package library from `runtime/`, the native permutation library from `analysis/lib/`, and MAGeCK RRA and Python from `mageck/`. The [runtime manifest](../../runtime-manifest.json) records their versions.
