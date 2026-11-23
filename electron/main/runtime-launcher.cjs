const path = require("node:path");
const { analysisRoot, runtimeRoot } = require("./paths.cjs");

function createWorkerLaunch({ app, workerFile }) {
  if (!app.isPackaged) {
    return {
      command: process.env.ATOP_RSCRIPT || "Rscript",
      args: ["--vanilla", workerFile],
    };
  }
  const root = runtimeRoot(app);
  if (process.platform === "darwin") {
    const rHome = path.join(root, "runtime", "R.framework", "Resources");
    return {
      command: path.join(rHome, "bin", "exec", "R"),
      args: ["--no-echo", "--no-restore", "--vanilla", `--file=${workerFile}`],
    };
  }
  return {
    command: path.join(root, "runtime", "R", "bin", "Rscript.exe"),
    args: ["--vanilla", workerFile],
  };
}

function createWorkerEnvironment(app) {
  const root = runtimeRoot(app);
  const rHome = process.platform === "darwin"
    ? path.join(root, "runtime", "R.framework", "Resources")
    : path.join(root, "runtime", "R");
  const runtimeLibrary = path.join(rHome, "library");
  const cppLibrary = app.isPackaged
    ? path.join(analysisRoot(app), "lib", process.platform === "win32" ? "atop_permutation_engine.dll" : "atop_permutation_engine.so")
    : path.join(root, "analysis", "lib", process.platform === "win32" ? "atop_permutation_engine.dll" : "atop_permutation_engine.so");
  const localeEnvironment = process.platform === "darwin"
    ? { LANG: "en_US.UTF-8", LC_ALL: "en_US.UTF-8", LC_CTYPE: "en_US.UTF-8" }
    : {};
  return {
    ...process.env,
    ...localeEnvironment,
    ATOP_ANALYSIS_ROOT: analysisRoot(app),
    ATOP_RUNTIME_ROOT: root,
    ATOP_CPP_ENGINE_LIBRARY: cppLibrary,
    R_HOME: rHome,
    R_LIBS: runtimeLibrary,
    R_LIBS_USER: runtimeLibrary,
    R_SHARE_DIR: path.join(rHome, "share"),
    R_INCLUDE_DIR: path.join(rHome, "include"),
    R_DOC_DIR: path.join(rHome, "doc"),
    PATH: [path.join(rHome, "bin"), process.env.PATH || ""].filter(Boolean).join(path.delimiter),
    DYLD_LIBRARY_PATH: process.platform === "darwin"
      ? [path.join(rHome, "lib"), process.env.DYLD_LIBRARY_PATH || ""].filter(Boolean).join(path.delimiter)
      : process.env.DYLD_LIBRARY_PATH,
  };
}

module.exports = { createWorkerEnvironment, createWorkerLaunch };
