const path = require("node:path");

function projectRoot() {
  // __dirname is electron/main.
  return path.resolve(__dirname, "../..");
}

function runtimePlatformDirectory() {
  if (process.platform === "darwin" && process.arch === "arm64") return "macos-arm64";
  if (process.platform === "win32" && process.arch === "x64") return "windows-x64";
  return `${process.platform}-${process.arch}`;
}

function analysisRoot(app) {
  return app.isPackaged
    ? path.join(process.resourcesPath, "atop", "analysis")
    : path.join(projectRoot(), "packages", "analysis-engine");
}

function runtimeRoot(app) {
  return app.isPackaged
    ? path.join(process.resourcesPath, "atop")
    : path.join(projectRoot(), "packaging", "runtime", runtimePlatformDirectory());
}

function workerFile(app) {
  return path.join(analysisRoot(app), "inst", "worker", "analysis_worker.R");
}

function rendererFile() {
  return path.join(projectRoot(), "dist", "index.html");
}

module.exports = {
  projectRoot,
  runtimePlatformDirectory,
  analysisRoot,
  runtimeRoot,
  workerFile,
  rendererFile,
};
