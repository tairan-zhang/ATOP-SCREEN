const { execFileSync } = require("node:child_process");
const path = require("node:path");

module.exports = function afterPack(context) {
  if (context.electronPlatformName !== "darwin") return;
  const app = path.join(context.appOutDir, context.packager.appInfo.productFilename + ".app");
  execFileSync("python3", [path.join(__dirname, "prepare-macos-bundle.py"), app], { stdio: "inherit" });
};
