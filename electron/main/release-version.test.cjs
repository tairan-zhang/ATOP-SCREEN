const assert = require("node:assert/strict");
const { readFileSync } = require("node:fs");
const path = require("node:path");
const test = require("node:test");

test("ATOP-SCREEN components use the desktop release version", () => {
  const root = path.resolve(__dirname, "../..");
  const version = require("../../package.json").version;
  const description = readFileSync(path.join(root, "packages/analysis-engine/DESCRIPTION"), "utf8");
  assert.equal(/^Version: (.+)$/m.exec(description)?.[1], version, "analysis engine");
  assert.equal(require("../../web/package.json").version, version, "renderer");
  assert.equal(require("../../packages/ipc-contracts/package.json").version, version, "IPC contracts");
  assert.equal(readFileSync(path.join(root, "VERSION"), "utf8").trim(), version, "VERSION");
});
