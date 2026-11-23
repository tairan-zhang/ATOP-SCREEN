const assert = require("node:assert/strict");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const test = require("node:test");
const { attachArtifactDirectory } = require("./artifact-store.cjs");

test("only artifact-producing commands receive an isolated output directory", () => {
  const userData = fs.mkdtempSync(path.join(os.tmpdir(), "atop-artifacts-"));
  const app = { getPath: () => userData };
  try {
    const health = attachArtifactDirectory(app, "request-1", { command: "health" });
    assert.equal(health.artifactDir, undefined);
    const plot = attachArtifactDirectory(app, "request-2", { command: "visualize" });
    assert.equal(plot.artifactDir, path.join(userData, "artifacts", "request-2"));
    assert.equal(fs.statSync(plot.artifactDir).isDirectory(), true);
  } finally {
    fs.rmSync(userData, { recursive: true, force: true });
  }
});
