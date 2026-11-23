const assert = require("node:assert/strict");
const test = require("node:test");
const { createWorkerEnvironment } = require("./runtime-launcher.cjs");

test("macOS analysis workers use a UTF-8 locale for non-ASCII paths", { skip: process.platform !== "darwin" }, () => {
  const environment = createWorkerEnvironment({ isPackaged: false });
  assert.equal(environment.LANG, "en_US.UTF-8");
  assert.equal(environment.LC_ALL, "en_US.UTF-8");
  assert.equal(environment.LC_CTYPE, "en_US.UTF-8");
});
