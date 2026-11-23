const assert = require("node:assert/strict");
const test = require("node:test");
const { validateAnalysisRequest, validateTaskId } = require("./validation.cjs");

test("analysis request validation only retains contract fields", () => {
  assert.deepEqual(
    validateAnalysisRequest({ command: "visualize", parameters: { plotType: "gene_volcano" }, ignored: true }),
    { command: "visualize", parameters: { plotType: "gene_volcano" } },
  );
});

test("analysis validation rejects unknown commands and task ids", () => {
  assert.throws(() => validateAnalysisRequest({ command: "unknown" }), /unsupported command/);
  assert.throws(() => validateTaskId("../../worker"), /Invalid analysis task id/);
  assert.equal(validateTaskId("request-42"), "request-42");
});
