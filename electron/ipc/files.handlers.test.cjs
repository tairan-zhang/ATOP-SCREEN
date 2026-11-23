const assert = require("node:assert/strict");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const test = require("node:test");
const { isPathWithin, readAssetBytes, readAssetDataUrl } = require("./files.handlers.cjs");

test("isPathWithin rejects siblings with a shared prefix", () => {
  assert.equal(isPathWithin("/tmp/atop/artifacts", "/tmp/atop/artifacts/plot.svg"), true);
  assert.equal(isPathWithin("/tmp/atop/artifacts", "/tmp/atop/artifacts-other/plot.svg"), false);
});

test("readAssetBytes returns a typed payload suitable for Blob previews", () => {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "atop-pdf-preview-"));
  const asset = path.join(directory, "plot.pdf");
  const content = Buffer.from("%PDF-1.4\n%%EOF\n");
  fs.writeFileSync(asset, content);
  try {
    const result = readAssetBytes(asset, { allowedRoots: [directory] });
    assert.equal(result.mimeType, "application/pdf");
    assert.deepEqual(Buffer.from(result.bytes), content);
  } finally {
    fs.rmSync(directory, { recursive: true, force: true });
  }
});

test("readAssetDataUrl requires an explicit file or root grant", () => {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "atop-file-policy-"));
  const asset = path.join(directory, "plot.svg");
  fs.writeFileSync(asset, "<svg xmlns=\"http://www.w3.org/2000/svg\"/>");
  try {
    assert.throws(() => readAssetDataUrl(asset), /not granted/);
    const result = readAssetDataUrl(asset, { allowedRoots: [directory] });
    assert.match(result, /^data:image\/svg\+xml;base64,/);
  } finally {
    fs.rmSync(directory, { recursive: true, force: true });
  }
});
