const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");

const artifactCommands = new Set(["process_data", "gsea", "visualize"]);

function createArtifactDirectory(app, taskId) {
  const userData = typeof app.getPath === "function"
    ? app.getPath("userData")
    : path.join(os.tmpdir(), "atop-screen");
  const directory = path.join(userData, "artifacts", taskId);
  fs.mkdirSync(directory, { recursive: true });
  return directory;
}

function attachArtifactDirectory(app, taskId, request) {
  if (!artifactCommands.has(request.command)) return { ...request, id: taskId };
  return {
    ...request,
    id: taskId,
    artifactDir: createArtifactDirectory(app, taskId),
  };
}

module.exports = { attachArtifactDirectory, createArtifactDirectory };
