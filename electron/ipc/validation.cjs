const COMMANDS = new Set(["health", "inspect_file", "process_data", "gsea", "visualize", "restore_state"]);

function isRecord(value) {
  return value !== null && typeof value === "object" && !Array.isArray(value);
}

function validateAnalysisRequest(payload) {
  if (!isRecord(payload)) throw new TypeError("Analysis request must be an object");
  if (typeof payload.command !== "string" || !COMMANDS.has(payload.command)) {
    throw new TypeError("Analysis request has an unsupported command");
  }
  const request = { command: payload.command };
  if (payload.path !== undefined) {
    if (typeof payload.path !== "string" || payload.path.length === 0) throw new TypeError("Analysis path must be a non-empty string");
    request.path = payload.path;
  }
  if (payload.gmtPath !== undefined) {
    if (typeof payload.gmtPath !== "string" || payload.gmtPath.length === 0) throw new TypeError("GMT path must be a non-empty string");
    request.gmtPath = payload.gmtPath;
  }
  if (payload.parameters !== undefined) {
    if (!isRecord(payload.parameters)) throw new TypeError("Analysis parameters must be an object");
    request.parameters = payload.parameters;
  }
  return request;
}

function validateTaskId(taskId) {
  if (typeof taskId !== "string" || !/^request-\d+$/.test(taskId)) throw new TypeError("Invalid analysis task id");
  return taskId;
}

module.exports = { validateAnalysisRequest, validateTaskId };
