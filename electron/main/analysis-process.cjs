const childProcess = require("node:child_process");
const fs = require("node:fs");
const readline = require("node:readline");
const { analysisRoot, workerFile } = require("./paths.cjs");
const { attachArtifactDirectory } = require("./artifact-store.cjs");
const { createWorkerEnvironment, createWorkerLaunch } = require("./runtime-launcher.cjs");

function createAnalysisProcess({ app, getWindow }) {
  let analysisWorker;
  let workerReader;
  let workerQueue = [];
  const pendingRequests = new Map();
  let workerReady = false;
  let requestCounter = 0;

  function flushWorkerQueue() {
    if (!workerReady || !analysisWorker || !analysisWorker.stdin.writable) return;
    while (workerQueue.length > 0) {
      const item = workerQueue.shift();
      pendingRequests.set(item.id, item);
      analysisWorker.stdin.write(`${JSON.stringify(item.request)}\n`);
    }
  }

  function startWorker() {
    if (analysisWorker && !analysisWorker.killed) return;
    const file = workerFile(app);
    if (!fs.existsSync(file)) throw new Error(`Analysis worker not found: ${file}`);
    const launch = createWorkerLaunch({ app, workerFile: file });
    analysisWorker = childProcess.spawn(launch.command, launch.args, {
      cwd: analysisRoot(app),
      env: createWorkerEnvironment(app),
      stdio: ["pipe", "pipe", "pipe"],
    });
    analysisWorker.stdin.setDefaultEncoding("utf8");
    workerReady = false;
    workerReader = readline.createInterface({ input: analysisWorker.stdout });
    workerReader.on("line", (line) => {
      let message;
      try {
        message = JSON.parse(line);
      } catch (_) {
        const window = getWindow();
        if (window && !window.isDestroyed()) window.webContents.send("analysis:progress", { message: line });
        return;
      }
      if (message.event === "ready") {
        workerReady = true;
        flushWorkerQueue();
        return;
      }
      if (message.event === "progress") {
        const window = getWindow();
        if (window && !window.isDestroyed()) window.webContents.send("analysis:progress", message);
        return;
      }
      const pending = pendingRequests.get(message.id);
      if (!pending) return;
      pendingRequests.delete(message.id);
      message.ok ? pending.resolve(message.result) : pending.reject(new Error(message.error || "Analysis request failed"));
    });
    readline.createInterface({ input: analysisWorker.stderr }).on("line", (line) => {
      const window = getWindow();
      if (!window || window.isDestroyed() || !line.trim()) return;
      let message = line.trim();
      try {
        const event = JSON.parse(line);
        if (event.event === "progress" && typeof event.message === "string") message = event.message;
      } catch { /* R package diagnostics are plain text. */ }
      window.webContents.send("analysis:progress", { message });
    });
    analysisWorker.once("exit", (code) => {
      workerReady = false;
      const error = new Error(`Analysis worker exited with code ${code}`);
      for (const item of workerQueue.splice(0)) item.reject(error);
      for (const item of pendingRequests.values()) item.reject(error);
      pendingRequests.clear();
      analysisWorker = undefined;
      workerReader?.close();
      workerReader = undefined;
    });
  }

  function requestAnalysis(request) {
    startWorker();
    const id = `request-${++requestCounter}`;
    const enrichedRequest = attachArtifactDirectory(app, id, request);
    return new Promise((resolve, reject) => {
      workerQueue.push({ id, request: enrichedRequest, resolve, reject });
      flushWorkerQueue();
    });
  }

  function cancelAnalysis(taskId) {
    const pending = pendingRequests.get(taskId);
    const queuedIndex = workerQueue.findIndex((item) => item.id === taskId);
    if (queuedIndex >= 0) {
      const [item] = workerQueue.splice(queuedIndex, 1);
      item.reject(new Error("Analysis request cancelled"));
      return true;
    }
    if (!pending) return false;
    pendingRequests.delete(taskId);
    pending.reject(new Error("Analysis request cancelled"));
    if (analysisWorker && !analysisWorker.killed) analysisWorker.kill("SIGTERM");
    return true;
  }

  function start() {
    startWorker();
  }

  function stop() {
    if (analysisWorker && !analysisWorker.killed) analysisWorker.kill("SIGTERM");
    workerReader?.close();
    analysisWorker = undefined;
    workerReady = false;
  }

  return { start, requestAnalysis, cancelAnalysis, stop };
}

module.exports = { createAnalysisProcess };
