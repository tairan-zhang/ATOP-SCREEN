const { validateAnalysisRequest, validateTaskId } = require("./validation.cjs");

function registerAnalysisHandlers({ ipcMain, analysisProcess }) {
  ipcMain.handle("analysis:start", async (_event, payload) => {
    return analysisProcess.requestAnalysis(validateAnalysisRequest(payload));
  });
  ipcMain.handle("analysis:cancel", async (_event, taskId) => {
    return analysisProcess.cancelAnalysis(validateTaskId(taskId));
  });
}

module.exports = { registerAnalysisHandlers };
