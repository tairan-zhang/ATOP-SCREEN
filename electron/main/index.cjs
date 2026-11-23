const { app, BrowserWindow, dialog, ipcMain } = require("electron");
const path = require("node:path");
const { createAnalysisProcess } = require("./analysis-process.cjs");
const { createMainWindow, loadRenderer } = require("./window-manager.cjs");
const { registerAnalysisHandlers } = require("../ipc/analysis.handlers.cjs");
const { registerFileHandlers } = require("../ipc/files.handlers.cjs");
const { registerProjectHandlers } = require("../ipc/projects.handlers.cjs");
const { createProjectStore } = require("./project-store.cjs");

let mainWindow;
const analysisProcess = createAnalysisProcess({ app, getWindow: () => mainWindow });

app.whenReady().then(async () => {
  app.setName("ATOP SCREEN");
  mainWindow = createMainWindow({
    BrowserWindow,
    app,
    preload: path.join(__dirname, "..", "preload", "index.cjs"),
  });
  const projectStore = createProjectStore({ app });
  registerFileHandlers({ ipcMain, dialog, getWindow: () => mainWindow, app, projectStore });
  registerProjectHandlers({ ipcMain, projectStore });
  registerAnalysisHandlers({ ipcMain, analysisProcess });
  mainWindow.on("closed", () => {
    mainWindow = undefined;
    app.quit();
  });
  try {
    analysisProcess.start();
    await loadRenderer({ app, window: mainWindow });
  } catch (error) {
    dialog.showErrorBox("ATOP-SCREEN startup error", error instanceof Error ? error.message : String(error));
    app.quit();
  }
});

app.on("window-all-closed", () => app.quit());

app.on("before-quit", () => analysisProcess.stop());
