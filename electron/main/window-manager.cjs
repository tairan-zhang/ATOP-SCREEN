const { rendererFile } = require("./paths.cjs");

function createMainWindow({ BrowserWindow, app, preload }) {
  const window = new BrowserWindow({
    width: 1480,
    height: 960,
    minWidth: 1120,
    minHeight: 700,
    show: false,
    titleBarStyle: process.platform === "darwin" ? "hiddenInset" : "hidden",
    ...(process.platform === "darwin" ? { trafficLightPosition: { x: 16, y: 12 } } : {}),
    titleBarOverlay: process.platform === "win32" ? { color: "#111315", symbolColor: "#e6e9ef", height: 36 } : undefined,
    backgroundColor: "#101214",
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      sandbox: true,
      plugins: true,
      preload,
    },
  });
  window.once("ready-to-show", () => {
    window.show();
    window.focus();
  });
  window.webContents.setWindowOpenHandler(() => ({ action: "deny" }));
  return window;
}

async function loadRenderer({ app, window }) {
  const devUrl = process.env.ELECTRON_RENDERER_URL || "http://127.0.0.1:5173";
  if (!app.isPackaged && process.env.ELECTRON_RENDERER_URL) {
    await window.loadURL(devUrl);
    if (process.env.ATOP_OPEN_DEVTOOLS === "1") window.webContents.openDevTools({ mode: "detach" });
  } else {
    await window.loadFile(rendererFile());
  }
}

module.exports = { createMainWindow, loadRenderer };
