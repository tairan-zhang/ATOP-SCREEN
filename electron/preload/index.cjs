const { contextBridge, ipcRenderer, webUtils } = require("electron");

contextBridge.exposeInMainWorld("atop", {
  platform: process.platform,
  isDesktop: true,
  files: {
    chooseInputFile: () => ipcRenderer.invoke("files:choose-input-file"),
    chooseGmtFile: () => ipcRenderer.invoke("files:choose-gmt-file"),
    chooseDirectory: () => ipcRenderer.invoke("files:choose-directory"),
    readAssetDataUrl: (filePath) => ipcRenderer.invoke("files:read-asset-data-url", filePath),
    readAssetBytes: (filePath) => ipcRenderer.invoke("files:read-asset-bytes", filePath),
    downloadAsset: (filePath) => ipcRenderer.invoke("files:download-asset", filePath),
    downloadArtifacts: (assets, baseName, projectDirectory) => ipcRenderer.invoke("files:download-artifacts", assets, baseName, projectDirectory),
    pathForFile: (file) => {
      try {
        return webUtils.getPathForFile(file);
      } catch (_) {
        return null;
      }
    },
  },
  analysis: {
    start: (payload) => ipcRenderer.invoke("analysis:start", payload),
    cancel: (taskId) => ipcRenderer.invoke("analysis:cancel", taskId),
  },
  projects: {
    list: () => ipcRenderer.invoke("projects:list"),
    inspectFolder: (folderPath) => ipcRenderer.invoke("projects:inspect-folder", folderPath),
    create: (payload) => ipcRenderer.invoke("projects:create", payload),
    open: (projectId, selection) => ipcRenderer.invoke("projects:open", projectId, selection),
    saveAnalysis: (payload) => ipcRenderer.invoke("projects:save-analysis", payload),
    saveGsea: (payload) => ipcRenderer.invoke("projects:save-gsea", payload),
    saveVisualization: (payload) => ipcRenderer.invoke("projects:save-visualization", payload),
    delete: (projectId) => ipcRenderer.invoke("projects:delete", projectId),
  },
  onProgress: (listener) => {
    const wrapped = (_event, payload) => listener(payload);
    ipcRenderer.on("analysis:progress", wrapped);
    return () => ipcRenderer.removeListener("analysis:progress", wrapped);
  },
});
