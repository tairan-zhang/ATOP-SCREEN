function requireObject(value, label) {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw new TypeError(`${label} must be an object`);
  return value;
}

function requireProjectId(value) {
  if (typeof value !== "string" || !/^[0-9a-f-]{20,}$/i.test(value)) throw new TypeError("Invalid project id");
  return value;
}

function registerProjectHandlers({ ipcMain, projectStore }) {
  ipcMain.handle("projects:list", () => projectStore.listProjects());
  ipcMain.handle("projects:inspect-folder", (_event, folderPath) => projectStore.inspectProjectFolder(folderPath));
  ipcMain.handle("projects:create", (_event, payload) => projectStore.createProject(requireObject(payload, "Project")));
  ipcMain.handle("projects:open", (_event, projectId, selection) => {
    const safeSelection = selection && typeof selection === "object" && !Array.isArray(selection) ? selection : {};
    return projectStore.openProject(requireProjectId(projectId), safeSelection);
  });
  ipcMain.handle("projects:save-analysis", (_event, payload) => {
    requireObject(payload, "Analysis snapshot");
    requireProjectId(payload.projectId);
    return projectStore.saveAnalysis(payload);
  });
  ipcMain.handle("projects:save-gsea", (_event, payload) => {
    requireObject(payload, "GSEA snapshot");
    requireProjectId(payload.projectId);
    return projectStore.saveGsea(payload);
  });
  ipcMain.handle("projects:save-visualization", (_event, payload) => {
    requireObject(payload, "Visualization snapshot");
    requireProjectId(payload.projectId);
    return projectStore.saveVisualization(payload);
  });
  ipcMain.handle("projects:delete", (_event, projectId) => projectStore.deleteProject(requireProjectId(projectId)));
}

module.exports = { registerProjectHandlers, requireProjectId };
