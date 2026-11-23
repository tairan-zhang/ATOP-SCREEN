const fs = require("node:fs");
const path = require("node:path");

const assetMimeTypes = {
  ".png": "image/png",
  ".jpg": "image/jpeg",
  ".jpeg": "image/jpeg",
  ".webp": "image/webp",
  ".gif": "image/gif",
  ".svg": "image/svg+xml",
  ".pdf": "application/pdf",
};

function selectedPath(result) {
  return result.canceled ? null : result.filePaths[0] || null;
}

function isPathWithin(rootPath, targetPath) {
  const relative = path.relative(rootPath, targetPath);
  return relative === "" || (!relative.startsWith("..") && !path.isAbsolute(relative));
}

function resolveAsset(assetPath, policy = {}) {
  if (typeof assetPath !== "string" || assetPath.trim().length === 0) {
    throw new TypeError("Asset path must be a non-empty string");
  }
  const resolvedPath = path.resolve(assetPath);
  const allowedFiles = policy.allowedFiles || new Set();
  const allowedRoots = policy.allowedRoots || [];
  const allowed = allowedFiles.has(resolvedPath)
    || allowedRoots.some((rootPath) => isPathWithin(rootPath, resolvedPath));
  if (!allowed) throw new Error("Preview access was not granted for this file");
  const extension = path.extname(resolvedPath).toLowerCase();
  const mimeType = assetMimeTypes[extension];
  if (!mimeType) throw new TypeError("Only image, SVG, and PDF previews are supported");
  const info = fs.statSync(resolvedPath);
  if (!info.isFile()) throw new TypeError("Preview target is not a file");
  const maxBytes = extension === ".pdf" ? 50 * 1024 * 1024 : 25 * 1024 * 1024;
  if (info.size > maxBytes) throw new RangeError("Preview file is too large");
  return { extension, mimeType, resolvedPath };
}

function readAssetDataUrl(assetPath, policy = {}) {
  const { mimeType, resolvedPath } = resolveAsset(assetPath, policy);
  return `data:${mimeType};base64,${fs.readFileSync(resolvedPath).toString("base64")}`;
}

function readAssetBytes(assetPath, policy = {}) {
  const { mimeType, resolvedPath } = resolveAsset(assetPath, policy);
  return { mimeType, bytes: new Uint8Array(fs.readFileSync(resolvedPath)) };
}

function downloadFileName(value) {
  return String(value || "ATOP_SCREEN_figure")
    .trim()
    .replace(/[\\/:*?"<>|]+/g, "-")
    .replace(/\s+/g, "_")
    .replace(/_+/g, "_")
    .replace(/^_+|_+$/g, "")
    .slice(0, 160) || "ATOP_SCREEN_figure";
}

function registerFileHandlers({ ipcMain, dialog, getWindow, app, projectStore }) {
  const assetPolicy = () => ({
    allowedRoots: [
      path.resolve(app.getPath("userData"), "artifacts"),
      ...(projectStore?.listWorkspaceRoots() || []),
    ],
  });

  ipcMain.handle("files:choose-input-file", async () => {
    const result = await dialog.showOpenDialog(getWindow(), {
      title: "Select a CRISPR screen file",
      properties: ["openFile"],
      filters: [
        { name: "Screen data", extensions: ["csv", "txt", "tsv", "xlsx"] },
        { name: "All files", extensions: ["*"] },
      ],
    });
    return selectedPath(result);
  });

  ipcMain.handle("files:choose-gmt-file", async () => {
    const result = await dialog.showOpenDialog(getWindow(), {
      title: "Select a GMT gene-set file",
      properties: ["openFile"],
      filters: [
        { name: "Gene-set files", extensions: ["gmt"] },
        { name: "All files", extensions: ["*"] },
      ],
    });
    return selectedPath(result);
  });

  ipcMain.handle("files:choose-directory", async () => {
    const result = await dialog.showOpenDialog(getWindow(), {
      title: "Choose a project folder",
      properties: ["openDirectory", "createDirectory"],
    });
    return selectedPath(result);
  });

  ipcMain.handle("files:read-asset-data-url", async (_event, assetPath) => {
    try {
      return readAssetDataUrl(assetPath, assetPolicy());
    } catch (error) {
      throw new Error(error instanceof Error ? error.message : String(error));
    }
  });

  ipcMain.handle("files:read-asset-bytes", async (_event, assetPath) => {
    try {
      return readAssetBytes(assetPath, assetPolicy());
    } catch (error) {
      throw new Error(error instanceof Error ? error.message : String(error));
    }
  });

  ipcMain.handle("files:download-asset", async (_event, assetPath) => {
    try {
      const asset = resolveAsset(assetPath, assetPolicy());
      const result = await dialog.showSaveDialog(getWindow(), {
        title: "Download figure",
        defaultPath: path.basename(asset.resolvedPath),
        filters: [{ name: asset.extension.slice(1).toUpperCase(), extensions: [asset.extension.slice(1)] }],
      });
      if (result.canceled || !result.filePath) return null;
      const destination = path.resolve(result.filePath);
      if (destination !== asset.resolvedPath) fs.copyFileSync(asset.resolvedPath, destination);
      return destination;
    } catch (error) {
      throw new Error(error instanceof Error ? error.message : String(error));
    }
  });

  ipcMain.handle("files:download-artifacts", async (_event, assets, baseName, requestedProjectDirectory) => {
    try {
      const available = ["png", "svg", "pdf"].flatMap((format) => {
        const sourcePath = assets?.[format];
        if (!sourcePath) return [];
        const asset = resolveAsset(sourcePath, assetPolicy());
        return [{ format, asset }];
      });
      if (!available.length) throw new Error("No downloadable figure assets are available");
      const workspaceRoots = projectStore?.listWorkspaceRoots() || [];
      const projectDirectory = typeof requestedProjectDirectory === "string" && workspaceRoots.some((root) => path.resolve(root) === path.resolve(requestedProjectDirectory))
        ? path.join(path.resolve(requestedProjectDirectory), "exports")
        : null;
      const confirmation = await dialog.showMessageBox(getWindow(), {
        type: "question",
        buttons: projectDirectory ? ["Save to project", "Choose another folder", "Cancel"] : ["Choose folder", "Cancel"],
        defaultId: 0,
        cancelId: projectDirectory ? 2 : 1,
        message: projectDirectory ? "Save PNG, SVG, and PDF to this project?" : "Choose a folder for PNG, SVG, and PDF",
        detail: projectDirectory || "The figure will be saved to the folder you choose.",
      });
      let directory = confirmation.response === 0 && projectDirectory ? projectDirectory : null;
      if (!directory && !(confirmation.response === (projectDirectory ? 2 : 1))) {
        const result = await dialog.showOpenDialog(getWindow(), {
          title: "Choose a folder for PNG, SVG, and PDF",
          properties: ["openDirectory", "createDirectory"],
        });
        directory = selectedPath(result);
      }
      if (!directory) return [];
      fs.mkdirSync(directory, { recursive: true });
      const prefix = downloadFileName(baseName);
      return available.map(({ format, asset }) => {
        const destination = path.join(directory, `${prefix}_${format.toUpperCase()}.${format}`);
        fs.copyFileSync(asset.resolvedPath, destination);
        return destination;
      });
    } catch (error) {
      throw new Error(error instanceof Error ? error.message : String(error));
    }
  });
}

module.exports = { isPathWithin, readAssetBytes, readAssetDataUrl, registerFileHandlers, resolveAsset };
