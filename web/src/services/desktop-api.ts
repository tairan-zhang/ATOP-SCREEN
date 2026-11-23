import type {
  AnalysisProgress,
  AnalysisRequest,
} from "@atop-screen/ipc-contracts";
import type { ProjectFolderInspection } from "../types/atop";

export type { AnalysisProgress, AnalysisRequest, FileInspection, GseaResult, HealthResult, ProcessResult, VisualizationResult } from "@atop-screen/ipc-contracts";

export const desktopApi = {
  isDesktop() {
    return Boolean(window.atop?.isDesktop);
  },

  async pickScreenFile() {
    return window.atop?.files.chooseInputFile() ?? null;
  },

  async pickGmtFile() {
    return window.atop?.files.chooseGmtFile() ?? null;
  },

  async pickDirectory() {
    return window.atop?.files.chooseDirectory() ?? null;
  },

  async readAssetDataUrl(filePath: string) {
    if (!window.atop) throw new Error("Preview files are available from the Electron desktop shell.");
    return window.atop.files.readAssetDataUrl(filePath);
  },

  async readAssetBytes(filePath: string) {
    if (!window.atop) throw new Error("Preview files are available from the Electron desktop shell.");
    return window.atop.files.readAssetBytes(filePath);
  },

  async downloadAsset(filePath: string) {
    if (!window.atop) throw new Error("Figure downloads are available from the Electron desktop shell.");
    return window.atop.files.downloadAsset(filePath);
  },

  async downloadArtifacts(paths: Partial<Record<"svg" | "png" | "pdf", string>>, baseName: string, projectDirectory?: string) {
    if (!window.atop) throw new Error("Figure downloads are available from the Electron desktop shell.");
    return window.atop.files.downloadArtifacts(paths, baseName, projectDirectory);
  },

  async request<T>(payload: AnalysisRequest): Promise<T> {
    if (window.atop) return window.atop.analysis.start<T>(payload);
    if (payload.command === "health") {
      return {
        status: "browser",
        rVersion: "Open the Electron desktop shell to connect the R runtime",
        loadedModules: [],
        functions: {},
        packages: {},
      } as T;
    }
    throw new Error("The ATOP analysis runtime is available from the Electron desktop shell.");
  },

  subscribeProgress(listener: (progress: AnalysisProgress) => void) {
    return window.atop?.onProgress(listener) ?? (() => undefined);
  },

  projects: {
    async list() {
      return window.atop?.projects.list() ?? [];
    },
    async inspectFolder(folderPath: string): Promise<ProjectFolderInspection> {
      if (!window.atop) throw new Error("Projects are saved by the Electron desktop shell.");
      return window.atop.projects.inspectFolder(folderPath);
    },
    async create(payload: { name?: string; folderPath: string }) {
      if (!window.atop) throw new Error("Projects are saved by the Electron desktop shell.");
      return window.atop.projects.create(payload);
    },
    async open(projectId: string, selection?: { analysisId?: string; gseaId?: string }) {
      if (!window.atop) throw new Error("Projects are saved by the Electron desktop shell.");
      return window.atop.projects.open(projectId, selection);
    },
    async saveAnalysis(payload: Record<string, unknown>) {
      if (!window.atop) throw new Error("Projects are saved by the Electron desktop shell.");
      return window.atop.projects.saveAnalysis(payload);
    },
    async saveGsea(payload: Record<string, unknown>) {
      if (!window.atop) throw new Error("Projects are saved by the Electron desktop shell.");
      return window.atop.projects.saveGsea(payload);
    },
    async saveVisualization(payload: Record<string, unknown>) {
      if (!window.atop) throw new Error("Projects are saved by the Electron desktop shell.");
      return window.atop.projects.saveVisualization(payload);
    },
    async delete(projectId: string) {
      if (!window.atop) throw new Error("Projects are saved by the Electron desktop shell.");
      return window.atop.projects.delete(projectId);
    },
  },
};

export type { AnalysisRequest as AtopAnalysisRequest };
