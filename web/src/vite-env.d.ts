/// <reference types="vite/client" />

import type { AnalysisProgress, AnalysisRequest, FileInspection, GseaResult, ProcessResult, VisualizationResult } from "@atop-screen/ipc-contracts";
import type { Artifact, HistoryEntry, Project, ProjectFolderInspection } from "./types/atop";

type AtopDesktopApi = {
  platform: string;
  isDesktop: boolean;
  files: {
    chooseInputFile: () => Promise<string | null>;
    chooseGmtFile: () => Promise<string | null>;
    chooseDirectory: () => Promise<string | null>;
    readAssetDataUrl: (filePath: string) => Promise<string>;
    readAssetBytes: (filePath: string) => Promise<{ mimeType: string; bytes: Uint8Array }>;
    downloadAsset: (filePath: string) => Promise<string | null>;
    downloadArtifacts: (paths: Partial<Record<"svg" | "png" | "pdf", string>>, baseName: string, projectDirectory?: string) => Promise<string[]>;
    pathForFile: (file: File) => string | null;
  };
  analysis: {
    start: <T = unknown>(payload: AnalysisRequest) => Promise<T>;
    cancel: (taskId: string) => Promise<boolean>;
  };
  projects: {
    list: () => Promise<Project[]>;
    inspectFolder: (folderPath: string) => Promise<ProjectFolderInspection>;
    create: (payload: { name?: string; folderPath: string }) => Promise<Project>;
    open: (projectId: string, selection?: { analysisId?: string; gseaId?: string }) => Promise<{
      project: Project;
      workspace: {
        activeAnalysisId: string | null;
        activeGseaId: string | null;
        screenPath: string;
        inspection: FileInspection | null;
        analysis: ProcessResult | null;
        gsea: GseaResult | null;
  gseaRun?: {gmtPath: string; config: Record<string, unknown>} | null;
        visualization: VisualizationResult | null;
        artifact: Artifact | null;
        history: HistoryEntry[];
      };
      restore: { analysisStateRds: string | null; gseaStateRds: string | null };
    }>;
    saveAnalysis: (payload: Record<string, unknown>) => Promise<{ project: Project; analysisId: string; analysis: ProcessResult; inspection: FileInspection }>;
    saveGsea: (payload: Record<string, unknown>) => Promise<{ project: Project; analysisId: string; gseaId: string; gsea: GseaResult }>;
    saveVisualization: (payload: Record<string, unknown>) => Promise<{ project: Project; visualization: VisualizationResult; artifact: Artifact }>;
    delete: (projectId: string) => Promise<boolean>;
  };
  onProgress: (listener: (progress: AnalysisProgress) => void) => () => void;
};

declare global {
  interface Window {
    atop?: AtopDesktopApi;
  }
}

export {};
