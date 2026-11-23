export type AnalysisCommand = "health" | "inspect_file" | "process_data" | "gsea" | "visualize" | "restore_state";

export type AnalysisRequest = {
  command: AnalysisCommand;
  path?: string;
  gmtPath?: string;
  parameters?: Record<string, unknown>;
};

export type AnalysisProgress = {
  event?: "progress" | "ready";
  message?: string;
};

export type AnalysisTaskId = string;

export type FileInspection = {
  path: string;
  fileName: string;
  columns: string[];
  preview: Record<string, unknown>[];
  previewRows: number;
  fileSize: number;
};

export type ProcessResult = {
  status: string;
  summary: {
    processedSgrnaRows: number;
    geneRows: number;
    filteredRows: number;
    parameters: Record<string, unknown>;
  };
  tables: {
    mageckGenes?: Record<string, unknown>[] | null;
    mageckGuides?: Record<string, unknown>[] | null;
    ncsExcluded?: Record<string, unknown>[] | null;
    ncsQc?: Record<string, unknown>[] | null;
    sgrna: Record<string, unknown>[];
    genes: Record<string, unknown>[];
    filtered: Record<string, unknown>[] | null;
  };
  artifacts?: AnalysisArtifacts;
};

export type GseaResult = {
  status: string;
  positive: Record<string, unknown>[];
  negative: Record<string, unknown>[];
  artifacts?: AnalysisArtifacts;
};

export type VisualizationResult = {
  status: string;
  plotType: string;
  title: string;
  assets: {
    svg?: string;
    png?: string;
    pdf?: string;
  };
};

export type AnalysisArtifacts = {
  directory?: string;
  normalizedCountsRds?: string;
  analysisStateRds?: string;
  gseaPositiveRds?: string;
  gseaNegativeRds?: string;
  gseaStateRds?: string;
};

export type HealthResult = {
  status: string;
  rVersion: string;
  loadedModules: string[];
  functions: Record<string, boolean>;
  packages: Record<string, boolean>;
};
