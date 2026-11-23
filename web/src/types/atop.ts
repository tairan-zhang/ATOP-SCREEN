import type {
  FileInspection,
  GseaResult,
  ProcessResult,
  VisualizationResult,
} from "@atop-screen/ipc-contracts";

export type { FileInspection, GseaResult, HealthResult, ProcessResult, VisualizationResult } from "@atop-screen/ipc-contracts";

export type View = "home" | "analysis" | "gsea" | "gseaResults" | "visualize" | "results" | "settings";
export type Theme = "dark" | "light";
export type PlotSection = "gene" | "sgrna" | "pathway";

export type ProjectGseaRun = {
  id: string;
  name: string;
  completedAt: string | null;
  pathwayRows: number;
};

export type ProjectAnalysis = {
  id: string;
  name: string;
  completedAt: string | null;
  geneRows: number;
  gseaRuns: ProjectGseaRun[];
};

export type Project = {
  id: string;
  name: string;
  rootPath: string;
  workspacePath: string;
  engine: "atop";
  status: "draft" | "running" | "completed" | "failed" | "cancelled";
  inputFileName: string;
  createdAt: string;
  updatedAt: string;
  completedAt: string | null;
  analysisCount: number;
  analyses?: ProjectAnalysis[];
};

export type ProjectFolderInspection = {
  analysisCount: number;
  existing: boolean;
  folderPath: string;
  hasOtherFiles: boolean;
  historyCount: number;
  name: string;
  project?: Project;
};

export type HistoryEntry = {
  id: string;
  label: string;
  detail: string;
  view: View;
  createdAt: string;
};

export type Artifact = {
  id: string;
  title: string;
  downloadLabel?: string;
  paths: Partial<Record<"svg" | "png" | "pdf", string>>;
};

export type AnalysisValues = {
  comparison_mode: string;
  permutation_seed: number;
  gRNA_col: string;
  gene_col: string;
  sequence_col: string;
  condition1_replicate_cols: string;
  condition2_replicate_cols: string;
  ncs_delta0: number;
  ncs_delta1: number;
  ncs_delta2: number;
  ncs_threshold: number;
  N_perm: number;
  pseudo_count_lfc: number;
  min_sgrna_threshold: number;
  user_engine_choice: string;
};

export type WorkspaceState = {
  activeAnalysisId: string | null;
  activeGseaId: string | null;
  screenPath: string;
  inspection: FileInspection | null;
  analysis: ProcessResult | null;
  gsea: GseaResult | null;
  gseaRun?: {gmtPath: string; config: Record<string, unknown>} | null;
  visualization: VisualizationResult | null;
};

export type SavedProjectWorkspace = WorkspaceState & {
  project: Project;
  artifact: Artifact | null;
  history: HistoryEntry[];
  restore: {
    analysisStateRds: string | null;
    gseaStateRds: string | null;
  };
};
