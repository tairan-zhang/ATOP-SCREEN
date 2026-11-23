import type { AnalysisValues, HistoryEntry, Project, Theme } from "../types/atop";

export const storageKeys = {
  project: "atop-screen.project",
  history: "atop-screen.history",
  theme: "atop-screen.theme",
};

export const defaultProject: Project = {
  id: "default-project",
  name: "Untitled screen project",
  rootPath: "",
  workspacePath: "",
  engine: "atop",
  status: "draft",
  inputFileName: "",
  createdAt: new Date().toISOString(),
  updatedAt: new Date().toISOString(),
  completedAt: null,
  analysisCount: 0,
};

export const defaultAnalysisValues: AnalysisValues = {
  comparison_mode: "auto",
  permutation_seed: 42,
  gRNA_col: "gRNA",
  gene_col: "Gene",
  sequence_col: "",
  condition1_replicate_cols: "Day7T_1",
  condition2_replicate_cols: "Day7C_1",
  ncs_delta0: 0.1,
  ncs_delta1: 0.3,
  ncs_delta2: 0.2,
  ncs_threshold: 0.1,
  N_perm: 100000,
  pseudo_count_lfc: 1,
  min_sgrna_threshold: 3,
  user_engine_choice: "cpp",
};

export function readStorage<T>(key: string, fallback: T): T {
  try {
    const raw = window.localStorage.getItem(key);
    return raw ? (JSON.parse(raw) as T) : fallback;
  } catch {
    return fallback;
  }
}

export function writeStorage<T>(key: string, value: T) {
  try {
    window.localStorage.setItem(key, JSON.stringify(value));
  } catch {
    // Storage is a convenience; analysis state must not depend on it.
  }
}

export type StoredWorkspace = {
  project: Project;
  history: HistoryEntry[];
  theme: Theme;
};
