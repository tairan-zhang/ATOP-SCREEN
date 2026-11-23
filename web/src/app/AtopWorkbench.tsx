import { useCallback, useEffect, useState } from "react";
import { AnimatedSidebarProvider, AnimatedSidebarInset } from "@atop/ui/motion/animated-sidebar";
import { ParameterPane } from "../components/ParameterPane";
import { DataTable } from "../components/DataTable";
import { Icon } from "../components/Icon";
import type { FileInspection, GseaResult, HealthResult, ProcessResult, VisualizationResult } from "@atop-screen/ipc-contracts";
import { desktopApi } from "../services/desktop-api";
import { defaultProject, readStorage, storageKeys, writeStorage } from "../domain/storage";
import type { Artifact, HistoryEntry, PlotSection, Project, SavedProjectWorkspace, Theme, View } from "../types/atop";
import { ArtifactPanel } from "../components/ArtifactPanel";
import { ProjectDialog, type ProjectDraft } from "../components/ProjectDialog";
import { ProjectLibrary } from "../components/ProjectLibrary";
import { Sidebar } from "../components/Sidebar";
import { WorkspaceTitlebar } from "../components/WorkspaceTitlebar";
import { AnalysisView, type AnalysisRunMetadata } from "../features/analysis/AnalysisView";
import { GseaResultsView, GseaView, type GseaRunMetadata } from "../features/gsea/GseaView";
import { HomeView } from "../features/home/HomeView";
import { ResultPicker, ResultsView } from "../features/results/ResultsView";
import { SettingsView } from "../features/settings/SettingsView";
import { VisualizationView } from "../features/visualization/VisualizationView";

const pageTitles: Record<View, string> = {
  home: "Workspace", analysis: "Screen analysis", gsea: "Pathway enrichment", gseaResults: "Pathway results",
  visualize: "Visualization studio", results: "Screen results", settings: "Runtime settings",
};
const plotTitles: Record<PlotSection, string> = {gene: "Gene plots", sgrna: "sgRNA plots", pathway: "Pathway plots"};
const analysisRequiredViews = new Set<View>(["gsea", "gseaResults", "visualize", "results"]);

function newHistoryEntry(label: string, detail: string, view: View): HistoryEntry {
  return { id: `${Date.now()}-${Math.random().toString(36).slice(2, 7)}`, label, detail, view, createdAt: new Date().toISOString() };
}

export function AtopWorkbench() {
  const [view, setView] = useState<View>("analysis");
  const [theme, setTheme] = useState<Theme>(() => readStorage<Theme>(storageKeys.theme, "light"));
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
  const [sidebarWidth, setSidebarWidth] = useState(196);
  const [artifactOpen, setArtifactOpen] = useState(false);
  const [artifact, setArtifact] = useState<Artifact | null>(null);
  const [gseaRun, setGseaRun] = useState<GseaRunMetadata | null>(null);
  const [project, setProject] = useState<Project>(defaultProject);
  const [projects, setProjects] = useState<Project[]>([]);
  const [projectBusyId, setProjectBusyId] = useState("");
  const [projectLibraryOpen, setProjectLibraryOpen] = useState(false);
  const [history, setHistory] = useState<HistoryEntry[]>([]);
  const [screenPath, setScreenPath] = useState("");
  const [activeAnalysisId, setActiveAnalysisId] = useState<string | null>(null);
  const [activeGseaId, setActiveGseaId] = useState<string | null>(null);
  const [inspection, setInspection] = useState<FileInspection | null>(null);
  const [analysis, setAnalysis] = useState<ProcessResult | null>(null);
  const [gsea, setGsea] = useState<GseaResult | null>(null);
  const [visualization, setVisualization] = useState<VisualizationResult | null>(null);
  const [resultTab, setResultTab] = useState("genes");
  const [ncsChart, setNcsChart] = useState("density");
  const [plotSection, setPlotSection] = useState<PlotSection>("gene");
  const [health, setHealth] = useState<HealthResult | null>(null);
  const [progress, setProgress] = useState("Runtime standing by");
  const [projectDialogOpen, setProjectDialogOpen] = useState(false);

  useEffect(() => {
    document.documentElement.classList.toggle("dark", theme === "dark");
    document.documentElement.style.colorScheme = theme;
    writeStorage(storageKeys.theme, theme);
  }, [theme]);
  useEffect(() => {
    try {
      window.localStorage.removeItem(storageKeys.project);
      window.localStorage.removeItem(storageKeys.history);
    } catch {
      // Project state is read from the selected folder, never browser storage.
    }
  }, []);
  const refreshHealth = useCallback(async () => {
    setProgress("Checking the analysis runtime…");
    try {
      const result = await desktopApi.request<HealthResult>({ command: "health" });
      setHealth(result);
      setProgress(result.status === "ready" ? "Analysis runtime ready" : result.rVersion);
    } catch (error) { setProgress(error instanceof Error ? error.message : String(error)); }
  }, []);

  useEffect(() => {
    void refreshHealth();
    let timer: number | null = null;
    let pendingMessage = "";
    const flush = () => {
      timer = null;
      setProgress((current) => current === pendingMessage ? current : pendingMessage);
    };
    const unsubscribe = desktopApi.subscribeProgress((event) => {
      if (!event.message) return;
      pendingMessage = event.message;
      if (timer === null) timer = window.setTimeout(flush, 120);
    });
    return () => {
      unsubscribe();
      if (timer !== null) window.clearTimeout(timer);
    };
  }, [refreshHealth]);

  const refreshProjects = useCallback(async () => {
    if (!desktopApi.isDesktop()) return [];
    const savedProjects = await desktopApi.projects.list();
    setProjects(savedProjects);
    return savedProjects;
  }, []);

  const applyWorkspace = useCallback((opened: SavedProjectWorkspace) => {
    setResultTab("genes");
    setNcsChart("density");
    setProject(opened.project);
    setScreenPath(opened.screenPath);
    setActiveAnalysisId(opened.activeAnalysisId);
    setActiveGseaId(opened.activeGseaId);
    setInspection(opened.inspection);
    setAnalysis(opened.analysis);
    setGsea(opened.gsea);
    setGseaRun(opened.gseaRun ?? null);
    setVisualization(opened.visualization);
    setHistory(opened.history);
    setArtifact(opened.artifact);
    setArtifactOpen(false);
    setView(opened.analysis ? "results" : "analysis");
  }, []);

  const openSavedProject = useCallback(async (projectId: string, closeLibrary = true, selection?: { analysisId?: string; gseaId?: string }) => {
    if (!desktopApi.isDesktop()) return;
    setProjectBusyId(projectId);
    setProgress("Opening the saved project…");
    try {
      const opened = await desktopApi.projects.open(projectId, selection);
      applyWorkspace({ project: opened.project, ...opened.workspace, restore: opened.restore });
      if (opened.restore.analysisStateRds || opened.restore.gseaStateRds) {
        await desktopApi.request({ command: "restore_state", parameters: opened.restore });
      }
      if (closeLibrary) setProjectLibraryOpen(false);
      setProgress("Project opened.");
    } catch (error) { setProgress(error instanceof Error ? error.message : String(error)); }
    finally { setProjectBusyId(""); }
  }, [applyWorkspace]);

  useEffect(() => {
    if (!desktopApi.isDesktop()) return;
    let cancelled = false;
    void (async () => {
      try {
        const savedProjects = await refreshProjects();
        if (cancelled) return;
        const remembered = savedProjects.find((candidate) => candidate.id === project.id);
        if (remembered) return void await openSavedProject(remembered.id, false);
        if (savedProjects[0]) return void await openSavedProject(savedProjects[0].id, false);
        if (!cancelled) {
          setProject(defaultProject);
          setProjects([]);
          setHistory([]);
          setScreenPath("");
          setInspection(null);
          setAnalysis(null);
          setGsea(null); setGseaRun(null);
          setVisualization(null);
          setArtifact(null);
          setArtifactOpen(false);
          setProjectDialogOpen(true);
          setProgress("Choose a project folder to begin");
        }
      } catch (error) { if (!cancelled) setProgress(error instanceof Error ? error.message : String(error)); }
    })();
    return () => { cancelled = true; };
    // Initial desktop workspace resolution intentionally runs only once.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const openArtifact = useCallback((nextArtifact: Artifact) => { setArtifact(nextArtifact); setArtifactOpen(true); }, []);
  const navigate = useCallback((nextView: View) => {
    if (analysisRequiredViews.has(nextView) && !analysis) { setProgress("Complete Screen Analysis first"); return; }
    setArtifactOpen(nextView === "visualize");
    setView(nextView);
  }, [analysis]);

  const startNewAnalysis = useCallback(() => {
    if (!project.workspacePath) {
      setProjectDialogOpen(true);
      setProgress("Choose a project folder before adding a Screen analysis.");
      return;
    }
    setActiveAnalysisId(null);
    setActiveGseaId(null);
    setScreenPath("");
    setInspection(null);
    setAnalysis(null);
    setGsea(null); setGseaRun(null);
    setVisualization(null);
    setArtifact(null);
    setArtifactOpen(false);
    setPlotSection("gene");
    setView("analysis");
    setProgress("New analysis ready.");
  }, [project.workspacePath]);

  const clearWorkspace = useCallback((nextProject: Project) => {
    setProject(nextProject); setHistory([]); setScreenPath(""); setInspection(null); setAnalysis(null); setGsea(null); setGseaRun(null);
    setActiveAnalysisId(null); setActiveGseaId(null);
    setVisualization(null); setArtifact(null); setArtifactOpen(false); setPlotSection("gene"); setView("analysis");
  }, []);

  async function createProject({ name, rootPath }: ProjectDraft) {
    if (!rootPath) {
      setProgress("Choose a project folder before continuing.");
      return;
    }
    setProgress("Checking the selected project folder…");
    try {
      const now = new Date().toISOString();
      const created = desktopApi.isDesktop()
        ? await desktopApi.projects.create({ name, folderPath: rootPath })
        : { ...defaultProject, id: `project-${Date.now()}`, name, rootPath, workspacePath: rootPath, createdAt: now, updatedAt: now };
      if (desktopApi.isDesktop()) {
        await refreshProjects();
        await openSavedProject(created.id, false);
      } else {
        clearWorkspace(created);
      }
      setProjectDialogOpen(false);
      setProgress(created.analysisCount > 0 ? "Project opened." : "Project ready.");
    } catch (error) { setProgress(error instanceof Error ? error.message : String(error)); }
  }

  async function handleAnalysis(result: ProcessResult, metadata: AnalysisRunMetadata) {
    const nextHistory = [newHistoryEntry("Screen analysis", `${result.summary.geneRows.toLocaleString()} gene rows`, "results"), ...history].slice(0, 20);
    setResultTab("genes");
    setNcsChart("density");
    setAnalysis(result); setArtifact(null); setArtifactOpen(false); setGsea(null); setGseaRun(null); setVisualization(null); setHistory(nextHistory); setView("results");
    if (!desktopApi.isDesktop() || project.id.startsWith("project-")) {
      const now = new Date().toISOString();
      setProject((current) => ({ ...current, status: "completed", analysisCount: current.analysisCount + 1, updatedAt: now, completedAt: now }));
      return;
    }
    try {
      setProgress("Saving analysis into the project workspace…");
      const saved = await desktopApi.projects.saveAnalysis({ projectId: project.id, analysisId: activeAnalysisId || undefined, sourcePath: screenPath, inspection, analysis: result, analysisName: metadata.analysisName, config: metadata.config, history: nextHistory });
      setProject(saved.project); setAnalysis(saved.analysis); setInspection(saved.inspection); setScreenPath(saved.inspection.path);
      setActiveAnalysisId(saved.analysisId); setActiveGseaId(null);
      await refreshProjects(); setProgress("Analysis saved.");
    } catch (error) { setProgress(`Analysis completed, but project saving failed: ${error instanceof Error ? error.message : String(error)}`); }
  }

  async function handleGsea(result: GseaResult, metadata: GseaRunMetadata) {
    const nextHistory = [newHistoryEntry("GSEA enrichment", `${result.positive.length + result.negative.length} pathway rows`, "gsea"), ...history].slice(0, 20);
    setGsea(result); setGseaRun(metadata); setArtifactOpen(false); setView("gseaResults"); setHistory(nextHistory);
    if (!desktopApi.isDesktop() || project.id.startsWith("project-")) return;
    try {
      setProgress("Saving GSEA output into the project workspace…");
      const saved = await desktopApi.projects.saveGsea({ projectId: project.id, analysisId: activeAnalysisId || undefined, gseaId: activeGseaId || undefined, gsea: result, gmtPath: metadata.gmtPath, config: metadata.config, history: nextHistory });
      setProject(saved.project); setGsea(saved.gsea); await refreshProjects(); setProgress("GSEA results saved.");
      setActiveAnalysisId(saved.analysisId); setActiveGseaId(saved.gseaId);
    } catch (error) { setProgress(`GSEA completed, but project saving failed: ${error instanceof Error ? error.message : String(error)}`); }
  }

  async function handleVisualization(result: VisualizationResult, nextArtifact: Artifact) {
    setVisualization(result);
    if (!desktopApi.isDesktop() || project.id.startsWith("project-")) return;
    try {
      const saved = await desktopApi.projects.saveVisualization({ projectId: project.id, analysisId: activeAnalysisId || undefined, visualization: result, artifact: nextArtifact });
      setProject(saved.project); setVisualization(saved.visualization); if (saved.artifact) openArtifact(saved.artifact);
      await refreshProjects(); setProgress("Figure saved.");
    } catch (error) { setProgress(`Figure rendered, but project saving failed: ${error instanceof Error ? error.message : String(error)}`); }
  }

  async function deleteSavedProject(target: Project) {
    setProjectBusyId(target.id);
    try {
      await desktopApi.projects.delete(target.id);
      const remaining = await refreshProjects();
      if (target.id === project.id) {
        if (remaining[0]) await openSavedProject(remaining[0].id, false);
        else { clearWorkspace(defaultProject); setProjects([]); setProjectDialogOpen(true); }
      }
      setProgress("Project removed from the app. Files remain in its project folder.");
    } catch (error) { setProgress(error instanceof Error ? error.message : String(error)); }
    finally { setProjectBusyId(""); }
  }

  async function showProjectLibrary() {
    try { await refreshProjects(); setProjectLibraryOpen(true); }
    catch (error) { setProgress(error instanceof Error ? error.message : String(error)); }
  }

  const screenModule = view === "analysis" || view === "results";
  const pathwayModule = view === "gsea" || view === "gseaResults";
  const analysisName = project.analyses?.find(item => item.id === activeAnalysisId)?.name ?? project.name;
  const parameters = view === "results" ? (
    <ResultPicker analysis={analysis} tab={resultTab} chart={ncsChart} onTab={setResultTab} onChart={setNcsChart}/>
  ) : view === "analysis" ? (
    <AnalysisView key={project.id + (activeAnalysisId ?? "new")} analysis={analysis} inspection={inspection}
      onAnalysis={(result, metadata) => void handleAnalysis(result, metadata)} onInspection={setInspection}
      onPath={path => { setScreenPath(path); setView("analysis"); setArtifactOpen(false); }}
      onProgress={setProgress} path={screenPath} projectName={analysisName}/>
  ) : pathwayModule ? (
    <GseaView key={project.id + activeAnalysisId + activeGseaId} analysis={analysis} initialRun={gseaRun}
      onGsea={(result, metadata) => void handleGsea(result, metadata)} onProgress={setProgress}/>
  ) : view === "visualize" ? (
    <VisualizationView analysis={analysis} gsea={gsea} inspection={inspection} onArtifact={openArtifact}
      onProgress={setProgress} onVisualization={(result, figure) => void handleVisualization(result, figure)} plotSection={plotSection}/>
  ) : null;

  const plotPrefix = {gene: "gene_", sgrna: "sgrna_", pathway: "gsea_"}[plotSection];
  const matchingArtifact = view !== "visualize" || !visualization || visualization.plotType.startsWith(plotPrefix) ? artifact : null;
  let content;
  if (artifactOpen) {
    content = <ArtifactPanel artifact={matchingArtifact} projectName={project.name} projectPath={project.workspacePath}/>;
  } else if (view === "home") {
    content = <HomeView projectName={project.name} analysisReady={Boolean(analysis)} onNavigate={navigate}/>;
  } else if (view === "settings") {
    content = <SettingsView health={health} onRefresh={() => void refreshHealth()} setTheme={setTheme} theme={theme}/>;
  } else if (view === "results") {
    content = <ResultsView analysis={analysis} tab={resultTab} chart={ncsChart}/>;
  } else if (pathwayModule) {
    content = <GseaResultsView gsea={gsea}/>;
  } else if (view === "analysis" && inspection) {
    content = <section className="input-preview"><header className="canvas-heading"><Icon name="database"/><h2>Input preview</h2><span title={inspection.fileName}>{inspection.fileName}</span></header><DataTable rows={inspection.preview} empty="No preview rows."/></section>;
  } else {
    content = <div className="output-empty"><Icon name={screenModule ? "database" : "image"} size={26}/>
      <span>{screenModule ? "Choose a count matrix to preview its rows." : "Figure preview hidden."}</span></div>;
  }

  return <AnimatedSidebarProvider className="atop-shell" motionPreset="smooth" style={{"--sidebar-width": sidebarWidth + "px"}}
    open={!sidebarCollapsed} onOpenChange={open => setSidebarCollapsed(!open)}>
    <div className="product-shell" data-product-mode="bio">
      <WorkspaceTitlebar projectName={project.name} pageTitle={view === "visualize" ? plotTitles[plotSection] : pageTitles[view]} sidebarOpen={!sidebarCollapsed} previewOpen={artifactOpen}
        onToggleSidebar={() => setSidebarCollapsed(open => !open)} onTogglePreview={() => setArtifactOpen(open => !open)} />
      <div className="product-base">
          <Sidebar open={!sidebarCollapsed} onWidth={setSidebarWidth} project={project} analysis={analysis} gsea={gsea} view={view} plotSection={plotSection}
            onCreateProject={() => setProjectDialogOpen(true)} onNavigate={navigate} onNewAnalysis={startNewAnalysis} onOpenLibrary={() => void showProjectLibrary()}
            onOpenAnalysis={analysisId => void openSavedProject(project.id, false, {analysisId})}
            onOpenGsea={(analysisId,gseaId) => void openSavedProject(project.id, false, {analysisId,gseaId}).then(() => setView("gseaResults"))} onPlotSection={setPlotSection}/>
        <div className="product-surface"><div className="product-bio-workspace">
          {parameters && <ParameterPane>{parameters}</ParameterPane>}
          <AnimatedSidebarInset className="workspace-main">
            <div className="product-page screen-page result-canvas" role="region" aria-label={parameters ? "Results" : pageTitles[view]}>{content}</div>
          </AnimatedSidebarInset>
        </div></div>
      </div>
      <footer role="status" className="screen-status">{progress}</footer>
      <ProjectDialog open={projectDialogOpen} onClose={() => setProjectDialogOpen(false)} onCreate={(draft) => void createProject(draft)} required={!project.workspacePath} />
      {projectLibraryOpen ? <ProjectLibrary busyId={projectBusyId} onClose={() => setProjectLibraryOpen(false)} onDelete={(target) => void deleteSavedProject(target)} onOpen={(target) => void openSavedProject(target.id)} onOpenAnalysis={(target, analysisId) => void openSavedProject(target.id, true, { analysisId }).then(() => setView("results"))} onOpenGsea={(target, analysisId, gseaId) => void openSavedProject(target.id, true, { analysisId, gseaId }).then(() => setView("gseaResults"))} projects={projects} /> : null}
    </div>
  </AnimatedSidebarProvider>;
}
