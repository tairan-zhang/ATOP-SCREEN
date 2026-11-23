const fs = require("node:fs");
const path = require("node:path");
const { randomUUID } = require("node:crypto");

const PROJECT_VERSION = 2;

function safeName(value, fallback = "screen-project") {
  const normalized = String(value || fallback).trim().replace(/[^\p{L}\p{N}._-]+/gu, "-");
  return normalized.replace(/^-+|-+$/g, "").slice(0, 80) || fallback;
}

function readJson(filePath, fallback) {
  try {
    return JSON.parse(fs.readFileSync(filePath, "utf8"));
  } catch (_) {
    return fallback;
  }
}

function writeJsonAtomic(filePath, value) {
  fs.mkdirSync(path.dirname(filePath), { recursive: true });
  const temporaryPath = `${filePath}.${process.pid}.tmp`;
  fs.writeFileSync(temporaryPath, `${JSON.stringify(value, null, 2)}\n`, "utf8");
  try {
    fs.renameSync(temporaryPath, filePath);
  } catch (error) {
    if (!error || !["EEXIST", "EPERM"].includes(error.code)) throw error;
    fs.rmSync(filePath, { force: true });
    fs.renameSync(temporaryPath, filePath);
  }
}

function copyFile(sourcePath, targetDirectory, preferredName) {
  if (!sourcePath || !fs.existsSync(sourcePath)) return null;
  if (!fs.statSync(sourcePath).isFile()) return null;
  fs.mkdirSync(targetDirectory, { recursive: true });
  const targetPath = path.join(targetDirectory, safeName(preferredName || path.basename(sourcePath), "artifact"));
  if (path.resolve(sourcePath) !== path.resolve(targetPath)) fs.copyFileSync(sourcePath, targetPath);
  return targetPath;
}

function copyArtifacts(artifacts, targetDirectory) {
  if (!artifacts || typeof artifacts !== "object") return undefined;
  const copied = {};
  for (const [key, sourcePath] of Object.entries(artifacts)) {
    if (typeof sourcePath !== "string") continue;
    const targetPath = copyFile(sourcePath, targetDirectory, path.basename(sourcePath));
    if (targetPath) copied[key] = targetPath;
  }
  return Object.keys(copied).length ? copied : undefined;
}

function createProjectStore({ app }) {
  const projectsRoot = path.join(app.getPath("userData"), "projects");
  const registryPath = path.join(projectsRoot, "registry.json");
  fs.mkdirSync(projectsRoot, { recursive: true });

  function readRegistry() {
    const registry = readJson(registryPath, { version: PROJECT_VERSION, projects: [] });
    return Array.isArray(registry.projects) ? registry.projects : [];
  }

  function writeRegistry(projects) {
    writeJsonAtomic(registryPath, { version: PROJECT_VERSION, projects });
  }

  function manifestPathFor(workspacePath) {
    return path.join(workspacePath, "manifest.json");
  }

  function requireProjectFolder(folderPath) {
    if (typeof folderPath !== "string" || !folderPath.trim()) throw new Error("Choose a project folder before continuing");
    const resolved = path.resolve(folderPath);
    if (!fs.existsSync(resolved) || !fs.statSync(resolved).isDirectory()) throw new Error("The selected project folder is not available");
    return resolved;
  }

  function projectDirectories(workspacePath) {
    return {
      input: path.join(workspacePath, "input"),
      screen: path.join(workspacePath, "screen"),
      gsea: path.join(workspacePath, "gsea"),
      analyses: path.join(workspacePath, "analyses"),
      plots: path.join(workspacePath, "plots"),
      temp: path.join(workspacePath, "temp"),
    };
  }

  function registerManifest(manifestPath, projectId) {
    const normalizedManifestPath = path.resolve(manifestPath);
    writeRegistry([
      { id: projectId, manifestPath: normalizedManifestPath },
      ...readRegistry().filter((entry) => entry.id !== projectId && path.resolve(entry.manifestPath || "") !== normalizedManifestPath),
    ]);
  }

  function rebaseWorkspacePaths(value, oldRoot, newRoot) {
    if (!oldRoot || path.resolve(oldRoot) === path.resolve(newRoot)) return value;
    if (typeof value === "string" && path.isAbsolute(value)) {
      const relative = path.relative(path.resolve(oldRoot), path.resolve(value));
      if (relative === "" || (!relative.startsWith("..") && !path.isAbsolute(relative))) return path.join(newRoot, relative);
      return value;
    }
    if (Array.isArray(value)) return value.map((item) => rebaseWorkspacePaths(item, oldRoot, newRoot));
    if (value && typeof value === "object") {
      for (const key of Object.keys(value)) value[key] = rebaseWorkspacePaths(value[key], oldRoot, newRoot);
    }
    return value;
  }

  function getEntry(projectId) {
    return readRegistry().find((entry) => entry.id === projectId) || null;
  }

  function requireManifest(projectId) {
    const entry = getEntry(projectId);
    if (!entry || !entry.manifestPath || !fs.existsSync(entry.manifestPath)) throw new Error("Project workspace was not found");
    const manifest = readJson(entry.manifestPath, null);
    if (!manifest || manifest.id !== projectId) throw new Error("Project manifest is invalid");
    ensureAnalysisRecords(manifest);
    return { entry, manifest };
  }

  function ensureAnalysisRecords(manifest) {
    if (Array.isArray(manifest.analyses)) return manifest.analyses;
    const legacyAnalysis = manifest.screen ? [{
      id: "analysis-legacy",
      name: manifest.analysisName || "Screen analysis",
      directory: "screen",
      screen: manifest.screen,
      input: manifest.input,
      completedAt: manifest.screen.completedAt || manifest.completedAt || null,
      gseaRuns: manifest.gsea ? [{ id: "gsea-legacy", name: "Pathway enrichment", ...manifest.gsea }] : [],
    }] : [];
    manifest.version = PROJECT_VERSION;
    manifest.analyses = legacyAnalysis;
    manifest.currentAnalysisId = legacyAnalysis[0]?.id || null;
    return manifest.analyses;
  }

  function inspectProjectFolder(folderPath) {
    const workspacePath = requireProjectFolder(folderPath);
    const manifestPath = manifestPathFor(workspacePath);
    const folderName = path.basename(workspacePath);
    const entries = fs.readdirSync(workspacePath).filter((entry) => entry !== ".DS_Store");
    if (!fs.existsSync(manifestPath)) {
      return { existing: false, folderPath: workspacePath, name: folderName, historyCount: 0, analysisCount: 0, hasOtherFiles: entries.length > 0 };
    }
    const manifest = readJson(manifestPath, null);
    if (!manifest || typeof manifest.id !== "string" || (manifest.engine && manifest.engine !== "atop")) {
      throw new Error("The selected folder contains an invalid ATOP project manifest");
    }
    const analyses = ensureAnalysisRecords(manifest);
    const historyPath = manifest.historyFile ? path.join(workspacePath, manifest.historyFile) : path.join(workspacePath, "history.json");
    const history = readJson(historyPath, []);
    return {
      existing: true,
      folderPath: workspacePath,
      name: manifest.name || folderName,
      historyCount: Array.isArray(history) ? history.length : 0,
      analysisCount: analyses.length,
      hasOtherFiles: true,
      project: projectSummary({ ...manifest, workspacePath }),
    };
  }

  function projectSummary(manifest) {
    const analyses = ensureAnalysisRecords(manifest).map((analysis) => ({
      id: analysis.id,
      name: analysis.name,
      completedAt: analysis.completedAt || null,
      geneRows: Number(analysis.geneRows || 0),
      gseaRuns: (analysis.gseaRuns || []).map((gsea) => ({ id: gsea.id, name: gsea.name, completedAt: gsea.completedAt || null, pathwayRows: Number(gsea.pathwayRows || 0) })),
    }));
    return {
      id: manifest.id,
      name: manifest.name,
      rootPath: manifest.workspacePath,
      workspacePath: manifest.workspacePath,
      engine: "atop",
      status: manifest.status || "draft",
      analysisCount: Number(manifest.analysisCount || 0),
      analyses,
      inputFileName: manifest.input?.fileName || "",
      createdAt: manifest.createdAt,
      updatedAt: manifest.updatedAt,
      completedAt: manifest.completedAt || null,
    };
  }

  function listProjects() {
    const registeredEntries = readRegistry();
    const validEntries = [];
    const projects = [];
    for (const entry of registeredEntries) {
      if (!entry.manifestPath || !fs.existsSync(entry.manifestPath)) continue;
      const manifest = readJson(entry.manifestPath, null);
      if (!manifest || manifest.id !== entry.id) continue;
      ensureAnalysisRecords(manifest);
      validEntries.push(entry);
      projects.push(projectSummary(manifest));
    }
    if (validEntries.length !== registeredEntries.length) writeRegistry(validEntries);
    return projects.sort((left, right) => String(right.updatedAt).localeCompare(String(left.updatedAt)));
  }

  function createProject(input = {}) {
    const workspacePath = requireProjectFolder(input.folderPath);
    const existingManifestPath = manifestPathFor(workspacePath);
    if (fs.existsSync(existingManifestPath)) {
      const manifest = readJson(existingManifestPath, null);
      if (!manifest || typeof manifest.id !== "string" || (manifest.engine && manifest.engine !== "atop")) {
        throw new Error("The selected folder contains an invalid ATOP project manifest");
      }
      const oldWorkspacePath = manifest.workspacePath;
      rebaseWorkspacePaths(manifest, oldWorkspacePath, workspacePath);
      manifest.workspacePath = workspacePath;
      manifest.directories = projectDirectories(workspacePath);
      for (const directory of Object.values(manifest.directories)) fs.mkdirSync(directory, { recursive: true });
      manifest.version = PROJECT_VERSION;
      ensureAnalysisRecords(manifest);
      writeJsonAtomic(existingManifestPath, manifest);
      registerManifest(existingManifestPath, manifest.id);
      return projectSummary(manifest);
    }
    const id = randomUUID();
    const folderName = path.basename(workspacePath);
    const name = String(input.name || folderName).trim().slice(0, 120) || folderName;
    const directories = projectDirectories(workspacePath);
    for (const directory of Object.values(directories)) fs.mkdirSync(directory, { recursive: true });
    const now = new Date().toISOString();
    const manifest = {
      version: PROJECT_VERSION,
      id,
      name,
      engine: "atop",
      status: "draft",
      analysisCount: 0,
      workspacePath,
      directories,
      input: null,
      screen: null,
      gsea: null,
      analyses: [],
      currentAnalysisId: null,
      plots: [],
      createdAt: now,
      updatedAt: now,
      completedAt: null,
    };
    const manifestPath = manifestPathFor(workspacePath);
    writeJsonAtomic(manifestPath, manifest);
    registerManifest(manifestPath, id);
    return projectSummary(manifest);
  }

  function openProject(projectId, selection = {}) {
    const { manifest } = requireManifest(projectId);
    const readStateFile = (relativePath) => relativePath ? readJson(path.join(manifest.workspacePath, relativePath), null) : null;
    const analyses = ensureAnalysisRecords(manifest);
    const selectedAnalysis = analyses.find((record) => record.id === selection.analysisId)
      || analyses.find((record) => record.id === manifest.currentAnalysisId)
      || analyses.at(-1)
      || null;
    const gseaRuns = selectedAnalysis?.gseaRuns || [];
    const selectedGsea = gseaRuns.find((record) => record.id === selection.gseaId)
      || gseaRuns.find((record) => record.id === selectedAnalysis?.currentGseaId)
      || gseaRuns.at(-1)
      || null;
    const currentGsea = selectedGsea || (!selectedAnalysis ? manifest.gsea : null);
    const selectedPlot = (manifest.plots || []).filter((plot) => !plot.analysisId || plot.analysisId === selectedAnalysis?.id).at(-1)
      || (!manifest.latestPlot?.analysisId || manifest.latestPlot.analysisId === selectedAnalysis?.id ? manifest.latestPlot : null)
      || null;
    return {
      project: projectSummary(manifest),
      workspace: {
        activeAnalysisId: selectedAnalysis?.id || null,
        activeGseaId: selectedGsea?.id || null,
        screenPath: selectedAnalysis?.input?.path || manifest.input?.path || "",
        inspection: readStateFile(selectedAnalysis?.screen?.inspectionFile || manifest.screen?.inspectionFile),
        analysis: readStateFile(selectedAnalysis?.screen?.resultFile || manifest.screen?.resultFile),
        gsea: readStateFile(currentGsea?.resultFile),
        gseaRun: currentGsea ? { gmtPath: currentGsea.sourcePath || "", config: readStateFile(currentGsea.configFile) || {} } : null,
        visualization: readStateFile(selectedPlot?.visualizationFile),
        artifact: readStateFile(selectedPlot?.artifactFile),
        history: readStateFile(manifest.historyFile) || [],
      },
      restore: {
        analysisStateRds: selectedAnalysis?.screen?.artifacts?.analysisStateRds || manifest.screen?.artifacts?.analysisStateRds || null,
        gseaStateRds: currentGsea?.artifacts?.gseaStateRds || null,
      },
    };
  }

  function saveAnalysis(payload) {
    const { manifest } = requireManifest(payload.projectId);
    const analysisId = payload.analysisId || `analysis-${Date.now()}-${randomUUID().slice(0, 6)}`;
    const analysisName = String(payload.analysisName || "Screen analysis").trim().slice(0, 120) || "Screen analysis";
    const analysesRoot = manifest.directories.analyses || path.join(manifest.workspacePath, "analyses");
    manifest.directories.analyses = analysesRoot;
    const analysisDirectory = path.join(analysesRoot, `${safeName(analysisName)}-${analysisId.slice(-6)}`);
    const inputDirectory = path.join(analysisDirectory, "input");
    const screenDirectory = path.join(analysisDirectory, "screen");
    const inputPath = copyFile(payload.sourcePath, inputDirectory, path.basename(payload.sourcePath || "screen-data"));
    if (!inputPath) throw new Error("The analysis input file could not be saved");
    const artifacts = copyArtifacts(payload.analysis?.artifacts, screenDirectory);
    const savedAnalysis = { ...payload.analysis, artifacts };
    const savedInspection = payload.inspection ? { ...payload.inspection, path: inputPath } : null;
    const resultFile = path.join(screenDirectory, "analysis-result.json");
    const inspectionFile = path.join(screenDirectory, "inspection.json");
    const configFile = path.join(screenDirectory, "config.json");
    writeJsonAtomic(resultFile, savedAnalysis);
    writeJsonAtomic(inspectionFile, savedInspection);
    const settings = { ...payload.config, ...savedAnalysis.summary?.parameters };
    writeJsonAtomic(configFile, settings);
    const now = new Date().toISOString();
    writeJsonAtomic(path.join(screenDirectory, "runs", randomUUID() + ".json"), { completedAt: now, analysisId, input: inputPath, settings });
    manifest.status = "completed";
    const screen = {
      resultFile: path.relative(manifest.workspacePath, resultFile),
      inspectionFile: path.relative(manifest.workspacePath, inspectionFile),
      configFile: path.relative(manifest.workspacePath, configFile),
      artifacts: artifacts || {},
      completedAt: now,
    };
    const existingAnalysis = ensureAnalysisRecords(manifest).find((record) => record.id === analysisId);
    const analysisRecord = {
      id: analysisId,
      name: analysisName,
      directory: path.relative(manifest.workspacePath, analysisDirectory),
      input: { path: inputPath, fileName: path.basename(inputPath) },
      screen,
      geneRows: Number(savedAnalysis.summary?.geneRows || 0),
      completedAt: now,
      gseaRuns: existingAnalysis?.gseaRuns || [],
      currentGseaId: existingAnalysis?.currentGseaId || null,
    };
    const records = ensureAnalysisRecords(manifest).filter((record) => record.id !== analysisId);
    manifest.analyses = [...records, analysisRecord];
    manifest.currentAnalysisId = analysisId;
    manifest.analysisCount = manifest.analyses.length;
    manifest.analysisName = analysisName;
    manifest.input = { path: inputPath, fileName: path.basename(inputPath) };
    manifest.screen = screen;
    manifest.updatedAt = now;
    manifest.completedAt = now;
    if (payload.history) {
      const historyFile = path.join(manifest.workspacePath, "history.json");
      writeJsonAtomic(historyFile, payload.history);
      manifest.historyFile = path.relative(manifest.workspacePath, historyFile);
    }
    writeJsonAtomic(manifestPathFor(manifest.workspacePath), manifest);
    return { project: projectSummary(manifest), analysisId, analysis: savedAnalysis, inspection: savedInspection };
  }

  function saveGsea(payload) {
    const { manifest } = requireManifest(payload.projectId);
    const analysis = ensureAnalysisRecords(manifest).find((record) => record.id === (payload.analysisId || manifest.currentAnalysisId));
    if (!analysis) throw new Error("Save a Screen analysis before saving pathway enrichment");
    const gseaId = payload.gseaId || `gsea-${Date.now()}-${randomUUID().slice(0, 6)}`;
    const gseaDirectory = path.join(manifest.workspacePath, analysis.directory, "gsea", gseaId);
    const sourcePath = copyFile(payload.gmtPath, gseaDirectory, path.basename(payload.gmtPath || "gene-sets.gmt"));
    const artifacts = copyArtifacts(payload.gsea?.artifacts, gseaDirectory);
    const savedGsea = { ...payload.gsea, artifacts };
    const resultFile = path.join(gseaDirectory, "gsea-result.json");
    const configFile = path.join(gseaDirectory, "config.json");
    writeJsonAtomic(resultFile, savedGsea);
    const settings = { ...payload.config, ...savedGsea.parameters };
    writeJsonAtomic(configFile, settings);
    writeJsonAtomic(path.join(gseaDirectory, "runs", randomUUID() + ".json"), { completedAt: new Date().toISOString(), gseaId, settings });
    const now = new Date().toISOString();
    const gseaRecord = {
      id: gseaId,
      name: String(payload.gseaName || "Pathway enrichment").trim().slice(0, 120) || "Pathway enrichment",
      resultFile: path.relative(manifest.workspacePath, resultFile),
      configFile: path.relative(manifest.workspacePath, configFile),
      artifacts: artifacts || {},
      sourcePath,
      completedAt: now,
      pathwayRows: Number((savedGsea.positive?.length || 0) + (savedGsea.negative?.length || 0)),
    };
    analysis.gseaRuns = [...(analysis.gseaRuns || []).filter((record) => record.id !== gseaId), gseaRecord];
    analysis.currentGseaId = gseaId;
    manifest.gsea = gseaRecord;
    manifest.updatedAt = now;
    if (payload.history) {
      const historyFile = path.join(manifest.workspacePath, "history.json");
      writeJsonAtomic(historyFile, payload.history);
      manifest.historyFile = path.relative(manifest.workspacePath, historyFile);
    }
    writeJsonAtomic(manifestPathFor(manifest.workspacePath), manifest);
    return { project: projectSummary(manifest), analysisId: analysis.id, gseaId, gsea: savedGsea };
  }

  function saveVisualization(payload) {
    const { manifest } = requireManifest(payload.projectId);
    const createdAt = new Date().toISOString();
    const plotId = `plot-${Date.now()}-${randomUUID().slice(0, 6)}`;
    const activeAnalysis = ensureAnalysisRecords(manifest).find((record) => record.id === (payload.analysisId || manifest.currentAnalysisId));
    const plotsRoot = activeAnalysis
      ? path.join(manifest.workspacePath, activeAnalysis.directory, "plots")
      : manifest.directories.plots;
    const plotDirectory = path.join(plotsRoot, plotId);
    const assets = copyArtifacts(payload.visualization?.assets, plotDirectory) || {};
    const savedVisualization = { ...payload.visualization, assets };
    const visualizationFile = path.join(plotDirectory, "visualization.json");
    const artifactFile = path.join(plotDirectory, "artifact.json");
    const savedArtifact = payload.artifact ? { ...payload.artifact, paths: assets } : null;
    writeJsonAtomic(visualizationFile, savedVisualization);
    writeJsonAtomic(artifactFile, savedArtifact);
    const plot = {
      id: plotId,
      type: savedVisualization.plotType,
      title: savedVisualization.title,
      assets,
      visualizationFile: path.relative(manifest.workspacePath, visualizationFile),
      artifactFile: path.relative(manifest.workspacePath, artifactFile),
      analysisId: activeAnalysis?.id || null,
      createdAt,
    };
    manifest.plots = [...(Array.isArray(manifest.plots) ? manifest.plots : []), plot];
    manifest.latestPlot = plot;
    manifest.updatedAt = createdAt;
    writeJsonAtomic(manifestPathFor(manifest.workspacePath), manifest);
    return { project: projectSummary(manifest), visualization: savedVisualization, artifact: savedArtifact };
  }

  function deleteProject(projectId) {
    const { entry } = requireManifest(projectId);
    writeRegistry(readRegistry().filter((candidate) => candidate.manifestPath !== entry.manifestPath));
    return true;
  }

  function listWorkspaceRoots() {
    return listProjects().map((project) => path.resolve(project.workspacePath));
  }

  return { createProject, inspectProjectFolder, listProjects, openProject, saveAnalysis, saveGsea, saveVisualization, deleteProject, listWorkspaceRoots };
}

module.exports = { createProjectStore, safeName, writeJsonAtomic };
