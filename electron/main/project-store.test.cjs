const assert = require("node:assert/strict");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const test = require("node:test");
const { createProjectStore } = require("./project-store.cjs");

function fixture() {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), "atop-project-store-"));
  const userData = path.join(root, "user-data");
  const inputs = path.join(root, "fixtures");
  fs.mkdirSync(inputs, { recursive: true });
  const sourcePath = path.join(inputs, "screen-数据.tsv");
  const analysisStateRds = path.join(inputs, "analysis_state.rds");
  const gmtPath = path.join(inputs, "pathways.gmt");
  fs.writeFileSync(sourcePath, "sgRNA\tGene\tT1\tC1\ng1\tTP53\t10\t2\n");
  fs.writeFileSync(analysisStateRds, "rds-fixture");
  fs.writeFileSync(gmtPath, "PATHWAY\tdescription\tTP53\n");
  let projectIndex = 0;
  return {
    root,
    sourcePath,
    analysisStateRds,
    gmtPath,
    projectFolder(name) {
      projectIndex += 1;
      const folder = path.join(root, `${name}-${projectIndex}`);
      fs.mkdirSync(folder, { recursive: true });
      return folder;
    },
    store: createProjectStore({ app: { getPath: () => userData } }),
  };
}

test("project store creates the reference workspace structure", (context) => {
  const data = fixture();
  context.after(() => fs.rmSync(data.root, { recursive: true, force: true }));
  const projectFolder = data.projectFolder("Primary screen");
  const project = data.store.createProject({ folderPath: projectFolder });

  assert.equal(project.engine, "atop");
  assert.equal(project.status, "draft");
  assert.equal(project.name, path.basename(projectFolder));
  assert.equal(project.workspacePath, projectFolder);
  for (const directory of ["input", "screen", "gsea", "plots", "temp"]) {
    assert.ok(fs.statSync(path.join(project.workspacePath, directory)).isDirectory());
  }
  assert.ok(fs.existsSync(path.join(project.workspacePath, "manifest.json")));
  assert.deepEqual(data.store.listProjects().map(({ id }) => id), [project.id]);
});

test("project store refuses to create a project without a selected folder", (context) => {
  const data = fixture();
  context.after(() => fs.rmSync(data.root, { recursive: true, force: true }));

  assert.throws(() => data.store.createProject({ name: "Missing folder" }), /project folder/i);
  assert.equal(data.store.listProjects().length, 0);
});

test("project store persists and reopens analysis, GSEA, history, and plot state", (context) => {
  const data = fixture();
  context.after(() => fs.rmSync(data.root, { recursive: true, force: true }));
  const project = data.store.createProject({ name: "Persistent screen", folderPath: data.projectFolder("persistent-folder") });
  const history = [{ id: "run-1", label: "Screen analysis", detail: "1 gene row", view: "results", createdAt: new Date().toISOString() }];
  const inspection = { path: data.sourcePath, fileName: path.basename(data.sourcePath), columns: ["sgRNA", "Gene", "T1", "C1"], preview: [], previewRows: 1, fileSize: 10 };
  const analysis = {
    status: "completed",
    summary: { processedSgrnaRows: 1, geneRows: 1, filteredRows: 0, parameters: {} },
    tables: { sgrna: [{ sgRNA: "g1" }], genes: [{ Gene: "TP53" }], filtered: [] },
    artifacts: { analysisStateRds: data.analysisStateRds },
  };
  const savedAnalysis = data.store.saveAnalysis({ projectId: project.id, sourcePath: data.sourcePath, inspection, analysis, analysisName: "Run one", config: { N_perm: 100000 }, history });

  assert.equal(savedAnalysis.project.status, "completed");
  assert.equal(savedAnalysis.project.analysisCount, 1);
  assert.ok(savedAnalysis.analysis.artifacts.analysisStateRds.startsWith(project.workspacePath));
  assert.ok(savedAnalysis.inspection.path.startsWith(project.workspacePath));

  const gsea = { status: "completed", positive: [{ ID: "PATHWAY" }], negative: [], artifacts: {} };
  const savedGsea = data.store.saveGsea({ projectId: project.id, analysisId: savedAnalysis.analysisId, gsea, gmtPath: data.gmtPath, config: { clp_n_permutations: 1000 }, history });
  const plotSource = path.join(data.root, "plot.svg");
  fs.writeFileSync(plotSource, "<svg xmlns=\"http://www.w3.org/2000/svg\"/>");
  data.store.saveVisualization({
    projectId: project.id,
    visualization: { status: "completed", plotType: "gene_volcano", title: "Volcano", assets: { svg: plotSource } },
    artifact: { id: "plot-1", title: "Volcano", paths: { svg: plotSource } },
  });

  const opened = data.store.openProject(project.id, { analysisId: savedAnalysis.analysisId, gseaId: savedGsea.gseaId });
  assert.equal(opened.workspace.analysis.tables.genes[0].Gene, "TP53");
  assert.equal(opened.workspace.gsea.positive[0].ID, "PATHWAY");
  assert.equal(opened.workspace.visualization.title, "Volcano");
  assert.ok(opened.workspace.artifact.paths.svg.startsWith(project.workspacePath));
  assert.deepEqual(opened.workspace.history, history);
  assert.ok(fs.existsSync(opened.restore.analysisStateRds));
  assert.ok(fs.existsSync(path.join(project.workspacePath, "analyses")));
  assert.ok(fs.existsSync(path.join(project.workspacePath, "analyses", fs.readdirSync(path.join(project.workspacePath, "analyses"))[0], "gsea", savedGsea.gseaId, "pathways.gmt")));
});

test("project store keeps GSEA runs beneath their own screen analysis", (context) => {
  const data = fixture();
  context.after(() => fs.rmSync(data.root, { recursive: true, force: true }));
  const project = data.store.createProject({ name: "Analysis tree", folderPath: data.projectFolder("analysis-tree") });
  const inspection = { path: data.sourcePath, fileName: path.basename(data.sourcePath), columns: [], preview: [], previewRows: 0, fileSize: 0 };
  const makeAnalysis = (gene) => ({ status: "completed", summary: { geneRows: 1 }, tables: { genes: [{ Gene: gene }] }, artifacts: {} });
  const first = data.store.saveAnalysis({ projectId: project.id, sourcePath: data.sourcePath, inspection, analysis: makeAnalysis("TP53"), analysisName: "Screen one" });
  const second = data.store.saveAnalysis({ projectId: project.id, sourcePath: data.sourcePath, inspection, analysis: makeAnalysis("KRAS"), analysisName: "Screen two" });
  data.store.saveGsea({ projectId: project.id, analysisId: first.analysisId, gsea: { status: "completed", positive: [{ ID: "P1" }], negative: [], artifacts: {} }, gmtPath: data.gmtPath });
  const secondGsea = data.store.saveGsea({ projectId: project.id, analysisId: second.analysisId, gsea: { status: "completed", positive: [], negative: [{ ID: "P2" }], artifacts: {} }, gmtPath: data.gmtPath });

  const overview = data.store.listProjects()[0];
  assert.equal(overview.analyses.length, 2);
  assert.equal(overview.analyses.find((analysis) => analysis.id === first.analysisId).gseaRuns.length, 1);
  assert.equal(overview.analyses.find((analysis) => analysis.id === second.analysisId).gseaRuns.length, 1);
  const reopened = data.store.openProject(project.id, { analysisId: second.analysisId, gseaId: secondGsea.gseaId });
  assert.equal(reopened.workspace.analysis.tables.genes[0].Gene, "KRAS");
  assert.equal(reopened.workspace.gsea.negative[0].ID, "P2");
});

test("project folders are inspected and reopened with their saved history", (context) => {
  const data = fixture();
  context.after(() => fs.rmSync(data.root, { recursive: true, force: true }));
  const projectFolder = data.projectFolder("history-project");
  const project = data.store.createProject({ folderPath: projectFolder });
  const history = [{ id: "history-1", label: "Screen analysis", detail: "saved", view: "results", createdAt: new Date().toISOString() }];
  const inspection = { path: data.sourcePath, fileName: path.basename(data.sourcePath), columns: [], preview: [], previewRows: 0, fileSize: 0 };
  data.store.saveAnalysis({
    projectId: project.id,
    sourcePath: data.sourcePath,
    inspection,
    analysis: { status: "completed", summary: { geneRows: 1 }, tables: { genes: [{ Gene: "TP53" }] }, artifacts: {} },
    analysisName: "Saved screen",
    history,
  });

  const secondUserData = path.join(data.root, "second-user-data");
  const secondStore = createProjectStore({ app: { getPath: () => secondUserData } });
  const folderInspection = secondStore.inspectProjectFolder(projectFolder);
  assert.equal(folderInspection.existing, true);
  assert.equal(folderInspection.analysisCount, 1);
  assert.equal(folderInspection.historyCount, 1);

  const reopenedProject = secondStore.createProject({ folderPath: projectFolder });
  const reopened = secondStore.openProject(reopenedProject.id);
  assert.deepEqual(reopened.workspace.history, history);
  assert.equal(reopened.workspace.analysis.tables.genes[0].Gene, "TP53");
});

test("removing a project unregisters it without deleting its folder", (context) => {
  const data = fixture();
  context.after(() => fs.rmSync(data.root, { recursive: true, force: true }));
  const projectFolder = data.projectFolder("keep-files");
  const project = data.store.createProject({ folderPath: projectFolder });

  data.store.deleteProject(project.id);

  assert.equal(data.store.listProjects().length, 0);
  assert.ok(fs.existsSync(path.join(projectFolder, "manifest.json")));
});
