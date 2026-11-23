import {Button} from "@atop/ui";
import { useRef, useState } from "react";
import type { ChangeEvent } from "react";
import type { GseaResult, ProcessResult } from "@atop-screen/ipc-contracts";
import { desktopApi } from "../../services/desktop-api";
import { fileName } from "../../domain/format";
import { DataTable } from "../../components/DataTable";
import { Icon } from "../../components/Icon";
import { ProcessingIndicator } from "../../components/ProcessingIndicator";
import { NumberField, SelectField } from "../../components/FormFields";
import { ErrorBanner, Prerequisite } from "../../components/WorkspacePrimitives";

export type GseaRunMetadata = { gmtPath: string; config: Record<string, unknown> };

export function GseaView({ analysis, initialRun, onGsea, onProgress }: { initialRun: GseaRunMetadata | null; analysis: ProcessResult | null; onGsea: (result: GseaResult, metadata: GseaRunMetadata) => void; onProgress: (message: string) => void }) {
  const fileInput = useRef<HTMLInputElement>(null);
  const config = initialRun?.config;
  const [gmtPath, setGmtPath] = useState(initialRun?.gmtPath ?? "");
  const [seed, setSeed] = useState(Number(config?.permutation_seed ?? 123));
  const [permutations, setPermutations] = useState(Number(config?.clp_n_permutations ?? 1000));
  const [pCutoff, setPCutoff] = useState(Number(config?.clp_pvalue_cutoff ?? 1));
  const [minSize, setMinSize] = useState(Number(config?.clp_min_gene_set_size ?? 15));
  const [maxSize, setMaxSize] = useState(Number(config?.max_gene_set_size ?? 500));
  const [algorithm, setAlgorithm] = useState(config?.clp_by === "permute" ? "permute" : "multilevel");
  const [advancedOpen, setAdvancedOpen] = useState(false);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");

  async function chooseGmt() {
    const selected = await desktopApi.pickGmtFile();
    if (selected) setGmtPath(selected);
    else if (!desktopApi.isDesktop()) fileInput.current?.click();
  }

  function browserGmt(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    const localPath = file ? window.atop?.files.pathForFile(file) : null;
    if (localPath) setGmtPath(localPath);
    else if (file) setError("Use the Electron desktop shell to grant the analysis runtime access to this GMT file.");
  }

  async function runGsea() {
    if (!analysis) {
      setError("Complete Screen analysis before running GSEA.");
      return;
    }
    if (!gmtPath) {
      setError("Select a GMT gene-set file first.");
      return;
    }
    setBusy(true);
    setError("");
    onProgress("Running clusterProfiler GSEA…");
    try {
      const result = await desktopApi.request<GseaResult>({ command: "gsea", gmtPath, parameters: { clp_n_permutations: permutations, clp_pvalue_cutoff: pCutoff, clp_min_gene_set_size: minSize, max_gene_set_size: maxSize, clp_by: algorithm, clp_seed_logical: true, permutation_seed: seed } });
      onGsea(result, { gmtPath, config: { clp_n_permutations: permutations, clp_pvalue_cutoff: pCutoff, clp_min_gene_set_size: minSize, max_gene_set_size: maxSize, clp_by: algorithm, clp_seed_logical: true, permutation_seed: seed } });
      onProgress("GSEA completed. Pathway results are available.");
    } catch (requestError) {
      setError(requestError instanceof Error ? requestError.message : String(requestError));
      onProgress("GSEA stopped with an error");
    } finally {
      setBusy(false);
    }
  }

  return (
    <div className="workspace-view">
      {!analysis && <Prerequisite>Run screen analysis before enrichment.</Prerequisite>}
      <section className="gsea-form">
        <h2>Gene sets</h2>
        <div className="actions">
          <Button variant="secondary" onClick={() => void chooseGmt()} type="button"><Icon name="file" />Choose local GMT</Button>
          {gmtPath && <span title={gmtPath}>{fileName(gmtPath)}</span>}
          <input accept=".gmt" hidden onChange={browserGmt} ref={fileInput} type="file" />
        </div>
        <Button variant="ghost" className="analysis-disclosure" aria-expanded={advancedOpen}
          onClick={() => setAdvancedOpen(current => !current)} type="button">
          GSEA parameters <Icon name="chevronDown" size={15} />
        </Button>
        {advancedOpen && <div className="field-grid">
          <SelectField label="Algorithm" value={algorithm} options={["multilevel", "permute"]} optionLabels={{multilevel:"Multilevel",permute:"Permutation"}} onChange={setAlgorithm} />
          <NumberField label="Permutation seed" max={2147483647} value={seed} onChange={setSeed}/>
          {algorithm === "permute" && <NumberField label="Permutations" min={1} value={permutations} onChange={setPermutations} />}
          <NumberField label="P-value cutoff" value={pCutoff} step="0.01" onChange={setPCutoff} />
          <NumberField label="Minimum gene-set size" min={1} value={minSize} onChange={setMinSize} />
          <NumberField label="Maximum gene-set size" min={1} value={maxSize} onChange={setMaxSize} />
        </div>}
        <div className="actions run-action"><Button disabled={busy || !analysis || !gmtPath || !Number.isInteger(seed) || seed < 0 || seed > 2147483647 || (algorithm === "permute" && (!Number.isInteger(permutations) || permutations < 1))} onClick={() => void runGsea()} type="button">
          {busy ? <ProcessingIndicator label="Running GSEA…" /> : <><Icon name="play" />Run GSEA</>}
        </Button></div>
      </section>
      {error && <ErrorBanner message={error} />}
    </div>
  );
}

export function GseaResultsView({ gsea }: { gsea: GseaResult | null }) {
  return <div className="workspace-view">{!gsea ? <Prerequisite>Run GSEA to see pathway results.</Prerequisite> : <div className="result-stack"><section className="panel-card"><div className="result-card-heading"><div><h2>Positive pathways</h2></div><span className="count-badge">{gsea.positive.length}</span></div><DataTable rows={gsea.positive} empty="No positive pathways were returned." /></section><section className="panel-card"><div className="result-card-heading"><div><h2>Negative pathways</h2></div><span className="count-badge">{gsea.negative.length}</span></div><DataTable rows={gsea.negative} empty="No negative pathways were returned." /></section></div>}</div>;
}
