import {Button, Input} from "@atop/ui";
import { useRef, useState } from "react";
import type { ChangeEvent, FormEvent } from "react";
import type { FileInspection, ProcessResult } from "@atop-screen/ipc-contracts";
import { desktopApi } from "../../services/desktop-api";
import { fileName, splitColumns } from "../../domain/format";
import { defaultAnalysisValues } from "../../domain/storage";
import type { AnalysisValues } from "../../types/atop";
import { Icon } from "../../components/Icon";
import { ProcessingIndicator } from "../../components/ProcessingIndicator";
import { ReplicatePairingEditor } from "../../components/ReplicatePairingEditor";
import { CheckboxField, MultiSelectField, NumberField, SelectField } from "../../components/FormFields";
import { ErrorBanner } from "../../components/WorkspacePrimitives";
import { inferMatrixColumnMapping } from "./columnMapping";


export type AnalysisRunMetadata = { analysisName: string; config: Record<string, unknown> };

export function AnalysisView({ analysis, inspection, onAnalysis, onInspection, onPath, onProgress, path, projectName }: { analysis: ProcessResult | null; inspection: FileInspection | null; onAnalysis: (result: ProcessResult, metadata: AnalysisRunMetadata) => void; onInspection: (inspection: FileInspection | null) => void; onPath: (path: string) => void; onProgress: (message: string) => void; path: string; projectName: string }) {
  const fileInput = useRef<HTMLInputElement>(null);
  const [values, setValues] = useState<AnalysisValues>(() => {
    if (!inspection) return defaultAnalysisValues;
    const mapping = inferMatrixColumnMapping(inspection.columns);
    const saved = analysis?.summary.parameters;
    const values = {
      ...defaultAnalysisValues,
      gRNA_col: mapping.gRNA, gene_col: mapping.gene, sequence_col: mapping.sequence,
      condition1_replicate_cols: mapping.treatment.join(","), condition2_replicate_cols: mapping.control.join(","),
    };
    if (!saved) return values;
    for (const key of Object.keys(values) as (keyof AnalysisValues)[]) {
      const value = saved[key];
      if (typeof values[key] === "number" && typeof value === "number" && Number.isFinite(value)) {
        Object.assign(values, {[key]: value});
      } else if (typeof values[key] === "string" && (typeof value === "string" || Array.isArray(value))) {
        Object.assign(values, {[key]: Array.isArray(value) ? value.join(",") : value});
      }
    }
    if (saved.ncs_threshold === undefined) values.ncs_threshold = 0;
    return values;
  });
  const [withMageck, setWithMageck] = useState(analysis?.summary.parameters.run_mageck === true);
  const [advancedOpen, setAdvancedOpen] = useState(false);
  const [analysisName, setAnalysisName] = useState("");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");

  function updateValue<K extends keyof AnalysisValues>(key: K, value: AnalysisValues[K]) {
    setValues((current) => ({ ...current, [key]: value }));
  }

  async function inspect(selectedPath: string) {
    if (!selectedPath) return;
    onPath(selectedPath);
    onInspection(null);
    setError("");
    onProgress("Reading columns and previewing the selected file…");
    try {
      const result = await desktopApi.request<FileInspection>({ command: "inspect_file", path: selectedPath });
      onInspection(result);
      const mapping = inferMatrixColumnMapping(result.columns);
      setValues((current) => ({
        ...current,
        gRNA_col: mapping.gRNA,
        gene_col: mapping.gene,
        condition1_replicate_cols: mapping.treatment.join(","),
        condition2_replicate_cols: mapping.control.join(","),
        sequence_col: mapping.sequence,
      }));
      onProgress(`Loaded ${result.fileName} · ${result.columns.length} columns`);
    } catch (requestError) {
      setError(requestError instanceof Error ? requestError.message : String(requestError));
      onProgress("Unable to inspect this file");
    }
  }

  async function chooseFile() {
    const selected = await desktopApi.pickScreenFile();
    if (selected) await inspect(selected);
    else if (!desktopApi.isDesktop()) fileInput.current?.click();
  }

  function handleBrowserFile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.target.files?.[0];
    const localPath = file ? window.atop?.files.pathForFile(file) : null;
    if (localPath) void inspect(localPath);
    else if (file) setError("Use the Electron desktop shell to grant the analysis runtime access to this local file.");
  }

  async function runAnalysis(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    if (!path) {
      setError("Select a screen matrix before starting the analysis.");
      return;
    }
    const treatment = splitColumns(values.condition1_replicate_cols);
    const control = splitColumns(values.condition2_replicate_cols);
    const availableColumns = new Set(inspection?.columns ?? []);
    if (!inspection || !availableColumns.has(values.gRNA_col) || !availableColumns.has(values.gene_col)) {
      setError("Inspect the screen matrix and select its sgRNA and gene columns first.");
      return;
    }
    if (!treatment.length || !control.length || [...treatment, ...control].some((column) => !availableColumns.has(column))) {
      setError("Select at least one treatment and one control column from the inspected matrix.");
      return;
    }
    if (treatment.some((column) => control.includes(column))) {
      setError("Treatment and control selections cannot contain the same count column.");
      return;
    }
    setBusy(true);
    setError("");
    onProgress("Running the canonical ATOP analysis pipeline…");
    try {
      const result = await desktopApi.request<ProcessResult>({
        command: "process_data",
        path,
        parameters: {
          ...values,
          run_mageck: withMageck,
          sequence_col: values.sequence_col || null,
          condition1_replicate_cols: treatment,
          condition2_replicate_cols: control,
          N_perm: Number(values.N_perm),
          pseudo_count_lfc: Number(values.pseudo_count_lfc),
          min_sgrna_threshold: Number(values.min_sgrna_threshold),
        },
      });
      onProgress("Analysis completed. Results and downstream modules are ready.");
      onAnalysis(result, {
        analysisName: analysisName.trim() || projectName,
        config: {
          ...values,
          run_mageck: withMageck,
          condition1_replicate_cols: treatment,
          condition2_replicate_cols: control,
          sequence_col: values.sequence_col || null,
        },
      });
    } catch (requestError) {
      setError(requestError instanceof Error ? requestError.message : String(requestError));
      onProgress("The analysis stopped with an error");
    } finally {
      setBusy(false);
    }
  }

  const columns = inspection?.columns ?? [];
  const countColumns = columns.filter((column) => ![values.gRNA_col, values.gene_col, values.sequence_col].includes(column));
  const treatmentColumns = splitColumns(values.condition1_replicate_cols);
  const controlColumns = splitColumns(values.condition2_replicate_cols);
  let mode = values.comparison_mode;
  if (mode === "auto") mode = treatmentColumns.length === controlColumns.length ? "paired" : controlColumns.length === 1 ? "shared_baseline" : "group_means";
  const paired = mode === "paired";
  const groupMeans = mode === "group_means";
  const mappingReady = Boolean(
    inspection
    && columns.includes(values.gRNA_col)
    && columns.includes(values.gene_col)
    && treatmentColumns.length
    && controlColumns.length
    && treatmentColumns.every((column) => countColumns.includes(column))
    && controlColumns.every((column) => countColumns.includes(column))
    && !treatmentColumns.some(column => controlColumns.includes(column))
    && (!paired || treatmentColumns.length === controlColumns.length)
    && (mode !== "shared_baseline" || controlColumns.length === 1),
  );

  return <div className="workspace-view">
    <section><h2>Count matrix</h2>
      <div className="actions"><Button variant="secondary" disabled={busy} onClick={() => void chooseFile()}><Icon name="upload"/>{inspection ? "Change file" : "Choose data file"}</Button>
        <span className="muted">{inspection ? fileName(path) : "CSV, TSV, TXT or XLSX"}</span></div>
      <input hidden type="file" accept=".csv,.tsv,.txt,.xlsx" onChange={handleBrowserFile} ref={fileInput}/>
    </section>
    {inspection && <form className="analysis-form" onSubmit={event => void runAnalysis(event)}>
      <section><h2>Column mapping</h2><div className="field-grid id-mapping">
        <SelectField label="sgRNA identifier" value={values.gRNA_col} options={columns} onChange={v => updateValue("gRNA_col",v)}/>
        <SelectField label="Gene identifier" value={values.gene_col} options={columns} onChange={v => updateValue("gene_col",v)}/>
        <SelectField label="Guide sequence" value={values.sequence_col || "__none__"} options={["__none__",...columns]} optionLabels={{__none__:"Not supplied"}} onChange={v => updateValue("sequence_col",v === "__none__" ? "" : v)}/>
      </div></section>
      <section><h2>Sample groups</h2><div className="field-grid">
        <SelectField label="LFC comparison" value={values.comparison_mode} onChange={v => updateValue("comparison_mode", v)}
          options={["auto","paired","shared_baseline","group_means"]}
          optionLabels={{auto:"Automatic",paired:"Matched pairs",shared_baseline:"Shared baseline",group_means:"Ratio of group means"}}/>
        <MultiSelectField label="Endpoint / treatment counts" values={treatmentColumns} options={countColumns} onChange={v => updateValue("condition1_replicate_cols",v.join(","))}/>
        <MultiSelectField label="Baseline / control counts" values={controlColumns} options={countColumns} onChange={v => updateValue("condition2_replicate_cols",v.join(","))}/>
      </div>{paired ? <ReplicatePairingEditor control={controlColumns} treatment={treatmentColumns} disabled={busy}
        onTreatmentChange={v => updateValue("condition1_replicate_cols",v.join(","))} onControlChange={v => updateValue("condition2_replicate_cols",v.join(","))}/> : <p>{groupMeans ? "LFC uses group means. NCS requires replicate contrasts." : "Each endpoint is compared with the same baseline."}</p>}
      </section>
      {!groupMeans && <section><h2>N-group Concordance Score</h2>
        
        <div className="field-grid">
          <NumberField label="NCS cutoff" value={values.ncs_threshold} max={1} step="any" onChange={v => updateValue("ncs_threshold",v)}/>
        </div><details className="parameter-details"><summary>NCS scales<Icon name="chevronDown" size={14}/></summary><div className="field-grid">
          <NumberField label="Direction dead zone (δ₀)" value={values.ncs_delta0} step="any" onChange={v => updateValue("ncs_delta0",v)}/>
          <NumberField label="Magnitude scale (δ₁)" value={values.ncs_delta1} step="any" onChange={v => updateValue("ncs_delta1",v)}/>
          <NumberField label="Stability scale (δ₂)" value={values.ncs_delta2} step="any" onChange={v => updateValue("ncs_delta2",v)}/>
        </div></details>
      </section>}
      <section><Input label="Analysis name" placeholder={projectName} maxLength={120} value={analysisName} onChange={setAnalysisName}/>
        <Button className="analysis-disclosure" variant="ghost" type="button" aria-expanded={advancedOpen} onClick={() => setAdvancedOpen(!advancedOpen)}>Permutation settings <Icon name="chevronDown"/></Button>
        {advancedOpen && <div className="field-grid">
          <SelectField label="Permutation engine" options={["cpp","r_parallel"]} optionLabels={{cpp:"C++17",r_parallel:"R parallel"}} value={values.user_engine_choice} onChange={v => updateValue("user_engine_choice",v)}/>
          <NumberField label="Permutation seed" max={2147483647} value={values.permutation_seed} onChange={v => updateValue("permutation_seed",v)}/>
          <NumberField label="Permutations" value={values.N_perm} onChange={v => updateValue("N_perm",v)}/>
          <NumberField label="LFC pseudo-count" value={values.pseudo_count_lfc} step="any" onChange={v => updateValue("pseudo_count_lfc",v)}/>
          <NumberField label="Minimum sgRNAs per gene" min={1} value={values.min_sgrna_threshold} onChange={v => updateValue("min_sgrna_threshold",v)}/>
        </div>}
        <CheckboxField checked={withMageck} label="Run MAGeCK RRA" onChange={setWithMageck}/>
        <div className="actions run-action"><Button type="submit" disabled={busy || !mappingReady || (!groupMeans && (values.ncs_delta1 <= 0 || values.ncs_delta2 <= 0))}>{busy ? <ProcessingIndicator label="Running analysis…"/> : "Start analysis"}</Button></div>
      </section></form>}
    {error && <ErrorBanner message={error}/>}
  </div>;
}
