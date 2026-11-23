import {Button} from "@atop/ui";
import { useEffect, useMemo, useState } from "react";
import type { FileInspection, GseaResult, ProcessResult, VisualizationResult } from "@atop-screen/ipc-contracts";
import { desktopApi } from "../../services/desktop-api";
import type { Artifact, PlotSection } from "../../types/atop";
import { Icon } from "../../components/Icon";
import { ProcessingIndicator } from "../../components/ProcessingIndicator";
import { CheckboxField, ColorField, NumberField, SelectField, TextField } from "../../components/FormFields";
import { SearchableMultiSelectField } from "../../components/SearchableMultiSelectField";
import { SearchableSelectField } from "../../components/SearchableSelectField";
import { ErrorBanner, Prerequisite } from "../../components/WorkspacePrimitives";

type Props = { analysis: ProcessResult | null; gsea: GseaResult | null; inspection: FileInspection | null; onArtifact: (artifact: Artifact) => void; onProgress: (message: string) => void; onVisualization: (result: VisualizationResult, artifact: Artifact) => void; plotSection: PlotSection };

const enrichmentCurvePalette = ["#C02635", "#00609D", "#D05C68", "#4088B6", "#E0949A", "#80B0CE", "#EFC8CC", "#BFD7E6"];
const plotDimensionsCm: Record<string, readonly [number, number]> = {
  gene_volcano: [6.5, 6.5],
  gene_ranking: [11, 8.5],
  sgrna_paired: [6.5, 6.5],
  gsea_lollipop: [13, 11.5],
  gsea_running: [10.5, 7],
};

function columns(value: unknown) {
  if (Array.isArray(value)) return value.map(String).filter(Boolean);
  return typeof value === "string" ? value.split(",").map((item) => item.trim()).filter(Boolean) : [];
}

function unique(values: string[]) { return [...new Set(values.filter(Boolean))]; }

export function VisualizationView({ analysis, gsea, inspection, onArtifact, onProgress, onVisualization, plotSection }: Props) {
  const [plotType, setPlotType] = useState("gene_volcano");
  const [direction, setDirection] = useState("positive");
  const [topN, setTopN] = useState(10);
  const [pValueCutoff, setPValueCutoff] = useState(0.05);
  const [scoreCutoff, setScoreCutoff] = useState(0.2);
  const [targetGene, setTargetGene] = useState("");
  const [condition1Column, setCondition1Column] = useState("");
  const [condition2Column, setCondition2Column] = useState("");
  const [pathwayIds, setPathwayIds] = useState<string[]>([]);
  const [labelGenes, setLabelGenes] = useState<string[]>([]);
  const [knownGenes, setKnownGenes] = useState<string[]>([]);
  const [highlightedGenes, setHighlightedGenes] = useState<string[]>([]);
  const [pathwayColors, setPathwayColors] = useState<Record<string, string>>({});
  const [subplotType, setSubplotType] = useState("3");
  const [showNes, setShowNes] = useState(true);
  const [showPvalue, setShowPvalue] = useState(true);
  const [pvalueType, setPvalueType] = useState("p.adjust");
  const [legendX, setLegendX] = useState(0.85);
  const [legendY, setLegendY] = useState(0.8);
  const [showRankExplanation, setShowRankExplanation] = useState(true);
  const [rankPositiveLabel, setRankPositiveLabel] = useState("Top-ranked genes (positively correlated)");
  const [rankNegativeLabel, setRankNegativeLabel] = useState("Bottom-ranked genes (negatively correlated)");
  const [volcanoUpColor, setVolcanoUpColor] = useState("#C02635");
  const [volcanoDownColor, setVolcanoDownColor] = useState("#00609D");
  const [volcanoNeutralColor, setVolcanoNeutralColor] = useState("#CECECE");
  const [rankingLowColor, setRankingLowColor] = useState("#EFC8CC");
  const [rankingHighColor, setRankingHighColor] = useState("#C02635");
  const [positivePathwayColor, setPositivePathwayColor] = useState("#EFC8CC");
  const [negativePathwayColor, setNegativePathwayColor] = useState("#BFD7E6");
  const [condition1Color, setCondition1Color] = useState("#00609D");
  const [condition2Color, setCondition2Color] = useState("#C02635");
  const [connectingLineColor, setConnectingLineColor] = useState("#7F7F7F");
  const [widthCm, setWidthCm] = useState(6.5);
  const [heightCm, setHeightCm] = useState(6.5);
  const [fontSizePt, setFontSizePt] = useState(8);
  const [pngDpi, setPngDpi] = useState(300);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState("");
  const parameters = analysis?.summary.parameters as Record<string, unknown> | undefined;
  const geneColumn = typeof parameters?.gene_col === "string" ? parameters.gene_col : "Gene";
  const geneOptions = useMemo(() => unique((analysis?.tables.genes ?? []).map((row) => String(row[geneColumn] ?? ""))), [analysis, geneColumn]);
  const countColumns = useMemo(() => {
    const configured = unique([...columns(parameters?.condition1_replicate_cols), ...columns(parameters?.condition2_replicate_cols)]);
    if (configured.length) return configured;
    const identity = new Set([parameters?.gRNA_col, parameters?.gene_col, parameters?.sequence_col].filter((value): value is string => typeof value === "string"));
    return (inspection?.columns ?? []).filter((column) => !identity.has(column));
  }, [inspection, parameters]);
  const pathwayOptions = useMemo(() => unique((direction === "negative" ? gsea?.negative : gsea?.positive)?.map((row) => String(row.ID ?? row.Description ?? "")).filter(Boolean) ?? []), [direction, gsea]);
  const isSgrna = plotSection === "sgrna";
  const isPathway = plotSection === "pathway";

  useEffect(() => { setPlotType(isSgrna ? "sgrna_paired" : isPathway ? "gsea_lollipop" : "gene_volcano"); }, [isPathway, isSgrna, plotSection]);
  useEffect(() => {
    const [defaultWidth, defaultHeight] = plotDimensionsCm[plotType] ?? [18, 14];
    setWidthCm(defaultWidth);
    setHeightCm(defaultHeight);
  }, [plotType]);
  useEffect(() => {
    setCondition1Column((current) => countColumns.includes(current) ? current : countColumns[0] || "");
    setCondition2Column((current) => countColumns.includes(current) && current !== countColumns[0] ? current : countColumns[1] || countColumns[0] || "");
  }, [countColumns]);
  useEffect(() => { setPathwayIds((current) => current.filter((pathway) => pathwayOptions.includes(pathway))); }, [pathwayOptions]);
  useEffect(() => {
    setPathwayColors((current) => {
      const next: Record<string, string> = {};
      pathwayIds.forEach((pathway, index) => { next[pathway] = current[pathway] ?? enrichmentCurvePalette[index % enrichmentCurvePalette.length]; });
      return next;
    });
  }, [pathwayIds]);

  const changeDirection = (nextDirection: string) => {
    setDirection(nextDirection);
    setPathwayIds([]);
    setPathwayColors({});
    if (plotType === "gene_ranking") {
      if (nextDirection === "negative") {
        setRankingLowColor("#BFD7E6"); setRankingHighColor("#00609D");
      } else {
        setRankingLowColor("#EFC8CC"); setRankingHighColor("#C02635");
      }
    }
  };

  async function createVisualization() {
    if (!analysis) return setError("Complete Screen analysis before creating a visualization.");
    if (plotType.startsWith("gsea_") && !gsea) return setError("Run GSEA before creating a pathway visualization.");
    if (isSgrna && (!targetGene || !condition1Column || !condition2Column || condition1Column === condition2Column)) return setError("Select a gene and two different count columns.");
    if (plotType === "gsea_running" && pathwayIds.length === 0) return setError("Select at least one pathway for the enrichment curve.");
    setBusy(true); setError(""); onProgress("Rendering the selected figure…");
    try {
      const result = await desktopApi.request<VisualizationResult>({ command: "visualize", parameters: {
        plotType, direction, topN, pValueCutoff, scoreCutoff,
        targetGene: targetGene || null,
        condition1Column: condition1Column || null,
        condition2Column: condition2Column || null,
        pathwayIds,
        pathwayColors: pathwayIds.map((pathway, index) => pathwayColors[pathway] ?? enrichmentCurvePalette[index % enrichmentCurvePalette.length]),
        subplotType: Number(subplotType), showNes, showPvalue, pvalueType, legendX, legendY,
        showRankExplanation, rankPositiveLabel, rankNegativeLabel, gseaTermWidth: 20,
        labelGenes, knownGenes, highlightedGenes,
        colUp: volcanoUpColor, colDown: volcanoDownColor, colNs: volcanoNeutralColor,
        gradientLow: rankingLowColor, gradientHigh: rankingHighColor,
        positiveColor: positivePathwayColor, negativeColor: negativePathwayColor,
        condition1Color, condition2Color, connectingLineColor,
        title: null, widthCm, heightCm, fontSizePt, pngDpi,
      } });
      const pathwayScope = pathwayIds.length ? pathwayIds.slice(0, 2).join("-") : direction;
      const scope = targetGene ? `gene-${targetGene}` : isPathway ? `pathway-${pathwayScope}` : `genes-${direction}`;
      const artifact = { id: `visualization-${Date.now()}`, title: result.title, downloadLabel: `${plotType}_${scope}`, paths: result.assets };
      onVisualization(result, artifact); onArtifact(artifact); onProgress("Figure rendered. PNG preview is ready.");
    } catch (requestError) { const message = requestError instanceof Error ? requestError.message : String(requestError); setError(message); onProgress("Figure rendering stopped with an error"); }
    finally { setBusy(false); }
  }

  const plotTypes = isSgrna ? ["sgrna_paired"] : isPathway ? ["gsea_lollipop", "gsea_running"] : ["gene_volcano", "gene_ranking"];
  const disabled = busy || !analysis || (plotType.startsWith("gsea_") && !gsea) || (isSgrna && (!targetGene || !condition1Column || !condition2Column || condition1Column === condition2Column)) || (plotType === "gsea_running" && pathwayIds.length === 0);
  return <div className="workspace-view visualization-view">
    {!analysis && <Prerequisite>Run screen analysis before creating a figure.</Prerequisite>}
    {analysis && isPathway && !gsea && <Prerequisite>Run GSEA to plot pathways.</Prerequisite>}
    <section className="panel-card visualization-panel">
      <div className="visualization-settings-grid">
        <SelectField label="Figure type" onChange={setPlotType} optionLabels={{ gene_volcano: "Gene volcano", gene_ranking: "Top gene score barplot", gsea_lollipop: "GSEA pathway ranking", gsea_running: "GSEA enrichment curve", sgrna_paired: "sgRNA paired" }} options={plotTypes} value={plotType} />
        {isSgrna ? <>
          <SearchableSelectField disabled={!analysis} label="Target gene" onChange={setTargetGene} options={geneOptions} placeholder="Select a gene" value={targetGene} />
          <SelectField disabled={!analysis} label="First group" onChange={setCondition1Column} options={countColumns} placeholder="Select a count column" value={condition1Column} />
          <SelectField disabled={!analysis} label="Second group" onChange={setCondition2Column} options={countColumns} placeholder="Select a different count column" value={condition2Column} />
        </> : null}
        {plotType === "gene_volcano" ? <>
          <NumberField label="P-value cutoff" onChange={setPValueCutoff} step="0.01" value={pValueCutoff} />
          <NumberField label="Score cutoff" onChange={setScoreCutoff} step="0.05" value={scoreCutoff} />
          <SearchableMultiSelectField label="Gene labels" onChange={setLabelGenes} options={geneOptions} placeholder="Select genes to label" searchPlaceholder="Search gene symbols…" values={labelGenes} />
        </> : null}
        {plotType === "gene_ranking" ? <>
          <SelectField label="Direction" onChange={changeDirection} optionLabels={{ positive: "Positive score", negative: "Negative score" }} options={["positive", "negative"]} value={direction} />
          <NumberField label="Number of genes" onChange={setTopN} value={topN} />
          <SearchableMultiSelectField label="Genes to mark (*)" onChange={setKnownGenes} options={geneOptions} placeholder="Select genes to mark" searchPlaceholder="Search gene symbols…" values={knownGenes} />
        </> : null}
        {plotType === "gsea_lollipop" ? <>
          <SelectField label="Direction" onChange={changeDirection} optionLabels={{ positive: "Positive pathways", negative: "Negative pathways" }} options={["positive", "negative"]} value={direction} />
          <NumberField label="Number of pathways" onChange={setTopN} value={topN} />
        </> : null}
        {plotType === "gsea_running" ? <>
          <SelectField label="Direction" onChange={changeDirection} optionLabels={{ positive: "Positive pathways", negative: "Negative pathways" }} options={["positive", "negative"]} value={direction} />
          <SearchableMultiSelectField disabled={!gsea} emptyMessage="No matching pathway" label="Pathways" onChange={setPathwayIds} options={pathwayOptions} placeholder="Select pathways" searchPlaceholder="Search pathways…" values={pathwayIds} />
          <SearchableMultiSelectField disabled={!gsea} emptyMessage="No matching gene" label="Highlight genes" onChange={setHighlightedGenes} options={geneOptions} placeholder="Select genes to highlight" searchPlaceholder="Search gene symbols…" values={highlightedGenes} />
          <SelectField label="Plot content" onChange={setSubplotType} optionLabels={{ "1": "Enrichment score", "2": "Score + pathway hits", "3": "Score + hits + ranked list" }} options={["1", "2", "3"]} value={subplotType} />
        </> : null}
      </div>
      <details className="visualization-disclosure" open>
        <summary><span><strong>Colors</strong></span><Icon name="chevronDown" size={15} /></summary>
        <div className="visualization-appearance-grid">
          {plotType === "gene_volcano" ? <>
            <ColorField label="Enriched genes" onChange={setVolcanoUpColor} value={volcanoUpColor} />
            <ColorField label="Depleted genes" onChange={setVolcanoDownColor} value={volcanoDownColor} />
            <ColorField label="Neutral genes" onChange={setVolcanoNeutralColor} value={volcanoNeutralColor} />
          </> : null}
          {plotType === "gene_ranking" ? <>
            <ColorField label="Shorter-bar color" onChange={setRankingLowColor} value={rankingLowColor} />
            <ColorField label="Longer-bar color" onChange={setRankingHighColor} value={rankingHighColor} />
          </> : null}
          {plotType === "gsea_lollipop" ? <>
            <ColorField label="Positive pathway color" onChange={setPositivePathwayColor} value={positivePathwayColor} />
            <ColorField label="Negative pathway color" onChange={setNegativePathwayColor} value={negativePathwayColor} />
          </> : null}
          {plotType === "gsea_running" ? pathwayIds.length ? pathwayIds.map((pathway, index) => <ColorField key={pathway} label={pathway} onChange={(color) => setPathwayColors((current) => ({ ...current, [pathway]: color }))} value={pathwayColors[pathway] ?? enrichmentCurvePalette[index % enrichmentCurvePalette.length]} />) : <p className="visualization-color-empty">Select a pathway to configure its curve color.</p> : null}
          {plotType === "sgrna_paired" ? <>
            <ColorField label="Condition 1 point color" onChange={setCondition1Color} value={condition1Color} />
            <ColorField label="Condition 2 point color" onChange={setCondition2Color} value={condition2Color} />
            <ColorField label="Connecting line color" onChange={setConnectingLineColor} value={connectingLineColor} />
          </> : null}
        </div>
      </details>
      {plotType === "gsea_running" ? <details className="visualization-disclosure" open>
        <summary><span><strong>Enrichment annotations</strong></span><Icon name="chevronDown" size={15} /></summary>
        <div className="visualization-appearance-grid">
          <CheckboxField checked={showNes} label="Show NES" onChange={setShowNes} />
          <CheckboxField checked={showPvalue} label="Show P value" onChange={setShowPvalue} />
          <SelectField disabled={!showPvalue} label="P-value type" onChange={setPvalueType} optionLabels={{ "p.adjust": "Adjusted P", pvalue: "Raw P" }} options={["p.adjust", "pvalue"]} value={pvalueType} />
          {pathwayIds.length > 1 ? <>
            <NumberField label="Legend X (0–1)" onChange={setLegendX} step="0.05" value={legendX} />
            <NumberField label="Legend Y (0–1)" onChange={setLegendY} step="0.05" value={legendY} />
          </> : null}
          {subplotType === "3" ? <>
            <CheckboxField checked={showRankExplanation} label="Explain ranked list" onChange={setShowRankExplanation} />
            <TextField label="Positive-rank label" onChange={setRankPositiveLabel} value={rankPositiveLabel} />
            <TextField label="Negative-rank label" onChange={setRankNegativeLabel} value={rankNegativeLabel} />
          </> : null}
        </div>
      </details> : null}
      {isSgrna && !countColumns.length ? <div className="inline-notice"><Icon name="info" size={14} /> No count columns were restored. Reopen this Screen result, then try again.</div> : null}
    </section>
    <section className="panel-card visualization-panel visualization-output-panel"><header className="visualization-panel__header"><div><h2>Output</h2></div></header><div className="visualization-output-grid"><NumberField label="Width (cm)" onChange={setWidthCm} step="0.1" value={widthCm} /><NumberField label="Height (cm)" onChange={setHeightCm} step="0.1" value={heightCm} /><NumberField label="Font size (pt)" onChange={setFontSizePt} step="0.5" value={fontSizePt} /><NumberField label="PNG resolution (dpi)" onChange={setPngDpi} value={pngDpi} /></div><div className="reference-primary-action"><Button className="primary-button primary-button-large" disabled={disabled} onClick={() => void createVisualization()} type="button">{busy ? <ProcessingIndicator label="Rendering figure…" /> : <><Icon name="image" />Render figure</>}</Button></div></section>
    {error ? <ErrorBanner message={error} /> : null}
  </div>;
}
