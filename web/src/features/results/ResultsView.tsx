import { Button } from "@atop/ui";
import type { ProcessResult } from "@atop-screen/ipc-contracts";
import { DataTable } from "../../components/DataTable";
import { Icon } from "../../components/Icon";
import { Prerequisite } from "../../components/WorkspacePrimitives";
import { NcsCharts, NcsChartSelector } from "./NcsCharts";

type Output = { id: string; label: string; group: string; rows?: ProcessResult["tables"]["genes"] };
type ResultProps = { analysis: ProcessResult | null; tab: string; chart: string };

function outputs(analysis: ProcessResult): Output[] {
  const { tables, summary: { parameters } } = analysis;
  const items: Output[] = [
    { id: "genes", label: "Gene summary", group: "ATOP", rows: tables.genes },
    { id: "sgrna", label: "sgRNA data", group: "ATOP", rows: tables.sgrna },
    { id: "filtered", label: "Below minimum guides", group: "ATOP", rows: tables.filtered ?? [] },
  ];
  if (parameters.ncs_threshold !== undefined && parameters.ncs_enabled !== false && !parameters.skip_normalization)
    items.push({ id: "ncsCharts", label: "NCS plots", group: "Quality control" });
  if (tables.ncsQc) items.push({ id: "ncs", label: "NCS retention by gene", group: "Quality control", rows: tables.ncsQc });
  if (tables.ncsExcluded?.length) items.push({ id: "ncsExcluded", label: "Below NCS cutoff", group: "Quality control", rows: tables.ncsExcluded });
  if (tables.mageckGenes) items.push(
    { id: "mageckGenes", label: "MAGeCK genes", group: "MAGeCK RRA", rows: tables.mageckGenes },
    { id: "mageckGuides", label: "MAGeCK sgRNAs", group: "MAGeCK RRA", rows: tables.mageckGuides ?? [] },
  );
  return items;
}

export function ResultPicker({ analysis, tab, chart, onTab, onChart }: ResultProps & {
  onTab: (id: string) => void; onChart: (id: string) => void;
}) {
  if (!analysis) return null;
  const items = outputs(analysis);
  const selected = items.find(item => item.id === tab) ?? items[0];
  return <div className="workspace-view result-picker">
    {Array.from(new Set(items.map(item => item.group))).map(group =>
      <section key={group}>
        <h2>{group}</h2>
        <div className="result-options">
          {items.filter(item => item.group === group).map(item =>
            <Button key={item.id} type="button" variant="ghost" aria-pressed={item.id === selected.id}
              className="output-option" onClick={() => onTab(item.id)}>
              <Icon name={item.id === "ncsCharts" ? "bar" : "database"}/><span>{item.label}</span>
              {item.rows && <small>{item.rows.length.toLocaleString()}</small>}
            </Button>
          )}
        </div>
        {selected.id === "ncsCharts" && group === selected.group && <NcsChartSelector view={chart} onChange={onChart}/>}
      </section>
    )}
  </div>;
}

export function ResultsView({ analysis, tab, chart }: ResultProps) {
  if (!analysis) return <Prerequisite>Run screen analysis to see results.</Prerequisite>;
  const items = outputs(analysis);
  const active = items.find(item => item.id === tab) ?? items[0];
  return <section className="results-table-panel">
    <header className="canvas-heading"><Icon name={active.id === "ncsCharts" ? "bar" : "database"}/><h2>{active.label}</h2></header>
    {active.id === "ncsCharts" ? <div className="ncs-sheet"><NcsCharts analysis={analysis} view={chart}/></div> :
      <DataTable key={active.id} rows={active.rows} empty={active.id === "filtered" ? "All genes met the minimum guide count." : "No rows."}/>}
  </section>;
}
