import {useEffect, useState} from "react";
import type {ProcessResult} from "@atop-screen/ipc-contracts";
import {SelectField} from "../../components/SelectField";
import ChartWorker from "./ncs/charts.worker.js?worker&inline";
import labels from "./ncs/labels.json";

const titles = {density:labels.ncs_density, removal:labels.ncs_removal, scatter:labels.ncs_scatter, components:labels.ncs_components, signal:labels.ncs_signal};

export function NcsChartSelector({view, onChange}: {view: string; onChange: (next: string) => void}) {
  return <SelectField label="NCS chart" value={view} onChange={onChange} options={Object.keys(titles)} optionLabels={titles}/>;
}

export function NcsCharts({analysis, view}: {analysis: ProcessResult; view: string}) {
  const [result, setResult] = useState<{charts?:Record<string,string>; error?:string}>({});
  useEffect(() => {
    setResult({});
    const worker = new ChartWorker();
    worker.onmessage = event => setResult(event.data);
    worker.onerror = event => setResult({error:event.message || "Could not draw NCS charts."});
    const {sgrna, ncsExcluded} = analysis.tables;
    worker.postMessage({tables:{sgrna,ncsExcluded}, parameters:analysis.summary.parameters});
    return () => worker.terminate();
  }, [analysis]);
  const svg = result.charts?.[view];
  return <section>
    {result.error ? <p role="alert">{result.error}</p> : svg ? <img className="ncs-chart" alt={titles[view as keyof typeof titles]} src={"data:image/svg+xml;charset=UTF-8," + encodeURIComponent(svg)}/> : <p role="status">Preparing NCS plots…</p>}
  </section>;
}
