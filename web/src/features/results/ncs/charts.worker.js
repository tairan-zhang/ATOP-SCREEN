import { prepareModel, snapshot } from "./model.js";
import { renderNcsSvg, views } from "./svg.js";
import labels from "./labels.json";

self.onmessage = ({data: {tables, parameters}}) => {
  try {
    const rows = [...tables.sgrna, ...(tables.ncsExcluded ?? [])].map(row => ({
      guide: row[parameters.gRNA_col], gene: row[parameters.gene_col],
      mean: row.diff_score, ncs: row.NCS, dir: row.dir_score, mag: row.mag_score, stab: row.stab_score,
    }));
    const params = {threshold: parameters.ncs_threshold, topN: 50};
    const model = prepareModel(rows, parameters.min_sgrna_threshold, params);
    const data = snapshot(model, params);
    self.postMessage({charts: Object.fromEntries(views.map(view => [view, renderNcsSvg(data, view, labels)]))});
  } catch (error) {
    self.postMessage({error: String(error)});
  }
};
