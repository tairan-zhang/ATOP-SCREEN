export type MatrixColumnMapping = {
  gRNA: string;
  gene: string;
  sequence: string;
  treatment: string[];
  control: string[];
};

function normalizedName(value: string) {
  return value.toLowerCase().replace(/[^a-z0-9]/g, "");
}

function findColumn(columns: string[], aliases: string[], excluded: string[] = []) {
  const blocked = new Set(excluded);
  const normalizedAliases = new Set(aliases.map(normalizedName));
  return columns.find((column) => !blocked.has(column) && normalizedAliases.has(normalizedName(column))) ?? "";
}

export function inferMatrixColumnMapping(columns: string[]): MatrixColumnMapping {
  const gRNA = findColumn(columns, ["gRNA", "sgRNA", "guide", "guide_id"]) || columns[0] || "";
  const gene = findColumn(columns, ["Gene", "gene_id", "symbol", "gene_symbol"], [gRNA])
    || columns.find((column) => column !== gRNA)
    || "";
  const sequence = findColumn(columns, ["Sequence", "seq", "sgRNA_sequence", "guide_sequence"], [gRNA, gene]);
  const identifiers = new Set([gRNA, gene, sequence].filter(Boolean));
  const countColumns = columns.filter((column) => !identifiers.has(column));
  const controls = countColumns.filter(column => /ctrl|control|vehicle|untreat|day\d+c(?:_|$)/i.test(column));
  const baselines = countColumns.filter(column => /baseline|initial|(?:^|[_ .-])(?:t0|day0)(?:$|[_ .-])/i.test(column));
  const control = controls.length ? controls : baselines;
  const treatment = countColumns.filter(column => !control.includes(column) && /treat|endpoint|final|case|day\d+t(?:_|$)/i.test(column));
  const unassigned = countColumns.filter((column) => !treatment.includes(column) && !control.includes(column) && !baselines.includes(column));

  return {
    gRNA,
    gene,
    sequence,
    treatment: treatment.length ? treatment : control.length ? unassigned : unassigned.slice(0, Math.max(1, Math.ceil(unassigned.length / 2))),
    control: control.length ? control : treatment.length ? unassigned : unassigned.slice(Math.max(1, Math.ceil(unassigned.length / 2))),
  };
}
