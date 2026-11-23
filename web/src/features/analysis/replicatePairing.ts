export type ReplicatePairStatus = "matched" | "review" | "mismatch" | "missing";

export type ReplicatePair = {
  control: string | null;
  index: number;
  status: ReplicatePairStatus;
  treatment: string | null;
};

export function reorderReplicates(values: string[], fromIndex: number, toIndex: number) {
  if (fromIndex === toIndex || fromIndex < 0 || toIndex < 0 || fromIndex >= values.length || toIndex >= values.length) return values;
  const reordered = [...values];
  const [moved] = reordered.splice(fromIndex, 1);
  reordered.splice(toIndex, 0, moved);
  return reordered;
}

export function replicateNumber(value: string) {
  const match = value.trim().match(/(?:rep(?:licate)?[\s._-]*)?(\d+)$/i);
  return match ? Number(match[1]) : null;
}

export function buildReplicatePairs(treatment: string[], control: string[]): ReplicatePair[] {
  const rowCount = Math.max(treatment.length, control.length);
  return Array.from({ length: rowCount }, (_, index) => {
    const treatmentColumn = treatment[index] ?? null;
    const controlColumn = control[index] ?? null;
    if (!treatmentColumn || !controlColumn) return { control: controlColumn, index, status: "missing", treatment: treatmentColumn };
    const treatmentNumber = replicateNumber(treatmentColumn);
    const controlNumber = replicateNumber(controlColumn);
    const status = treatmentNumber === null || controlNumber === null
      ? "review"
      : treatmentNumber === controlNumber ? "matched" : "mismatch";
    return { control: controlColumn, index, status, treatment: treatmentColumn };
  });
}

export function validateReplicatePairing(treatment: string[], control: string[]) {
  const overlap = treatment.filter((column) => control.includes(column));
  return {
    overlap,
    ready: treatment.length > 0 && treatment.length === control.length && overlap.length === 0,
    sameCount: treatment.length === control.length,
  };
}
