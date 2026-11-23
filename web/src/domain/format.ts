export function fileName(filePath: string) {
  return filePath.split(/[\\/]/).pop() || filePath;
}

export function fileExtension(filePath: string) {
  const dot = filePath.lastIndexOf(".");
  return dot >= 0 ? filePath.slice(dot).toLowerCase() : "";
}

export function splitColumns(value: string) {
  return value.split(",").map((entry) => entry.trim()).filter(Boolean);
}

export function formatValue(value: unknown) {
  if (value === null || value === undefined || value === "") return "—";
  if (typeof value === "number") {
    if (!Number.isFinite(value)) return "—";
    if (value !== 0 && Math.abs(value) < 0.00001) return value.toExponential(2);
    return value.toLocaleString(undefined, { maximumFractionDigits: 5 });
  }
  if (typeof value === "object") {
    try {
      return JSON.stringify(value);
    } catch {
      return String(value);
    }
  }
  return String(value);
}

export function clamp(value: number, minimum: number, maximum: number) {
  return Math.max(minimum, Math.min(maximum, value));
}
