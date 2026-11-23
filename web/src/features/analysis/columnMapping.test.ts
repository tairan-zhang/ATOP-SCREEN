import { describe, expect, it } from "vitest";
import { inferMatrixColumnMapping } from "./columnMapping";

describe("inferMatrixColumnMapping", () => {
  it("maps the real ATOP count-matrix header without free-text input", () => {
    expect(inferMatrixColumnMapping(["gRNA", "Gene", "Sequence", "Day0", "Day7C", "Day7T"])).toEqual({
      gRNA: "gRNA",
      gene: "Gene",
      sequence: "Sequence",
      treatment: ["Day7T"],
      control: ["Day7C"],
    });
  });

  it("uses a single T0 for all unlabelled endpoint replicates", () => {
    const sample = inferMatrixColumnMapping(["guide", "gene", "E1", "T0", "E2", "E3"]);
    expect(sample.control).toEqual(["T0"]);
    expect(sample.treatment).toEqual(["E1", "E2", "E3"]);
  });

  it("falls back to distinct count-column groups when names have no condition hints", () => {
    expect(inferMatrixColumnMapping(["sgRNA", "symbol", "replicate_1", "replicate_2"])).toEqual({
      gRNA: "sgRNA",
      gene: "symbol",
      sequence: "",
      treatment: ["replicate_1"],
      control: ["replicate_2"],
    });
  });
});
