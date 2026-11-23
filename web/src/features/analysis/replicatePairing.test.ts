import { describe, expect, it } from "vitest";
import { buildReplicatePairs, reorderReplicates, validateReplicatePairing } from "./replicatePairing";

describe("replicate pairing", () => {
  it("keeps the exact user-defined order", () => {
    expect(reorderReplicates(["T1", "T2", "T3"], 2, 0)).toEqual(["T3", "T1", "T2"]);
    expect(reorderReplicates(["C1", "C2", "C3"], 0, 2)).toEqual(["C2", "C3", "C1"]);
  });

  it("reports matched, mismatched, and missing replicate rows", () => {
    expect(buildReplicatePairs(["Day7T_1", "Day7T_2", "Day7T_3"], ["Day7C_1", "Day7C_3"]).map((pair) => pair.status))
      .toEqual(["matched", "mismatch", "missing"]);
  });

  it("requires equal, non-overlapping replicate groups", () => {
    expect(validateReplicatePairing(["T1", "T2"], ["C1", "C2"]).ready).toBe(true);
    expect(validateReplicatePairing(["T1"], ["C1", "C2"]).ready).toBe(false);
    expect(validateReplicatePairing(["T1"], ["T1"]).ready).toBe(false);
  });
});
