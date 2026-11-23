import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it, vi } from "vitest";
import { formatValue } from "../domain/format";
import { DataTable } from "./DataTable";

// The virtualized library table needs a measured browser viewport. Exercise
// the page and column contract here; the real table is mounted in Electron.
vi.mock("@atop/ui/motion/table", () => ({
  Table: ({data, columns}: {data:Record<string,string>[]; columns:{key:string}[]}) =>
    <table><tbody>{data.map((row,i) => <tr key={i}>{columns.map(({key}) => <td key={key}>{key} {row[key]}</td>)}</tr>)}</tbody></table>,
}));

describe("DataTable", () => {
  it("distinguishes small nonzero values from zero", () => {
    expect(formatValue(2.3e-7)).toBe("2.30e-7");
    expect(formatValue(0)).toBe("0");
  });

  it("keeps every result column and renders the first paginated page", () => {
    const rows = Array.from({ length: 30 }, (_, rowIndex) => Object.fromEntries(
      Array.from({ length: 12 }, (_, columnIndex) => [`column_${columnIndex + 1}`, `row-${rowIndex + 1}-${columnIndex + 1}`]),
    ));

    const html = renderToStaticMarkup(<DataTable rows={rows} />);

    expect(html).toContain("column_12");
    expect(html).toContain("Rows per page");
    expect(html).not.toContain("Show all");
    expect(html).toContain("1–25 of 30");
    expect(html).toContain("row-25-12");
    expect(html).not.toContain("row-26-1");
  });
});
