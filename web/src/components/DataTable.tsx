import {Table} from "@atop/ui/motion/table";
import {Button} from "@atop/ui";
import { useEffect, useMemo, useState, type CSSProperties } from "react";
import { formatValue } from "../domain/format";
import { Icon } from "./Icon";
import { TablePageSizeSelect } from "./TablePageSizeSelect";

export function DataTable({ rows, empty = "No result rows yet." }: { rows: ReadonlyArray<Record<string, unknown>> | null | undefined; empty?: string }) {
  const safeRows = rows ?? [];
  const columns = useMemo(() => {
    const seen = new Set<string>();
    for (const row of safeRows) {
      for (const column of Object.keys(row)) seen.add(column);
    }
    return Array.from(seen);
  }, [safeRows]);
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(25);
  const pageCount = Math.max(1, Math.ceil(safeRows.length / pageSize));

  useEffect(() => {
    setPage(1);
  }, [rows]);

  useEffect(() => {
    if (page > pageCount) setPage(pageCount);
  }, [page, pageCount]);

  if (!columns.length) return <div className="empty-table"><Icon name="database" size={18} /><span>{empty}</span></div>;
  const visibleRows = safeRows.slice((page - 1) * pageSize, page * pageSize);
  const firstRow = (page - 1) * pageSize + 1;
  const lastRow = Math.min(page * pageSize, safeRows.length);

  return (
    <div className="data-table-shell" style={{"--table-content-height": (visibleRows.length + 1) * 32 + 2 + "px"} as CSSProperties}>
      <div className="data-table-wrap">
        <Table<Record<string, unknown>> data={visibleRows} rowHeight={32} height={Math.min(400, visibleRows.length * 32 + 62)} overscan={pageSize}
          columns={columns.map(key => ({key, header:key, width:Math.max(110, key.length * 8 + 48) + "px", cell:row => <span title={typeof row[key] === "number" ? String(row[key]) : formatValue(row[key])}>{formatValue(row[key])}</span>}))}/>
      </div>
      {safeRows.length > 10 && <div className="data-table-pagination">
        <TablePageSizeSelect value={pageSize} onChange={(size) => { setPageSize(size); setPage(1); }} />
        {pageCount > 1 && (
          <div className="data-table-pagination__controls">
            <Button variant="secondary" aria-label="Previous page" disabled={page === 1} onClick={() => setPage((current) => Math.max(1, current - 1))} type="button"><Icon name="arrowLeft" size={14} /></Button>
            <span>{firstRow.toLocaleString()}–{lastRow.toLocaleString()} of {safeRows.length.toLocaleString()}</span>
            <Button variant="secondary" aria-label="Next page" disabled={page === pageCount} onClick={() => setPage((current) => Math.min(pageCount, current + 1))} type="button"><Icon name="arrowRight" size={14} /></Button>
          </div>
        )}
      </div>}
    </div>
  );
}
