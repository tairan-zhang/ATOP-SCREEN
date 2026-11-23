import {SelectField} from "./SelectField";
export function TablePageSizeSelect({value, onChange}: {value:number; onChange:(size:number)=>void}) {
  return <SelectField label="Rows per page" value={String(value)} options={["10","25","50","100"]} onChange={size => onChange(Number(size))}/>;
}
