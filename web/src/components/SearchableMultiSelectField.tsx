import {useMemo, useState} from "react";
import {Button, Input} from "@atop/ui";
import {Popover, PopoverContent, PopoverTrigger} from "@atop/ui/motion/popover";
type Props = {
  disabled?: boolean;
  emptyMessage?: string;
  hint?: string;
  label: string;
  onChange: (values: string[]) => void;
  options: string[];
  placeholder: string;
  searchPlaceholder?: string;
  values: string[];
};
export function SearchableMultiSelectField({label, options, values, onChange, placeholder, searchPlaceholder, emptyMessage, hint, disabled}: Props) {
  const [open, setOpen] = useState(false);
  const [query, setQuery] = useState("");
  const matches = useMemo(() => options.filter(item => item.toLowerCase().includes(query.trim().toLowerCase())).slice(0,100), [query, options]);
  return <div className="field"><span>{label}</span>
    <Popover open={open} onOpenChange={next => {setOpen(next); if (!next) setQuery("");}} side="bottom" align="start" motionPreset="none">
      <PopoverTrigger><Button variant="outline" type="button" aria-label={label} disabled={disabled || !options.length}>{values.join(", ") || placeholder}</Button></PopoverTrigger>
      <PopoverContent className="screen-picker"><Input autoFocus aria-label={"Search " + label} placeholder={searchPlaceholder ?? "Search…"} value={query} onChange={setQuery}/>
        <div className="screen-option-list">{matches.map(item => <Button key={item} type="button" variant="ghost" aria-pressed={values.includes(item)} onClick={() => {onChange(values.includes(item) ? values.filter(v => v !== item) : [...values,item]);}}>{item}</Button>)}
          {!matches.length && <small>{emptyMessage ?? "No matches"}</small>}</div>
      </PopoverContent>
    </Popover>{hint && <small>{hint}</small>}</div>;
}
