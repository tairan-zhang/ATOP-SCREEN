import {useMemo, useState} from "react";
import {Button, Input} from "@atop/ui";
import {Popover, PopoverContent, PopoverTrigger} from "@atop/ui/motion/popover";
type SearchableSelectFieldProps = {
  label: string;
  onChange: (value: string) => void;
  options: string[];
  placeholder: string;
  disabled?: boolean;
  emptyMessage?: string;
  hint?: string;
  searchPlaceholder?: string;
  summaryNoun?: string;
  value: string;
};
export function SearchableSelectField({label, options, value, onChange, placeholder, searchPlaceholder, emptyMessage, hint, disabled}: SearchableSelectFieldProps) {
  const [open, setOpen] = useState(false);
  const [query, setQuery] = useState("");
  const matches = useMemo(() => options.filter(item => item.toLowerCase().includes(query.trim().toLowerCase())).slice(0,100), [query, options]);
  return <div className="field"><span>{label}</span>
    <Popover open={open} onOpenChange={next => {setOpen(next); if (!next) setQuery("");}} side="bottom" align="start" motionPreset="none">
      <PopoverTrigger><Button variant="outline" type="button" aria-label={label} disabled={disabled || !options.length}>{value || placeholder}</Button></PopoverTrigger>
      <PopoverContent className="screen-picker"><Input autoFocus aria-label={"Search " + label} placeholder={searchPlaceholder ?? "Search…"} value={query} onChange={setQuery}/>
        <div className="screen-option-list">{matches.map(item => <Button key={item} type="button" variant="ghost" aria-pressed={value === item} onClick={() => {onChange(item); setOpen(false); setQuery("");}}>{item}</Button>)}
          {!matches.length && <small>{emptyMessage ?? "No matches"}</small>}</div>
      </PopoverContent>
    </Popover>{hint && <small>{hint}</small>}</div>;
}
